"""
Copasi Process
"""
from process_bigraph import Process, Composite, ProcessTypes
from basico import load_model, get_species, run_time_course
#from utils.basico_helper import _set_initial_concentrations, _get_transient_concentration
import matplotlib.pyplot as plt

import sys 
import os
import COPASI
import argparse
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)
import json



toy_model="BindingKa.cps"

#model="C:/Users/mabda/Desktop/VSC R3/multiscale/BindingKa.cps"


# this is how we set values in COPASI really fast
def _set_initial_concentrations(changes,dm):
    model = dm.getModel()
    assert(isinstance(model, COPASI.CModel))

    references = COPASI.ObjectStdVector()

    for name, value in changes.items():
        species = model.getMetabolite(name)
        # assert(isinstance(species, COPASI.CMetab))
        if species is None:
            print(f"Species {name} not found in model")
            continue
        species.setInitialConcentration(value)
        references.append(species.getInitialConcentrationReference())
#Apply all updates in one go
    model.updateInitialValues(references)

# this is how we pull values out of COPASI fast
def _get_transient_concentration(name, dm): #gets the current [] of a species from model after simulation
    model = dm.getModel()
    assert(isinstance(model, COPASI.CModel))
    species = model.getMetabolite(name)
    assert(isinstance(species, COPASI.CMetab))
    if species is None:
        print(f"Species {name} not found in model")
        return None
    return species.getConcentration()


#Custom Vivarium process - loads COPASI model, takes input [], runs a simulation, returns updated []
class CopasiModule(Process):

    defaults = {
        'model_file': 'string',
        # 'time_step': 1.0,
        'linkers' : 'list'
    }

    def initialize(self, parameters=None, core=None):
        self.copasi_model_object = load_model(self.config['model_file'])
        self.all_species = get_species(model=self.copasi_model_object).index.tolist() #saves all species names  
        self.ic_default = get_species(model=self.copasi_model_object)["initial_concentration"].values #gets their inital concentrations
        self.linker_species = self.config['linkers'] #what ever is passed as update 
        
        #Inputs and outputs --> tells vivarium how to pass data in and get data out (like ports)
    
    def inputs(self):
        return {
            'species': 'map[float]'  # TODO -- make this a narrower type
        }
        
    def outputs(self):
        return {
            'species':
                'map[float]'  # TODO -- make this a narrower type
        }
        


    def update(self, states, interval):

        # get the current linker species levels
        species_levels = states['species']
        #TODO: check that species level is {ME:value}

        # set this value within copasi
        _set_initial_concentrations(species_levels,self.copasi_model_object)

        # run a short copasi simulation
        timecourse = run_time_course(duration=interval, intervals=1, update_model=True, model=self.copasi_model_object)

        # get the values out
        results = { 
                   mol_id:_get_transient_concentration(name=mol_id,dm=self.copasi_model_object)
                   for mol_id in self.all_species}
        #TODO: check that results looks like {ME:Delta value}
        # return the values through the species port
        return {
            'species': results}
        

#class TransportProcess(Process):
    # defaults = 
    

        
def run_copasi_process(core): #set up and run the vivarium simulation
    print('this ran!')
    total_time = 10
    json_document_state = create_vivarium_file(model=toy_model, linkers=['a'], number_cells=2) 
    

    # creat the vivarium simulation from document
    sim = Composite({'state': json_document_state}, core=core) #,core=core
    
   # initial_state = {}
    # run the simulation
    sim.run(total_time)


def create_vivarium_file(
    model, 
    linkers=None, 
    output_json=None,
    transport_data=None,
    number_cells=None,
    ): 
    #this function creates and saves a vivarium file according to the specs of newfilename and output_json
    cell_linkers = [] 
    for i in range(0, number_cells): 
        for j in linkers: 
            linker_name = f"{i}_{j}"
            cell_linkers.append(linker_name)
    json_document_state = {
        'external_species': {
            linker: 1 for linker in cell_linkers
            },
        'cells': {
            f'cell{i}': {
                '_type': 'process',
                'address': 'local:copasi',
                'config': {
                    'model_file': model,
                    'linkers': linkers
                },
                'inputs': {
                    'species': ['..', 'external_species']
                },
                'outputs': {
                    'species': ['..', 'external_species']
                }
            }
             for i in range(0,number_cells)
        }
        # 'transport': {
        #     '_type': 'process',
        #     'address': 'local:transport',  #TODO: need to make this
        #     'config': {
            #   ########
            # },
        #     'inputs': {
                    #'species': 'external species'
            # },
        #     'outputs': {
            #       # 'species': 'external species'
            # },
        # }
    }
    
    #TODO: add initial state from copasi file
    with open(output_json, "w") as outfile:
        json.dump(json_document_state, outfile, indent=2)
        print(f'Vivarium model written to {output_json}')
    
    return json_document_state

def register_copasi_types(core):
    core.register_process('copasi', CopasiModule)


if __name__ == '__main__':
    core=ProcessTypes()
    register_copasi_types(core=core)
    run_copasi_process(core=core)
    
#import json 
#with open('copasi_state.json', 'w') as f: 
 #   json.dump(json_document_state, f, indent=2)
    
