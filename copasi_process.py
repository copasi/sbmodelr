"""
Copasi Process
"""

from process_bigraph import Process, Composite, ProcessTypes
from basico.model_io import load_model, save_model
from basico.model_info import get_species
from basico.task_timecourse import run_time_course
import matplotlib.pyplot as plt
import sys
import os
import COPASI
import argparse
import warnings
import json
from basico.model_info import add_species, add_reaction
from basico.model_info import add_function

warnings.filterwarnings("ignore", category=SyntaxWarning)

toy_model="BindingKa.cps"


# How we set values in COPASI really fast
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
# Apply all updates in one go
    model.updateInitialValues(references)

# How we pull values out of COPASI fast
def _get_transient_concentration(name, dm):
    model = dm.getModel()
    assert(isinstance(model, COPASI.CModel))
    species = model.getMetabolite(name)
    assert(isinstance(species, COPASI.CMetab))
    if species is None:
        print(f"Species {name} not found in model")
        return None
    return species.getConcentration()


# Custom Vivarium process
class CopasiModule(Process):

    defaults = {
        'model_file': 'string',
        'linkers' : 'list',
        'transport_rate': 'float',
        'cell_id': 'float'
    }

    def initialize(self, config):
        self.copasi_model_object = load_model(self.config['model_file']) 
        self.all_species = get_species(model=self.copasi_model_object).index.tolist()
        self.ic_default = get_species(model=self.copasi_model_object)["initial_concentration"].values
        self.linker_species = self.config['linkers']
        self.transport_rate = self.config['transport_rate']
        self.external_linkers = []
        self.cid = self.config['cell_id']
        
        # External species that does the transport
        for linker in self.linker_species:
            ext = f'{linker}_external'
            if ext not in get_species(model=self.copasi_model_object).index:
                add_species(name = ext,
                    initial_concentration=1.0,
                    model=self.copasi_model_object
                    )
                self.external_linkers.append(ext)
                add_reaction(
                    model=self.copasi_model_object,
                    name=f'{linker}exchange',
                    function='Mass action (reversible)',
                    mapping={
                        'substrate': linker,
                        'product': ext,
                        'k1': self.transport_rate,
                        'k2':'k1',
                    },
                )
        # Save model after external species and reaction is added 
        base,_=os.path.splitext(self.config['model_file'])
        outfile=base + f"{'cell_id'}.cps"
        save_model(outfile,
                    model=self.copasi_model_object,
                    type='copasi',
                    overwrite=True)
        print(f"New COPASI model saved with external species to {outfile}")
        self.config['model_file']=outfile

             
    def inputs(self):
        return {
            'species': 'map[float]'
        }
        
    def outputs(self):
        return {
            'species':
                'map[float]'
        }
        

    def update(self, states, interval): 
        external_set = {
            ext: states['species'][ext]
            for ext in self.external_linkers
        }
        
        _set_initial_concentrations(external_set,
                                    self.copasi_model_object)
        
        timecourse = run_time_course(duration=interval, 
                                     intervals=1,
                                     update_model=True,
                                     model=self.copasi_model_object)
        results = { 
                   mol_id:_get_transient_concentration(name=mol_id,dm=self.copasi_model_object)
                   for mol_id in self.all_species}
        return {
            'species': results}
        

# Sets up and runs the vivarium simulation
def run_copasi_process(core): 
    print('this ran!')
    total_time = 10
    json_document_state = create_vivarium_file(
        model=toy_model,
        linkers=['a'],
        edges=[('cell 0', 'cell 1')],
        )
    
    # Create the vivarium simulation from document
    sim = Composite({'state': json_document_state}, core=core)
    sim.run(total_time)


# Creating JSON function
def create_vivarium_file(
    model, 
    linkers=None, 
    output_json=None,
    transport_data=None,
    cell_ids=None,
    transport_rate=None,
    edges=None
    ): 
    edges = edges or [] 
    
    
    if cell_ids is None: 
        cell_ids = []
    for edge in edges:
        c1, c2 = edge
        if c1 not in cell_ids: 
            cell_ids.append(c1)
        if c2 not in cell_ids:
            cell_ids.append(c2)
  
    
    external_species = {
        cid: { f"{linker}_external": 1.0 for linker in linkers }
        for cid in cell_ids
    }    
    
    # Make the JSON document
    json_document_state = {
        'external_species': external_species,
        'cells': {
            cid: {
                '_type': 'process',
                'address': 'local:copasi',
                'config': {
                    'model_file': model,
                    'linkers': linkers,
                    'cell_id': cid,
                },
                'inputs': {
                    'species': ['..', 'external_species', cid]
                },
                'outputs': {
                    'species': ['..', 'external_species', cid]
                }
            }
            for cid in cell_ids    
        },
    }
    
    if output_json:
        with open(output_json, "w") as outfile:
            json.dump(json_document_state, outfile, indent=2)
            print(f'Vivarium model written to {output_json}')

    return json_document_state

def register_copasi_types(core):
    core.register_process('copasi', CopasiModule),


if __name__ == '__main__':
    core=ProcessTypes()
    register_copasi_types(core=core)
    run_copasi_process(core=core)
