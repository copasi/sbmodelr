#!/usr/bin/env python3
#
# sbmodelr
# Model Replicator
#
# This takes a SBML or COPASI file and replicates it as a set of sub-models
# which can exist just side-by-side or can be connected in different ways.
#
# This code is released under the Artistic 2.0 License
# Initially written March-August 2024 by Pedro Mendes <pmendes@uchc.edu>
# Additions by Maya Abdalla January-May 2025

# v. 1.0: baseline functionality
# v. 1.1: bugfix allowing models to have species in more than one
#         compartment with same name
# v. 1.2: rationalize all topologies through network representation

#TODO: v.1.3: incorporate Maya's additions (vivarium mode)
#TODO: v.1.4: gene regulation by adding corresponding protein?
#TODO: WISH LIST
#TODO: binding between units?  (would nedd to select A,B and bind each one on each side)
#TODO: modification between units? (would need to select species and reaction)
#TODO: hexagonal arrays? (please no!)

__version__ = "1.3"

import os
import re
import sys
import argparse
import shlex
import time
import random
from datetime import date, datetime
import json

import pandas as pd
from basico import *

#######################
# GLOBAL VARIABLES #

mparams = None
mcomps = None
mspecs = None
mreacts = None
seednparams = None
seedncomps = None
seednspecs = None
seednreacts = None
newname = None
newmodel = None


#######################
# AUXILIARY FUNCTIONS #

# function to test if string is a number, it's amazing this is not native to python...
def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

# function to check that value is positive integer, helper for argparse
def positive_integer(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid negative value" % value)
    return ivalue

# function to check that value is positive number, helper for argparse
def positive_float(value):
    ivalue = float(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid negative value" % value)
    return ivalue

#  class to verify that value is integer within a range, helper for argparse
def integer_range(astr, min=1, max=10):
    try:
        value = int(astr)
    except:
        raise argparse.ArgumentTypeError('value not an integer')
    if min<= value <= max:
        return value
    else:
        raise argparse.ArgumentTypeError('value not in range %s-%s'%(min,max))

# function to check if a string is an element in the model, optionally ignore compartments
def is_element(candidate, ignore_compartments):
    el=False
    if( (seednparams>0) and (candidate in mparams.index) ):
        el=True
    else:
        if( (seedncomps>0) and (candidate in mcomps.index) ):
            el=True
        else:
            if( (not ignore_compartments) and (seednspecs>0) and (candidate in mspecs.index) ):
                el=True
            else:
                if( (seednreacts>0) and (candidate in mreacts.index) ):
                    el=True
    return el

# function to add a diffusive term to an ode expression
def add_diffusive_term(expression, c, a, b):
    return expression + f' + Values[{c}] * ( {a} - {b} )'

# function to add an irreversible diffusive term to an ode expression
def add_irr_diffusive_term(expression, c, sign, a):
    return expression + f' {sign} Values[{c}] * {a}'

# function to add noise to a value
def addnoise( rate, level, dist):
    if( dist == 'uni' ):
        s = 1.0 + level * random.uniform(-1.0,1.0)
    elif( dist == 'norm' ):
        s = random.gauss(1.0, level)
    else:
        print( 'ERROR: dist must be \'uni\' or \'norm\'')
        exit()
    return rate * s

# function to create a rate law for a regulatory synthesis with k inputs
def add_regulatory_ratelaw(k):
    # initialize with V
    exp = 'V'
    rimap = {'V': 'parameter' }
    for i in range(1, k+1):
        exp = exp + f' * ( ( 1 + ( 1 + a{i} ) * M{i} ^ h{i} ) / ( 1 + M{i} ^ h{i} ) )'
        rimap[f'M{i}'] = 'modifier'
        rimap[f'a{i}'] = 'parameter'
        rimap[f'h{i}'] = 'parameter'
    add_function(name=f'regulated_by_{k}', infix=exp, type='general', mapping=rimap, model=newmodel)

# function to change expression, fixing all references to element names with the appropriate suffix, with option to ignore compartments
def fix_expression(exp, suff, ignore_compartments):
    expression = str(exp)
    #is the full expression an element?
    if( is_element(expression, ignore_compartments) ):
        # just process it and return
        expression = expression + suff
        return expression
    mname = re.findall(r'^CN=Root,Model=(.+?),', expression )
    if( mname ):
        expression = re.sub(r'^CN=Root,Model=(.+?),', fr'CN=Root,Model={newname},', expression )
    # find object names inside []
    evars = re.findall(r'\[(.+?)\]', expression )
    if( evars ):
        for el in evars:
            #check that the variable exists
            if( is_element(el, ignore_compartments) ):
                elnew = el + suff
                expression = re.sub(fr'\[{el}\]', fr'[{elnew}]', expression )
    # find object names inside ()
    evars = re.findall(r'\((.+?)\)', expression )
    if( evars ):
        for el in evars:
            #check that the variable exists
            if( is_element(el, ignore_compartments) ):
                elnew = el + suff
                expression = re.sub(fr'\({el}\)', f'({elnew})', expression )
    # find object names inside () special case of ( something(else) )
    evars = re.findall(r'\((.*\(.*\).*?)\)', expression )
    if( evars ):
        for el in evars:
            #check that the variable exists
            if( is_element(el, ignore_compartments) ):
                elnew = el + suff
                el = re.sub(r'\(', r'\\(', el)
                el = re.sub(r'\)', r'\\)', el)
                expression = re.sub(el,elnew, expression )
    # find object names like R1.Rate, I2.InitialParticleNumber, etc.
    evars = re.findall(r'([^\s\]\)]+?)\.\w', expression )
    if( evars ):
        for el in evars:
            #check that the variable exists
            if( is_element(el, ignore_compartments) ):
                elnew = el + suff
                expression = re.sub(fr'{el}\.(\w)', fr'{elnew}.\1', expression )
    return expression

# function to load and parse a graphviz network file, return list of links as tuples
def read_network(network_file):
    # read file
    try:
        with open(network_file) as f:
            netf = f.read()
    except FileNotFoundError:
        print(f'ERROR: file {network_file} not found')
        exit(-1)
    # create list for links
    ll = []
    # check if this is a digraph
    if( re.search(r'digraph\s*.*{', netf) ):
        ll.append((0,0))
    # check for the odd case of /* and */ appearing inside separate quotes and refuse to parse
    if( re.search(r'"[\s\S]*?\/\*[\s\S]*\*\/[\s\S]*?"', netf) ):
        print(fr' Warning: {network_file} contains identifiers with \* and *\, and cannot be parsed')
        return []
    # remove comments
    netf = re.sub(r'\s*//.*', '', netf)          # // comments
    netf = re.sub(r'\s*#.*', '', netf)           # # comments
    netf = re.sub(r'\/\*[\s\S]*?\*\/', '', netf) # /* */ comments
    # find all simple connections like 1--2 and 1--2--3
    # make sure that comments are ignored
    matches = re.findall( r'"*(\d+)"*(?:\:\w+)*\s*(?=(-[-|>]\s*"*(?:\d+)"*(?:\:\w+)*))', netf)
    if( matches ):
        for match in matches:
            e1 = int(match[0])
            se = re.search(r'(\d+)',match[1])
            if( se is None ):
                print(f' Warning: {network_file} contains invalid node {match[1]}, no connections added')
                return []
            e2 = int(se.group(1))
            if( (e1,e2) not in ll):
                ll.append((e1,e2))
            else:
                print(f' Warning: duplicate entry for connection {e1} to {e2}, ignored')
    # find connections to groups
    gmatches = re.findall( r'"*(\d+)"*(?:\:\w+)*\s*-[-|>]\s*(\{[\d\s;,"]+\})', netf)
    if( gmatches ):
        for gmatch in gmatches:
            e1 = int(gmatch[0])
            se = re.findall(r'(\d+)',gmatch[1])
            if( se is None ):
                print(f' Warning: {network_file} contains invalid group node {gmatch[1]}, no connections added')
                return []
            for s in se:
                e2 = int(s)
                if( e1 == e2 ):
                    print(f' Warning: self-connections not allowed, ignoring {e1} -> {e2}')
                else:
                    if( (e1,e2) not in ll):
                        ll.append((e1,e2))
                    else:
                        print(f' Warning: duplicate entry for connection {e1} to {e2}, ignored')
    # check if list is empty and issue warning
    if( (len(ll) == 0) or ( (len(ll) == 1) and ll[0] == (0,0) ) ):
        print(f' Warning: {network_file} did not contain any valid edges, no connections added')
        if( len(ll) == 1 ):
            del ll[0]
    return ll

############################
# MAIN PROGRAM STARTS HERE #

############
# Strategy:
#  1. parse command line and deal with options
#  2. read original model,
#  3. pre-process the items that will have noise
#  4. copy notes, annotations, and units
#  5. create new model
#  MAIN LOOP, iterating over each model element:
#   6. create parameters, compartments, species without expressions
#   7. create reactions (and fix mappings)
#   8. set expressions for compartments and species
#   9. create events that depende on variable
#  10. loop over remaining events that don't depend on variables
#  11. create medium unit if needed
#  12. create unit connections
#  13. copy task settings
#  14. save model
############

def main():
    global mparams, mcomps, mspecs, mreacts, seednparams, seedncomps, seednspecs, seednreacts, newname, newmodel
    
    #####
    #  1. parsing the command line
    #####

    parser = argparse.ArgumentParser(
                        prog='sbmodelr',
                        description='Convert one SBML/COPASI model into another one containing many replicas of the original.')
    # command line arguments
    parser.add_argument('filename', help='original model file')
    parser.add_argument('rows', type=positive_integer, default=2,
                        help='total number of units or number of rows in a 2D or 3D grid')
    parser.add_argument('columns', nargs='?', type=positive_integer, default=1,
                        help='number of columns in a 2D or 3D grid')
    parser.add_argument('layers', nargs='?', type=positive_integer, default=1,
                        help='number of layers in a 3D grid')
    parser.add_argument('-V', '--version', action='version', version=f'sbmodelr {__version__}')
    parser.add_argument('-q', '--quiet', action='store_true', help='supress information messages')
    parser.add_argument('-o', '--output', help='output filename', metavar='filename')
    parser.add_argument('-n', '--network', help='graphviz network file mapping unit connections', metavar='netfile')
    # types of connection
    parser.add_argument('-t', '--transport', action='append', help='species to be transported between units (mass action)', metavar='species')
    parser.add_argument('--Hill-transport', action='append', help='species to be transported between units (Hill kinetics)', metavar='species')
    parser.add_argument('-g', '--grn', action='append', help='species to be connected between units by regulation of synthesis (eg gene regulatory network)', metavar='species')
    parser.add_argument('-d', '--ode-diffusive', action='append', help='explicit ODE to be coupled between units by diffusive terms', metavar='ode')
    parser.add_argument('-s', '--ode-synaptic', action='append', help='explicit ODE to be coupled between units by chemical synaptic terms', metavar='ode')
    # options for transport and diffusive coupling
    parser.add_argument('-k', '--transport-k', dest='transport_rate', type=positive_float, help='value of rate constant for transport between units', default='1.0', metavar='value')
    parser.add_argument('-c', '--coupling-constant', dest='coupling_constant', type=positive_float, default=1.0, help='value for strength of ODE coupling between units', metavar='value')
    # options for Michaelis-Menten transport
    parser.add_argument('--transport-Km', dest='transport_Km', type=positive_float, help='value of Km for transport between units', default='1.0', metavar='value')
    parser.add_argument('--transport-Vmax', dest='transport_Vmax', type=positive_float, help='value of Vmax for transport between units', default='1.0', metavar='value')
    parser.add_argument('--transport-h', dest='transport_h', type=positive_float, help='value of Hill coefficient for transport between units', default='1.0', metavar='value')
    # options for regulatory synthesis interactions
    parser.add_argument('--grn-V', type=positive_float, help='value of basal rate for connection by regulation of synthesis', default=1.0, metavar='value')
    parser.add_argument('--grn-a', type=float, help='value of activation/inhibition strength for connection by regulation of synthesis', default=1.0, metavar='value')
    parser.add_argument('--grn-h', type=integer_range, metavar='INT[1:10]', help='value of Hill coefficient for connection by regulation of synthesis', default=2)
    # options for synaptic coupling
    parser.add_argument('--synapse-g', type=positive_float, help='value of g for synaptic connections between units', default=1.0, metavar='value')
    parser.add_argument('--synapse-V0', type=float, help='value of V0 for synaptic connections between units', default=-20.0, metavar='value')
    parser.add_argument('--synapse-Vsyn', type=float, help='value of Vsyn for synaptic connections between units', default=20.0, metavar='value')
    parser.add_argument('--synapse-tau-r', type=positive_float, help='value of tau_r for synaptic connections between units', default=0.5, metavar='value')
    parser.add_argument('--synapse-tau-d', type=positive_float, help='value of tau_d for synaptic connections between units', default=10,  metavar='value')
    parser.add_argument('--synapse-link-g', action='store_true', help='link all synapse g to a single value')
    # options for noise
    parser.add_argument('--pn', dest='noisy', help='add noise to parameter with level magnitude, dist: {uni,norm}', nargs=3, metavar=('parameter', 'level', 'dist'), action='append')
    parser.add_argument('--cn', dest='cnoise', help='add noise to all coupling parameters with level magnitude, dist: {uni,norm}', nargs=2, metavar=('level', 'dist'))
    # other options
    parser.add_argument('--add-medium', action='store_true', help='add a medium unit with all transported species')
    parser.add_argument('--medium-volume', type=positive_float, help='volume of medium unit', default='1.0', metavar='value')
    parser.add_argument('--ignore-compartments', action='store_true', help='do not replicate compartments')
    parser.add_argument('--ignore-tasks', action='store_true', help='do not copy over task settings')
    parser.add_argument('--sbml', choices=['l1v2', 'l2v3', 'l2v4', 'l2v5', 'l3v1', 'l3v2'], help='export in SBML format of this level and version.')
    #Vivarium argument option
    parser.add_argument('--vivarium', action='store_true', help='Export JSON and model for Vivarium multiscale simulation')

    # Parse the arguments
    args = parser.parse_args()

    seedmodelfile = args.filename
    r = args.rows
    c = args.columns
    l = args.layers

    ignc = args.ignore_compartments
    # value for mass action transport
    trate = args.transport_rate

    # value for ODE coupling constants
    coupleconst = args.coupling_constant

    # values for Michaelis-Menten transport
    tVmax = args.transport_Vmax
    tKm = args.transport_Km
    th = args.transport_h

    # values for ODE synaptic connections
    taurinit = args.synapse_tau_r
    taudinit = args.synapse_tau_d
    v0init = args.synapse_V0
    vsyninit = args.synapse_Vsyn
    gcinit = args.synapse_g
    linkg = args.synapse_link_g

    # values for GRN synthetic regulation
    grnV = args.grn_V
    grna = args.grn_a
    grnh = args.grn_h

    # medium volume
    mediumVol = args.medium_volume

    # check if level for --cn is a positive float and distribution is allowed
    if( args.cnoise ):
        (level, dist) = args.cnoise
        if( not is_float(level) ):
            print( 'ERROR: \'level\' must be a positive floating point number')
            exit()
        if( float(level) < 0.0 ):
            print( 'ERROR: \'level\' must be a positive floating point number')
            exit()
        if( (dist != 'uni') and (dist != 'norm') ):
            print( 'ERROR: dist must be \'uni\' or \'norm\'')
            exit()
        # check if user wanted g values to be linked...
        if( linkg ):
            print( 'ERROR: --cn and --synapse-link-g options cannot be used together, chose only one!')
            exit()

    # check if level for --pn are positive floats and distribution is allowed
    if( args.noisy ):
        for ns in args.noisy:
            (param,level, dist) = ns
            if( not is_float(level) ):
                print( 'ERROR: \'level\' must be a positive floating point number')
                exit()
            if( float(level) < 0.0 ):
                print( 'ERROR: \'level\' must be a positive floating point number')
                exit()
            if( (dist != 'uni') and (dist != 'norm') ):
                print( 'ERROR: dist must be \'uni\' or \'norm\'')
                exit()

    # check if value for --grn-a lies between -1 and +1
    if( (grna < -1.0) or (grna > 1.0) ):
        print( 'ERROR: \'--grn-a\' value must be in the interval [-1,1]')
        exit()

    # unify the lists of species to transport
    transported = []
    if( args.transport ):
        for sp in args.transport:
            # 'a' means mass action
            transported.append((sp,'a'))
    if( args.Hill_transport ):
        for sp in args.Hill_transport:
            # 'h' means Hill kinetics
            transported.append((sp,'h'))

    # unify the lists of ODEs to couple
    odelink = []
    if( args.ode_diffusive ):
        for ode in args.ode_diffusive:
            # 'd' means diffusive coupling
            odelink.append((ode,'d'))
    if( args.ode_synaptic ):
        for ode in args.ode_synaptic:
            # 's' means synaptic coupling
            odelink.append((ode,'s'))

    # get the base of the input model filename
    base,ext = os.path.splitext(seedmodelfile)

    # if sbml file, then we have to ignore tasks
    if( ext.lower() != '.cps' ):
        args.ignore_tasks = True

    # set sbml level and version (l2v4 is the default)
    sbmll = 2
    sbmlv = 4
    if( args.sbml ):
        slv = re.match( r'l(\d)v(\d)', args.sbml )
        if( slv ):
            sbmll = int(slv.group(1))
            sbmlv = int(slv.group(2))
        # if we are writing out sbml, then we also can ignore tasks
        args.ignore_tasks = True

    # sanity check
    nmodels = r*c*l

    if(nmodels==1):
        print("ERROR: Nothing to do, one copy only is the same as the original model!\nAt least one of rows, columns or layers must be larger than 1.\n")
        exit()

    # check dimensionality
    dim=0
    if(r>1):
        dim = dim + 1
    if(c>1):
        dim = dim + 1
    if(l>1):
        dim = dim + 1

    # strings to add to comments and titles, etc
    if(dim==1):
        fsuff = f"{nmodels}"
        desc = f"a set of {nmodels} replicas"
        apdx1 = '_1'
        gridr = nmodels
        gridc = 1
        gridl = 1
    else:
        if(dim==3):
            gridr = r
            gridc = c
            gridl = l
            fsuff = f"{gridr}x{gridc}x{gridl}"
            desc = f"a 3D set of {nmodels} replicas ({gridr}x{gridc}x{gridl})"
            apdx1 = '_1,1,1'
        else:
            if(r==1):
                gridr = c
                gridc = l
                gridl = 1
                fsuff = f"{gridr}x{gridc}"
                desc = f"a 2D set of {nmodels} replicas ({gridr}x{gridc})"
                apdx1 = '_1,1'
            else:
                if(c==1):
                    gridr = r
                    gridc = l
                    gridl = 1
                    fsuff = f"{gridr}x{gridc}"
                    desc = f"a 2D set of {nmodels} replicas ({gridr}x{gridc})"
                    apdx1 = '_1,1'
                else:
                    gridr = r
                    gridc = c
                    gridl = 1
                    fsuff = f"{gridr}x{gridc}"
                    desc = f"a 2D set of {nmodels} replicas ({gridr}x{gridc})"
                    apdx1 = '_1,1'

    # parse the network file if needed
    if( args.network ):
        # check dimension
        if( dim > 1 ):
            print(f'ERROR: network file is only relevant for dimension 1 but you chose dimension {dim} ({fsuff})')
            exit()
        # parse the network file
        links = read_network(args.network)
        digraph = False
        if( len(links) > 0 ):
            if( links[0]==(0,0) ):
                digraph = True
                del links[0]
            for link in links:
                if( (link[0] < 1) or (link[1] < 1) or (link[0] > nmodels) or (link[1] > nmodels) ):
                    print(f'ERROR: network file lists nodes with numbers outside [1,{nmodels}]')
                    exit(1)
        # if we have a grn, we need to process edges and nodes a bit more
        if( args.grn ):
            # if it is not a digraph we should stop right here!
            if not digraph:
                print(f'ERROR: regulatory connections require a directed graph')
                exit(2)
            # dictionary of nodes and their in-degree (doesn't include nodes without inward edges)
            regulated = dict()
            # dictionary of nodes and their regulators
            reglinks = dict()
            for link in links:
                # add one more value to the receiving node
                regulated[link[1]] = regulated.get(link[1], 0) + 1
                reglinks[link[1]] = reglinks.get(link[1], []) + [link[0]]
            # dictionary of in-degrees and nodes with that in-degree (useful for creating rate laws)
            indegrees = dict()
            for k, v in regulated.items():
                indegrees[v] = indegrees.get(v, []) + [k]
            # indegrees is the inverse of regulated...


    #####
    #  2. read the original model
    #####

    seedmodel = load_model(seedmodelfile, remove_user_defined_functions=True)
    if( seedmodel is None):
        print(f'ERROR: {seedmodelfile} failed to load.\n')
        exit(-1)

    # print some information about the model
    if( not args.quiet ):
        print(f'Processing {seedmodelfile}')

    #Get the species first to check for duplicate names and rename them
    # this is not nice but needs to be done to avoid all sorts of troubles...
    mspecs = get_species(model=seedmodel, exact=True)
    if mspecs is not None:
        dups = mspecs.groupby('name')
        for i, df in dups:
            count = len(df)
            if count > 1:
                for _, j in df.iterrows():
                    nn = i + '_' + j['compartment']
                    if count > 1:
                        set_species(j['display_name'], new_name=nn, exact=True, model=seedmodel)
                    else:
                        set_species(i, new_name=nn, exact=True, model=seedmodel)
                    count -= 1

    #Get the global quantities
    mparams = get_parameters(model=seedmodel, exact=True)
    if( mparams is None):
        seednparams = 0
        pfixed = 0
        passg = 0
        pode = 0
    else:
        seednparams = mparams.shape[0]
        # count subsets (fixed, assignment, ode)
        pfixed = (mparams['type']=='fixed').sum()
        passg = (mparams['type']=='assignment').sum()
        pode = (mparams['type']=='ode').sum()

    #Get the compartments
    mcomps = get_compartments(model=seedmodel, exact=True)
    if( mcomps is None):
        seedncomps = 0
        cfixed = 0
        cassg = 0
        code = 0
    else:
        seedncomps = mcomps.shape[0]
        # count subsets (fixed, assignment, ode)
        cfixed = (mcomps['type']=='fixed').sum()
        cassg = (mcomps['type']=='assignment').sum()
        code = (mcomps['type']=='ode').sum()

    #Get the species
    mspecs = get_species(model=seedmodel, exact=True)
    if( mspecs is None):
        seednspecs = 0
        sreact = 0
        sfixed = 0
        sassg = 0
        sode = 0
    else:
        seednspecs = mspecs.shape[0]
        # count subsets (fixed, assignment, ode)
        sreact = (mspecs['type']=='reactions').sum()
        sfixed = (mspecs['type']=='fixed').sum()
        sassg = (mspecs['type']=='assignment').sum()
        sode = (mspecs['type']=='ode').sum()

    # get the reactions
    mreacts = get_reactions(model=seedmodel, exact=True)
    if( mreacts is None):
        seednreacts = 0
    else:
        seednreacts = mreacts.shape[0]

    # get the events
    mevents = get_events(model=seedmodel, exact=True)
    if( mevents is None):
        seednevents = 0
    else:
        seednevents = mevents.shape[0]

    # create string for summary of base model
    if( not args.quiet ):
        base_model_summary = f"  Reactions:         {seednreacts}\n"
        base_model_summary = base_model_summary + f"  Species:           {seednspecs}\t(Reactions: {sreact}, Fixed: {sfixed}, Assignment: {sassg}, ODE: {sode})\n"
        base_model_summary = base_model_summary + f"  Compartments:      {seedncomps}\t(Fixed: {cfixed}, Assignment: {cassg}, ODE: {code})\n"
        base_model_summary = base_model_summary + f"  Global quantities: {seednparams}\t(Fixed: {pfixed}, Assignment: {passg}, ODE: {pode})\n"
        # we print the events later to be able to discriminate how many are only time dependent

    # read scan items
    # we need to retrieve then now before we create a new model due to a bug in COPASI/BasiCO (not sure which)
    scanitems = get_scan_items(model=seedmodel)

    #####
    #  3. pre-process the items that will have noise
    #####

    # deal with noisy parameters, put them in 3 separate lists
    # TODO: apply to classes of parameters ?
    # TODO:  b) reaction parameters, c) initial state, d) compartment volumes

    noisy_species = {}
    noisy_param = {}
    noisy_comp = {}
    if( args.noisy ):
        for (item, level, dist) in args.noisy:
            # species?
            if( (seednspecs>0) and (item in mspecs.index) ):
                noisy_species[item] = (level,dist)
            # global quantity?
            elif( (seednparams>0) and (item in mparams.index) ):
                noisy_param[item] = (level,dist)
            # compartment?
            elif( (seedncomps>0) and (item in mcomps.index) ):
                if( ignc ):
                    print( f' Warning: compartment {item} will not have noise added given ignore-compartments option')
                else:
                    noisy_comp[item] = (level,dist)
            else:
                print( f' Warning: {item} is not an item allowed to be noisy, ignored')

    #####
    #  4. copy notes, annotations, and units
    #####

    seedname = get_model_name(model=seedmodel)

    # edit the notes
    # get the command line
    cmd = " ".join(map(shlex.quote, sys.argv))
    nnotes = get_notes(model=seedmodel)
    # check if notes are empty
    if not nnotes:
        nnotes = f"<body xmlns=\"http://www.w3.org/1999/xhtml\"><p><br/></p><hr/><p>Processed with sbmodelr to produce {desc} of {seedmodelfile}</p><pre style=\"font-size:small\">{cmd}</pre></body>"
    else:
        # check if the notes are in HTML
        index = nnotes.find('<body')
        if( index == -1 ):
            # not HTML, so add a simple string
            nnotes = f"Processed with sbmodelr to produce {desc} of {seedmodelfile}\n{cmd}\n\n ----- \n" + nnotes
        else:
            # find end of body tag
            eobt = nnotes.find('>')
            if( eobt > -1):
                eobt = eobt+1
                # add info at the top of the body section
                nend = nnotes[eobt:]
                nnotes = nnotes[:eobt] + f"<p>Processed with sbmodelr to produce {desc} of {seedmodelfile}</p><pre style=\"font-size:small\">{cmd}</pre><p>notes of original file below:</p><hr/>" + nend
            else:
                # something went wrong, let's just write comments as if there was nothing
                nnotes = f"<body xmlns=\"http://www.w3.org/1999/xhtml\"><p><br/></p><hr/><p>Processed with sbmodelr to produce {desc} of {seedmodelfile}</p><pre style=\"font-size:small\">{cmd}</pre></body>"


    # get original model units
    munits = get_model_units(model=seedmodel)

    #####
    #  5. create new model
    #####

    # create new model filename (if sbml string was given force .xml extension)
    if( args.output ):
        newfilename = args.output
    else:
        base,ext = os.path.splitext(os.path.basename(seedmodelfile))
        if( args.sbml ):
            newfilename = f"{base}_{fsuff}.xml"
        else:
            newfilename = f"{base}_{fsuff}{ext}"

    # create the new model name
    newname = f"{desc} of {seedname}"

    # create the new model
    newmodel = new_model(name=newname,
                        notes=nnotes,
                        quantity_unit=munits['quantity_unit'],
                        time_unit=munits['time_unit'],
                        volume_unit=munits['volume_unit'],
                        area_unit=munits['area_unit'],
                        length_unit=munits['length_unit'])

    # set the intial time
    it= get_value('Time', initial=True, model=seedmodel)
    set_value('Time', it, initial=True, model=newmodel)

    # transfer the annotations
    miriam = get_miriam_annotation(model=seedmodel)
    if 'created' in miriam:
        set_miriam_annotation(model=newmodel, created=miriam['created'], replace=True)
    if 'creators' in miriam:
        set_miriam_annotation(model=newmodel, creators=miriam['creators'], replace=True)
    if 'references' in miriam:
        set_miriam_annotation(model=newmodel, references=miriam['references'], replace=True)
    if 'description' in miriam:
        set_miriam_annotation(model=newmodel, description=miriam['description'], replace=True)
    # add one modification now
    if 'modifications' in miriam:
        miriam['modifications'].append(datetime.datetime.now())
        set_miriam_annotation(model=newmodel, modifications=miriam['modifications'], replace=True)
    else:
        modf = []
        modf.append(datetime.datetime.now())
        set_miriam_annotation(model=newmodel, modifications=modf, replace=True)

    # if we are ignoring compartments, then simply copy the compartments to the new model
    if ignc:
        for p in mcomps.index:
            add_compartment(model=newmodel, name=p, status=mcomps.loc[p].at['type'], initial_size=mcomps.loc[p].at['initial_size'], unit=mcomps.loc[p].at['unit'], dimensionality=int(mcomps.loc[p].at['dimensionality']), expression=mcomps.loc[p].at['expression'], initial_expression=mcomps.loc[p].at['initial_expression'] )
            if( 'notes' in mcomps.loc[p] ):
                set_compartment(model=newmodel, name=p, notes=mcomps.loc[p].at['notes'])

    # if we want to have synapses with linked g, create the master g
    if( args.ode_synaptic and linkg ):
        linkedsyng = f'g_c_{ode}_synapse'
        add_parameter(linkedsyng, type='fixed', initial_value=gcinit, model=newmodel)

    #####
    #  MAIN LOOP FOR REPLICATION
    #####

    # we use "_i" as suffix for 1D, "_r,c" for 2D and "_r,c,l" for 3D
    i = 0
    for r in range(gridr):
        for c in range(gridc):
            for l in range(gridl):
                if(dim==1):
                    apdx = f"_{i+1}"
                else:
                    if(dim==2):
                        apdx = f"_{r+1},{c+1}"
                    else:
                        apdx = f"_{r+1},{c+1},{l+1}"

    #####
    #  6. create parameters, compartments and species
    #####
                # PARAMETERS
                if( seednparams>0 ):
                    for p in mparams.index:
                        nname = p + apdx
                        iv = mparams.loc[p].at['initial_value']
                        if( p in noisy_param ):
                            (level,dist) = noisy_param[p]
                            iv = addnoise(mparams.loc[p].at['initial_value'], float(level), dist)
                        add_parameter(model=newmodel, name=nname, status='fixed', initial_value=iv, unit=mparams.loc[p].at['unit'])
                        nt = get_notes(model=seedmodel, name=f'Values[{p}]')
                        if( nt is not None ):
                            set_parameters(model=newmodel, exact=True, name=nname, notes=nt)
                        an = get_miriam_annotation(model=seedmodel, name=f'Values[{p}]')
                        if( an ):
                            if( 'creators' in an ):
                                set_miriam_annotation(creators=an['creators'],model=newmodel, name=nname, replace=False)
                            if( 'references' in an ):
                                set_miriam_annotation(references=an['references'],model=newmodel, name=nname, replace=False)
                            if( 'descriptions' in an ):
                                set_miriam_annotation(descriptions=an['descriptions'],model=newmodel, name=nname, replace=False)
                            if( 'modifications' in an ):
                                set_miriam_annotation(modifications=an['modifications'],model=newmodel, name=nname, replace=False)
                            if( 'created' in an ):
                                set_miriam_annotation(created=an['created'],model=newmodel, name=nname, replace=False)
                # COMPARTMENTS
                # if we are ignore_compartments, then we already created the original ones, nothing done here
                if( (seedncomps > 0) and (not ignc) ):
                    for p in mcomps.index:
                        iv = mcomps.loc[p].at['initial_size']
                        if( p in noisy_comp ):
                            (level,dist) = noisy_comp[p]
                            iv = addnoise(mcomps.loc[p].at['initial_size'], float(level), dist)
                        nname = p + apdx
                        add_compartment(model=newmodel, name=nname, status=mcomps.loc[p].at['type'], initial_size=iv, unit=mcomps.loc[p].at['unit'], dimensionality=int(mcomps.loc[p].at['dimensionality']) )
                        nt = get_notes(model=seedmodel, name=f'Compartments[{p}]')
                        if( nt is not None ):
                            set_compartment(model=newmodel, name=nname, notes=nt)
                        an = get_miriam_annotation(model=seedmodel, name=f'Compartments[{p}]')
                        if( an ):
                            if( 'creators' in an ):
                                set_miriam_annotation(creators=an['creators'],model=newmodel, name=nname, replace=False)
                            if( 'references' in an ):
                                set_miriam_annotation(references=an['references'],model=newmodel, name=nname, replace=False)
                            if( 'descriptions' in an ):
                                set_miriam_annotation(descriptions=an['descriptions'],model=newmodel, name=nname, replace=False)
                            if( 'modifications' in an ):
                                set_miriam_annotation(modifications=an['modifications'],model=newmodel, name=nname, replace=False)
                            if( 'created' in an ):
                                set_miriam_annotation(created=an['created'],model=newmodel, name=nname, replace=False)
                # SPECIES
                if( seednspecs > 0):
                    for p in mspecs.index:
                        iv = mspecs.loc[p].at['initial_concentration']
                        if( p in noisy_species ):
                            (level,dist) = noisy_species[p]
                            iv = addnoise(mspecs.loc[p].at['initial_concentration'], float(level), dist)
                        nname = p + apdx
                        if ignc:
                            cp = mspecs.loc[p].at['compartment']
                        else:
                            cp = mspecs.loc[p].at['compartment'] + apdx
                        add_species(model=newmodel, name=nname, compartment_name=cp, status=mspecs.loc[p].at['type'], initial_concentration=iv, unit=mspecs.loc[p].at['unit'])
                        nt = get_notes(model=seedmodel, name=p)
                        if( nt is not None ):
                            set_species(model=newmodel, exact=True, name=nname, notes=nt)
                        an = get_miriam_annotation(model=seedmodel, name=p)
                        if( an ):
                            if( 'creators' in an ):
                                set_miriam_annotation(creators=an['creators'],model=newmodel, name=nname, replace=False)
                            if( 'references' in an ):
                                set_miriam_annotation(references=an['references'],model=newmodel, name=nname, replace=False)
                            if( 'descriptions' in an ):
                                set_miriam_annotation(descriptions=an['descriptions'],model=newmodel, name=nname, replace=False)
                            if( 'modifications' in an ):
                                set_miriam_annotation(modifications=an['modifications'],model=newmodel, name=nname, replace=False)
                            if( 'created' in an ):
                                set_miriam_annotation(created=an['created'],model=newmodel, name=nname, replace=False)

    #####
    #  7. create reactions
    #####

                # REACTIONS
                if( seednreacts > 0):
                    for p in mreacts.index:
                        nname = p + apdx
                        scheme = mreacts.loc[p].at['scheme']
                        tok = scheme.split(';')
                        tok2 = [shlex.split(sub, posix=False) for sub in tok]
                        # build the reaction string
                        rs = ""
                        for t in tok2[0]:
                            if( (t == '=') or (t == '->') or (t == '+') or is_float(t) or (t=="*")):
                                rs = rs + t + " "
                            else:
                                if re.match(r'\".+\"', t):
                                    t = re.sub( r'\"(.+)\"', f'"\\1{apdx}"', t )
                                    rs = rs + t + " "
                                else:
                                    rs = rs + t + apdx + " "
                        if( len(tok2) > 1 ):
                            # deal with the modifiers
                            rs = rs[:len(rs)-1] + "; "
                            for t in tok2[1]:
                                if re.match(r'\".+\"', t):
                                    t = re.sub( r'\"(.+)\"', f'"\\1{apdx}"', t )
                                    rs = rs + t + " "
                                else:
                                    rs = rs + t + apdx + " "
                        # fix the parameter mappings
                        mapp = mreacts.loc[p].at['mapping'].copy()
                        for key in mapp:
                            if( isinstance(mapp[key], str) ):
                                t = mapp[key]
                                if re.match(r'\".+\"', t):
                                    t = re.sub( r'\"(.+)\"', f'"\\1{apdx}"', t )
                                else:
                                    t = t + apdx
                                mapp[key] = t
                            else:
                                if( isinstance(mapp[key], list ) ):
                                    nmk = []
                                    for k2 in mapp[key]:
                                        if re.match(r'\".+\"', k2):
                                            k2 = re.sub( r'\"(.+)\"', f'"\\1{apdx}"', k2 )
                                        else:
                                            k2 = k2 + apdx
                                        nmk.append(k2)
                                    mapp[key] = nmk
                                    #mapp[key] = [k2 + apdx for k2 in mapp[key]]
                        add_reaction(model=newmodel, name=nname, scheme=rs, mapping=mapp, function=mreacts.loc[p].at['function'] )
                        nt = get_notes(model=seedmodel, name=p)
                        if( nt is not None ):
                            set_reaction(model=newmodel, exact=True, name=nname, notes=nt)
                        an = get_miriam_annotation(model=seedmodel, name=p)
                        if( an ):
                            if( 'creators' in an ):
                                set_miriam_annotation(creators=an['creators'],model=newmodel, name=nname, replace=False)
                            if( 'references' in an ):
                                set_miriam_annotation(references=an['references'],model=newmodel, name=nname, replace=False)
                            if( 'descriptions' in an ):
                                set_miriam_annotation(descriptions=an['descriptions'],model=newmodel, name=nname, replace=False)
                            if( 'modifications' in an ):
                                set_miriam_annotation(modifications=an['modifications'],model=newmodel, name=nname, replace=False)
                            if( 'created' in an ):
                                set_miriam_annotation(created=an['created'],model=newmodel, name=nname, replace=False)
    #####
    #  8. set expressions and initial_expressions
    #####

                # PARAMETERS
                if( seednparams > 0 ):
                    for p in mparams.index:
                        nname = p + apdx
                        if( mparams.loc[p].at['initial_expression'] ):
                            ie = fix_expression(mparams.loc[p].at['initial_expression'], apdx, ignc)
                            set_parameters(model=newmodel, name=nname, exact=True, initial_expression=ie )
                        if( mparams.loc[p].at['type']=='assignment' or mparams.loc[p].at['type']=='ode'):
                            ex = fix_expression(mparams.loc[p].at['expression'], apdx, ignc)
                            set_parameters(model=newmodel, name=nname, exact=True, status=mparams.loc[p].at['type'], expression=ex )
                # COMPARTMENTS
                if( (seedncomps > 0) and (not ignc) ):
                    for p in mcomps.index:
                        nname = p + apdx
                        if( mcomps.loc[p].at['initial_expression'] ):
                            ie = fix_expression(mcomps.loc[p].at['initial_expression'], apdx, ignc)
                            set_compartment(model=newmodel, name=nname, exact=True, initial_expression=ie )
                        if( mcomps.loc[p].at['type']=='assignment' or mcomps.loc[p].at['type']=='ode'):
                            ex = fix_expression(mcomps.loc[p].at['expression'], apdx, ignc)
                            set_compartment(model=newmodel, name=nname, exact=True, expression=ex )
                # SPECIES
                if( seednspecs > 0):
                    for p in mspecs.index:
                        nname = p + apdx
                        if( mspecs.loc[p].at['initial_expression'] ):
                            ie = fix_expression(mspecs.loc[p].at['initial_expression'], apdx, ignc)
                            set_species(model=newmodel, name=nname, exact=True, initial_expression=ie )
                        if( mspecs.loc[p].at['type']=='assignment' or mspecs.loc[p].at['type']=='ode'):
                            ex = fix_expression(mspecs.loc[p].at['expression'], apdx, ignc)
                            set_species(model=newmodel, name=nname, exact=True, expression=ex )

    #####
    #  9. create events
    #####

                # EVENTS
                timeonlyevents = []
                if( seednevents > 0):
                    for p in mevents.index:
                        # fix the trigger expression
                        tr = fix_expression(mevents.loc[p].at['trigger'], apdx, ignc)
                        # we skip events that have no elements in the trigger (time-dependent only)
                        if(tr != mevents.loc[p].at['trigger']):
                            # fix name
                            nm = p + apdx
                            # process the targets and expressions
                            assg = []
                            for a in mevents.loc[p].at['assignments']:
                                assg.append((fix_expression(a['target'],apdx, ignc),  fix_expression(a['expression'],apdx, ignc)))
                            # add the event
                            add_event(model=newmodel, name=nm, trigger=tr, assignments=assg, delay=fix_expression(mevents.loc[p].at['delay'],apdx, ignc), priority=fix_expression(mevents.loc[p].at['priority'],apdx, ignc), persistent=mevents.loc[p].at['persistent'], fire_at_initial_time=mevents.loc[p].at['fire_at_initial_time'], delay_calculation=mevents.loc[p].at['delay_calculation'])
                            nt = get_notes(model=seedmodel, name=p)
                            if( nt is not None ):
                                set_notes(model=newmodel, name=nm, notes=nt)
                            an = get_miriam_annotation(model=seedmodel, name=p)
                            if( an ):
                                if( 'creators' in an ):
                                    set_miriam_annotation(creators=an['creators'],model=newmodel, name=nm, replace=False)
                                if( 'references' in an ):
                                    set_miriam_annotation(references=an['references'],model=newmodel, name=nm, replace=False)
                                if( 'descriptions' in an ):
                                    set_miriam_annotation(descriptions=an['descriptions'],model=newmodel, name=nm, replace=False)
                                if( 'modifications' in an ):
                                    set_miriam_annotation(modifications=an['modifications'],model=newmodel, name=nm, replace=False)
                                if( 'created' in an ):
                                    set_miriam_annotation(created=an['created'],model=newmodel, name=nm, replace=False)
                        else:
                            # the trigger does not involve any model element other than time
                            # (or compartments when ignore_compartments is on)
                            # add it to the list to be dealt with later
                            timeonlyevents.append(p)

                i += 1

    #####
    #  10. create events not dependent on variables
    #####

    etd=len(timeonlyevents)
    entd = seednevents - etd

    # now we can print out the rest of the model information
    if( not args.quiet ):
        #how many events there are...
        base_model_summary = base_model_summary + f"  Events:            {seednevents}\t(Only time-dependent: {etd}, variable-dependent: {entd})\n"
        # and also information on network if one is used
        if( args.network ):
            if digraph:
                netype = 'directed'
            else:
                netype = 'undirected'
            base_model_summary = base_model_summary + f"  Network:           {netype}\t(Nodes: {nmodels}, Edges: {len(links)})\n"
        # let's print it
        print( base_model_summary )

    # let's go over the events again to process those that are only time dependent
    if( etd > 0 ):
        # loop over the time-only dependent events
        for p in timeonlyevents:
            # if the delay or priority expressions contain elements we use model_1
            dl = fix_expression(mevents.loc[p].at['delay'],apdx1, ignc)
            pr = fix_expression(mevents.loc[p].at['priority'],apdx1, ignc)
            if( not args.quiet ):
                if( dl != mevents.loc[p].at['delay'] ):
                    print(f' Warning: Event {p} contains a delay expression dependent on variables, it was set to the variables of unit {apdx1}')
                if( pr != mevents.loc[p].at['priority'] ):
                    print(f' Warning: Event {p} contains a priority expression dependent on variables, it was set to the variables of unit {apdx1}')
            # process the targets and expressions
            assg = []
            for a in mevents.loc[p].at['assignments']:
                # now loop over all replicates to duplicate the targets
                i = 0
                for r in range(gridr):
                    for c in range(gridc):
                        for l in range(gridl):
                            if(dim==1):
                                apdx = f"_{i+1}"
                            else:
                                if(dim==2):
                                    apdx = f"_{r+1},{c+1}"
                                else:
                                    apdx = f"_{r+1},{c+1},{l+1}"
                            # add the assignment
                            assg.append((fix_expression(a['target'],apdx, ignc), fix_expression(a['expression'],apdx, ignc)))
                            i = i + 1
            # add the event
            add_event(model=newmodel, name=p, trigger=mevents.loc[p].at['trigger'], assignments=assg, delay=dl, priority=pr, persistent=mevents.loc[p].at['persistent'], fire_at_initial_time=mevents.loc[p].at['fire_at_initial_time'], delay_calculation=mevents.loc[p].at['delay_calculation'] )
            nt = get_notes(model=seedmodel, name=p)
            if( nt is not None ):
                set_notes(model=newmodel, name=p, notes=nt)
            an = get_miriam_annotation(model=seedmodel, name=p)
            if( an ):
                if( 'creators' in an ):
                    set_miriam_annotation(creators=an['creators'],model=newmodel, name=p, replace=False)
                if( 'references' in an ):
                    set_miriam_annotation(references=an['references'],model=newmodel, name=p, replace=False)
                if( 'descriptions' in an ):
                    set_miriam_annotation(descriptions=an['descriptions'],model=newmodel, name=p, replace=False)
                if( 'modifications' in an ):
                    set_miriam_annotation(modifications=an['modifications'],model=newmodel, name=p, replace=False)
                if( 'created' in an ):
                    set_miriam_annotation(created=an['created'],model=newmodel, name=p, replace=False)

    #####
    # 11. create medium unit if needed
    #####
    if( args.add_medium ):
        if( transported or odelink ):
            # check if there is no compartment called "medium" already
            medium_name = 'medium'
            if((mcomps is not None) and (medium_name in mcomps.index)):
                medium_name = '_added_medium_'
            # create medium compartment
            add_compartment(model=newmodel, name=medium_name, status='fixed', initial_size=mediumVol, dimiensionality=3, notes="medium compartment added by sbmodlr" )
            # create the species that are transported
            if( transported ):
                for (sp,ttype) in transported:
                    nname = f'{sp}_medium'
                    add_species(model=newmodel, name=nname, compartment_name=medium_name, status='reactions', initial_concentration=mspecs.loc[sp].at['initial_concentration'], unit=mspecs.loc[sp].at['unit'] )
                    if( 'notes' in mspecs.loc[p] ):
                        set_species(model=newmodel, name=nname, notes=mspecs.loc[p].at['notes'])

            # it would be logic to create odes here too, but it is easier to create them further down
            # because they can be of one of three different types of entity
        else:
            print(' Warning: no medium unit created because no species are being transported or ODEs coupled')


    #####
    # 12. create unit connections
    #####

    # check if we need to add the Hill transport rate law
    if( args.Hill_transport ):
        # add Hill kinetics for molecule transport
        htmap = {'V': 'parameter', 'Km': 'parameter', 'h': 'parameter', 'S': 'substrate', 'P': 'product'}
        add_function(name='Hill Transport', infix='V * ( S ^ h - P ^ h ) / ( Km ^ h + S ^ h + P ^ h )', type='reversible', mapping=htmap, model=newmodel)

    # species to be transported
    if( transported ):
        for (sp,ttype) in transported:
            # check that the species exists
            if( (seednspecs>0) and (sp in mspecs.index) ):
                # check that the species depends on reactions
                if( mspecs.loc[sp].at['type'] != 'reactions' ):
                    print( f'ERROR: {sp} is a species that does not depend on reactions, no transport reactions can be added')
                    exit()
                if( ttype == 'a' ):
                    # add a rate constant for the transport reactions
                    rateconst = f'k_{sp}_transport'
                    if( not args.cnoise ):
                        add_parameter(name=rateconst, initial_value=trate, model=newmodel)
                elif( ttype == 'h' ):
                    # add parameters for the transport reactions
                    vmaxconst = f'Vmax_{sp}_transport'
                    if( not args.cnoise ):
                        add_parameter(name=vmaxconst, initial_value=tVmax, model=newmodel)
                    kmconst = f'Km_{sp}_transport'
                    add_parameter(name=kmconst, initial_value=tKm, model=newmodel)
                    hconst = f'h_{sp}_transport'
                    add_parameter(name=hconst, initial_value=th, model=newmodel)
                else:
                    print(f'ERROR: species transport of type \'{ttype}\' is not allowed.')
                    exit()
                # add a transport reaction for each neighbour
                if( dim == 1):
                    # add transport between species and the medium which is always reversible
                    if(args.add_medium):
                        for r in range(gridr):
                            suffa = f'{r+1}'
                            rname = f't_{sp}_{suffa}-medium'
                            rscheme = f'{sp}_{suffa} = {sp}_medium'
                            if( ttype == 'a' ):
                                thisrateconst = rateconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisrateconst = f'k_{sp}_transport_{suffa}-medium'
                                    v = addnoise(trate,float(level),dist)
                                    add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_medium'}
                                add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                            elif( ttype == 'h' ):
                                thisvmaxconst = vmaxconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 vmax per reaction
                                    (level, dist) = args.cnoise
                                    thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-medium'
                                    v = addnoise(tVmax,float(level),dist)
                                    add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_medium'}
                                add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )

                    if( args.network ):
                        for link in links:
                            suffa = f'{link[0]}'
                            suffb = f'{link[1]}'
                            rname = f't_{sp}_{suffa}-{suffb}'
                            # check whether reversible or irreversible
                            if digraph:
                                # check if we have a self-connection and ignore if type a or m
                                if( (suffa == suffb) and ( ttype == 'a' or ttype == 'h' ) ):
                                    print(f' Warning: transport on the same unit not allowed, ignoring {suffa} -> {suffb}')
                                    continue
                                rscheme = f'{sp}_{suffa} -> {sp}_{suffb}'
                                if( ttype == 'a' ):
                                    thisrateconst = rateconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisrateconst = f'k_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(trate,float(level),dist)
                                        add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                    rmap = {'k1': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_{suffb}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (irreversible)' )
                                elif( ttype == 'h' ):
                                    thisvmaxconst = vmaxconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 vmax per reaction
                                        (level, dist) = args.cnoise
                                        thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(tVmax,float(level),dist)
                                        add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                    rmap = {'V': thisvmaxconst, 'Shalve': kmconst, 'h': hconst, 'substrate': f'{sp}_{suffa}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Cooperativity' )

                            else:
                                # check if we have a self-connection and ignore if type a or m
                                if( (suffa == suffb) and ( ttype == 'a' or ttype == 'h' ) ):
                                    print(f' Warning: transport on the same unit not allowed, ignoring {suffa} -- {suffb}')
                                rscheme = f'{sp}_{suffa} = {sp}_{suffb}'
                                if( ttype == 'a' ):
                                    thisrateconst = rateconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisrateconst = f'k_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(trate,float(level),dist)
                                        add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                    rmap = {'k1': thisrateconst, 'k2': rateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_{suffb}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                elif( ttype == 'h' ):
                                    thisvmaxconst = vmaxconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 vmax per reaction
                                        (level, dist) = args.cnoise
                                        thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(tVmax,float(level),dist)
                                        add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                    rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_{suffb}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                elif( dim == 2 ):
                    for r in range(gridr):
                        for c in range(gridc):
                            suffa = f'{r+1},{c+1}'
                            if(args.add_medium):
                                rname = f't_{sp}_{suffa}-medium'
                                rscheme = f'{sp}_{suffa} = {sp}_medium'
                                if( ttype == 'a' ):
                                    thisrateconst = rateconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisrateconst = f'k_{sp}_transport_{suffa}-medium'
                                        v = addnoise(trate,float(level),dist)
                                        add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                    rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_medium'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                elif( ttype == 'h' ):
                                    thisvmaxconst = vmaxconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 vmax per reaction
                                        (level, dist) = args.cnoise
                                        thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-medium'
                                        v = addnoise(tVmax,float(level),dist)
                                        add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                    rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_medium'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                            if( c+1 < gridc ):
                                suffb = f'{r+1},{c+2}'
                                rname = f't_{sp}_{suffa}-{suffb}'
                                rscheme = f'{sp}_{suffa} = {sp}_{suffb}'
                                if( ttype == 'a' ):
                                    thisrateconst = rateconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisrateconst = f'k_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(trate,float(level),dist)
                                        add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                    rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_{suffb}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                elif( ttype == 'h' ):
                                    thisvmaxconst = vmaxconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 vmax per reaction
                                        (level, dist) = args.cnoise
                                        thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(tVmax,float(level),dist)
                                        add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                    rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_{suffb}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                            if( r+1 < gridr ):
                                suffb = f'{r+2},{c+1}'
                                rname = f't_{sp}_{suffa}-{suffb}'
                                rscheme = f'{sp}_{suffa} = {sp}_{suffb}'
                                if( ttype == 'a' ):
                                    thisrateconst = rateconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisrateconst = f'k_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(trate,float(level),dist)
                                        add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                    rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_{suffb}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                elif( ttype == 'h' ):
                                    thisvmaxconst = vmaxconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 vmax per reaction
                                        (level, dist) = args.cnoise
                                        thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-{suffb}'
                                        v = addnoise(tVmax,float(level),dist)
                                        add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                    rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_{suffb}'}
                                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                elif( dim == 3 ):
                    for r in range(gridr):
                        for c in range(gridc):
                            for l in range(gridl):
                                suffa = f'{r+1},{c+1},{l+1}'
                                if(args.add_medium):
                                    rname = f't_{sp}_{suffa}-medium'
                                    rscheme = f'{sp}_{suffa} = {sp}_medium'
                                    if( ttype == 'a' ):
                                        thisrateconst = rateconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 parameter per reaction
                                            (level, dist) = args.cnoise
                                            thisrateconst = f'k_{sp}_transport_{suffa}-medium'
                                            v = addnoise(trate,float(level),dist)
                                            add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                        rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_medium'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                    elif( ttype == 'h' ):
                                        thisvmaxconst = vmaxconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 vmax per reaction
                                            (level, dist) = args.cnoise
                                            thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-medium'
                                            v = addnoise(tVmax,float(level),dist)
                                            add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                        rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_medium'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                                if( c+1 < gridc ):
                                    suffb = f'{r+1},{c+2},{l+1}'
                                    rname = f't_{sp}_{suffa}-{suffb}'
                                    rscheme = f'{sp}_{suffa} = {sp}_{suffb}'
                                    if( ttype == 'a' ):
                                        thisrateconst = rateconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 parameter per reaction
                                            (level, dist) = args.cnoise
                                            thisrateconst = f'k_{sp}_transport_{suffa}-{suffb}'
                                            v = addnoise(trate,float(level),dist)
                                            add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                        rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_{suffb}'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                    elif( ttype == 'h' ):
                                        thisvmaxconst = vmaxconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 vmax per reaction
                                            (level, dist) = args.cnoise
                                            thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-{suffb}'
                                            v = addnoise(tVmax,float(level),dist)
                                            add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                        rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_{suffb}'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                                if( r+1 < gridr ):
                                    suffb = f'{r+2},{c+1},{l+1}'
                                    rname = f't_{sp}_{suffa}-{suffb}'
                                    rscheme = f'{sp}_{suffa} = {sp}_{suffb}'
                                    if( ttype == 'a' ):
                                        thisrateconst = rateconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 parameter per reaction
                                            (level, dist) = args.cnoise
                                            thisrateconst = f'k_{sp}_transport_{suffa}-{suffb}'
                                            v = addnoise(trate,float(level),dist)
                                            add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                        rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_{suffb}'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                    elif( ttype == 'h' ):
                                        thisvmaxconst = vmaxconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 vmax per reaction
                                            (level, dist) = args.cnoise
                                            thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-{suffb}'
                                            v = addnoise(tVmax,float(level),dist)
                                            add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                        rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_{suffb}'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                                if( l+1 < gridl ):
                                    suffb = f'{r+1},{c+1},{l+2}'
                                    rname = f't_{sp}_{suffa}-{suffb}'
                                    rscheme = f'{sp}_{suffa} = {sp}_{suffb}'
                                    if( ttype == 'a' ):
                                        thisrateconst = rateconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 parameter per reaction
                                            (level, dist) = args.cnoise
                                            thisrateconst = f'k_{sp}_transport_{suffa}-{suffb}'
                                            v = addnoise(trate,float(level),dist)
                                            add_parameter(name=thisrateconst, initial_value=v, model=newmodel)
                                        rmap = {'k1': thisrateconst, 'k2': thisrateconst, 'substrate': f'{sp}_{suffa}', 'product': f'{sp}_{suffb}'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Mass action (reversible)' )
                                    elif( ttype == 'h' ):
                                        thisvmaxconst = vmaxconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 vmax per reaction
                                            (level, dist) = args.cnoise
                                            thisvmaxconst = vmaxconst = f'Vmax_{sp}_transport_{suffa}-{suffb}'
                                            v = addnoise(tVmax,float(level),dist)
                                            add_parameter(name=thisvmaxconst, initial_value=v, model=newmodel)
                                        rmap = {'V': thisvmaxconst, 'Km': kmconst, 'h': hconst, 'S': f'{sp}_{suffa}', 'P': f'{sp}_{suffb}'}
                                        add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function='Hill Transport' )
                else:
                    # insanity error
                    print('ERROR: this should not happen, dim is not 1,2 or 3. Contact developer!')
                    exit()
            else:
                # error
                print( f'ERROR: Species {sp} does not exist in the model' )
                exit()

    # ODEs to be coupled by diffusive mechanism or chemical synapses
    if( odelink ):
        for (ode,linktype) in odelink:
            if( linktype=='d'):
                # add a rate constant for the diffusive connections
                diffconst = f'k_{ode}_coupling'
                # only add it if we don't have noise in coupling parameters
                if( not args.cnoise ):
                    add_parameter(diffconst, type='fixed', initial_value=coupleconst, model=newmodel)
            elif( linktype=='s'):
                # check dimension, only dim 1 can have synaptic connections
                if(dim>1):
                    print(f'ERROR: 2D or 3D grids cannot have synaptic connections')
                    exit()
                if( args.network and not digraph ):
                    print( f' Warning: network was defined as undirected, but synapses will be added as directed connections' )
                # add rate constants for the synaptic connections
                # tau_r is a time constant characteristic of pre-synaptic synthesis of neurotransmitter and its diffusion across the synapse
                syntaur = f'tau_r_{ode}_synapse'
                # only add it if we don't have noise in coupling parameters
                if( not args.cnoise ):
                    add_parameter(syntaur, type='fixed', initial_value=taurinit, model=newmodel)
                # tau_d is a time constant characteristic of post-synaptic degradation of neurotransmitter
                syntaud = f'tau_d_{ode}_synapse'
                # only add it if we don't have noise in coupling parameters
                if( not args.cnoise ):
                    add_parameter(syntaud, type='fixed', initial_value=taudinit, model=newmodel)
                # V_0 is the voltage at the pre-synaptic neuron that gives half maximal rate of release of neurotransmitter
                synv0 = f'V0_{ode}_synapse'
                add_parameter(synv0, type='fixed', initial_value=v0init, model=newmodel)
                # V_syn is the post-synaptic inversion potential characteristic
                synvsyn = f'Vsyn_{ode}_synapse'
                # only add it if we don't have noise in coupling parameters
                if( not args.cnoise ):
                    add_parameter(synvsyn, type='fixed', initial_value=vsyninit, model=newmodel)
            else:
                print(f'ERROR: ODE link of type \'{linktype}\' is not allowed.')
                exit()
            # check if the global variable is an ODE
            if( (mparams is not None) and (ode in mparams.index) ):
                if( mparams.loc[ode].at['type'] != 'ode' ):
                    print(f'ERROR: {ode} is a global variable that is not an ODE')
                    exit()
                if( dim == 1):
                    # add coupling between unit and the medium
                    if(args.add_medium and linktype=='d'):
                        # we first need to add the global quantity
                        mediumode = f'{ode}_medium'
                        medexp = ''
                        add_parameter(mediumode, type='ode', expression=medexp, initial_value=mparams.loc[ode].at['initial_value'], model=newmodel)
                        for r in range(gridr):
                            # name of ode in this unit
                            suffa = f'{r+1}'
                            oname = f'{ode}_{suffa}'
                            # get ode to get expression
                            tode = get_parameters(oname, exact=True, model=newmodel)
                            thisdiffconst = diffconst
                            if( args.cnoise ):
                                # if we add noise to couplings, then we need 1 parameter per reaction
                                (level, dist) = args.cnoise
                                thisdiffconst = f'k_{ode}_coupling_{suffa}-medium'
                                v = addnoise(coupleconst,float(level),dist)
                                add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                            # add term to medium ODE
                            medexp = add_diffusive_term(medexp, thisdiffconst,f'Values[{oname}]',f'Values[{mediumode}]')
                            # add term to target ODE
                            odexpr = add_diffusive_term(tode.loc[oname].at['expression'], thisdiffconst,f'Values[{mediumode}]', f'Values[{oname}]')
                            set_parameters(mediumode, exact=True, type='ode', expression=medexp, model=newmodel)
                            set_parameters(oname, exact=True, type='ode', expression=odexpr, model=newmodel)
                    if( args.network ):
                        for link in links:
                            # check if we have a self-connection and ignore it for diffusive connections
                            if( (link[0] == link[1]) and linktype=='d'):
                                print(f' Warning: diffusive coupling onto the same unit not allowed, ignoring {link[0]} -> {link[1]}')
                                continue
                            suffa = f'{link[0]}'
                            oaname = f'{ode}_{suffa}'
                            suffb = f'{link[1]}'
                            obname = f'{ode}_{suffb}'
                            # get ode expressions from a and b
                            aode = get_parameters(oaname, exact=True, model=newmodel)
                            bode = get_parameters(obname, exact=True, model=newmodel)
                            if( linktype=='d' ):
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                if digraph:
                                    # add term to ODE a
                                    odeaexpr = add_irr_diffusive_term(aode.loc[oaname].at['expression'], thisdiffconst,'-',f'Values[{oaname}]')
                                    # add term to ODE b
                                    odebexpr = add_irr_diffusive_term(bode.loc[obname].at['expression'], thisdiffconst,'+',f'Values[{oaname}]')
                                else:
                                    # add term to ODE a
                                    odeaexpr = add_diffusive_term(aode.loc[oaname].at['expression'], thisdiffconst, f'Values[{obname}]', f'Values[{oaname}]')
                                    # add term to ODE b
                                    odebexpr = add_diffusive_term(bode.loc[obname].at['expression'], thisdiffconst, f'Values[{oaname}]', f'Values[{obname}]')
                                    # we update oaname here because it is not affected in the other type of connections
                                    set_parameters(oaname, exact=True, expression=odeaexpr, model=newmodel)
                            elif( linktype=='s' ):
                                thissyntaur = syntaur
                                thissyntaud = syntaud
                                thissynvsyn = synvsyn
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need a few parameters per reaction
                                    (level, dist) = args.cnoise
                                    thissyntaur = f'tau_r_{ode}_synapse_{suffa}-{suffb}'
                                    thissyntaud = f'tau_d_{ode}_synapse_{suffa}-{suffb}'
                                    thissynvsyn = f'Vsyn_{ode}_synapse_{suffa}-{suffb}'
                                    r = addnoise(taurinit,float(level),dist)
                                    d = addnoise(taudinit,float(level),dist)
                                    v = addnoise(vsyninit,float(level),dist)
                                    add_parameter(name=thissyntaur, initial_value=r, model=newmodel)
                                    add_parameter(name=thissyntaud, initial_value=d, model=newmodel)
                                    add_parameter(name=thissynvsyn, initial_value=v, model=newmodel)
                                # add a new ODE to represent the proportion of bound post-synaptic receptor
                                brname = f'br_{ode}_{suffa},{suffb}'
                                brexp = f'( 1 / Values[{thissyntaur}] - 1 / Values[{thissyntaud}] ) * ( 1 - Values[{brname}] ) / ( 1 + exp( Values[{synv0}] - Values[{oaname}] ) ) -  Values[{brname}] / Values[{thissyntaud}]'
                                add_parameter(brname, type='ode', expression=brexp, initial_value=0, model=newmodel)
                                # add a synaptic maximum conductance parameter
                                syngc = f'g_c_{ode}_{suffa},{suffb}_synapse'
                                # link it the master g...
                                if( linkg ):
                                    gexpr = f'Values[{linkedsyng}]'
                                    add_parameter(syngc, type='assignment', expression=gexpr, model=newmodel)
                                # ...or otherwise set it to its own value
                                else:
                                    if( args.cnoise ):
                                        gc = addnoise(gcinit,float(level),dist)
                                    else:
                                        gc = gcinit
                                    add_parameter(syngc, type='fixed', initial_value=gc,
                                    model=newmodel)
                                # add a term to the postsynaptic ode corresponding to the voltage that is proportional to the bound receptor
                                odebexpr = bode.loc[obname].at['expression'] + f' + Values[{syngc}] * Values[{brname}] * ( Values[{thissynvsyn}] - Values[{oaname}] )'
                            # obname is affected by all types of connection so we can only update it here, after if/elif statements
                            set_parameters(obname, exact=True, expression=odebexpr, model=newmodel)
                elif( dim == 2 ):
                    if(args.add_medium and linktype=='d'):
                        # we first need to add the global quantity for the medium
                        mediumode = f'{ode}_medium'
                        medexp = ''
                        add_parameter(mediumode, type='ode', expression=medexp, initial_value=mparams.loc[ode].at['initial_value'], model=newmodel)
                    for r in range(gridr):
                        for c in range(gridc):
                            suffa = f'{r+1},{c+1}'
                            oaname = f'{ode}_{suffa}'
                            if(args.add_medium and linktype=='d'):
                                # get ode to get expression
                                aode = get_parameters(oaname, exact=True, model=newmodel)
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-medium'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                # add term to medium ODE
                                medexp = add_diffusive_term(medexp, thisdiffconst,f'Values[{oaname}]',f'Values[{mediumode}]')
                                # add term to target ODE
                                odeaexpr = add_diffusive_term(aode.loc[oaname].at['expression'], thisdiffconst,f'Values[{mediumode}]', f'Values[{oaname}]')
                                set_parameters(mediumode, exact=True, type='ode', expression=medexp, model=newmodel)
                                set_parameters(oaname, exact=True, type='ode', expression=odeaexpr, model=newmodel)
                            if( c+1 < gridc ):
                                suffb = f'{r+1},{c+2}'
                                obname = f'{ode}_{suffb}'
                                # get ode to get expression
                                aode = get_parameters(oaname, exact=True, model=newmodel)
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                # add term to target ODE
                                odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Values[{obname}]', f'Values[{oaname}]')
                                set_parameters(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                # get ode to get expression
                                bode = get_parameters(obname, exact=True, model=newmodel)
                                # add term to target ODE
                                odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Values[{oaname}]', f'Values[{obname}]')
                                set_parameters(obname, exact=True, expression=odebexpr, model=newmodel)
                            if( r+1 < gridr ):
                                suffb = f'{r+2},{c+1}'
                                obname = f'{ode}_{suffb}'
                                # get ode to get expression
                                aode = get_parameters(oaname, exact=True, model=newmodel)
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                # add term to target ODE
                                odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Values[{obname}]', f'Values[{oaname}]')
                                set_parameters(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                # get ode to get expression
                                bode = get_parameters(obname, exact=True, model=newmodel)
                                # add term to target ODE
                                odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Values[{oaname}]', f'Values[{obname}]')
                                set_parameters(obname, exact=True, expression=odebexpr, model=newmodel)
                elif( dim == 3 ):
                    if(args.add_medium and linktype=='d'):
                        # we first need to add the global quantity for the medium
                        mediumode = f'{ode}_medium'
                        medexp = ''
                        add_parameter(mediumode, type='ode', expression=medexp, initial_value=mparams.loc[ode].at['initial_value'], model=newmodel)
                    for r in range(gridr):
                        for c in range(gridc):
                            for l in range(gridl):
                                suffa = f'{r+1},{c+1},{l+1}'
                                oaname = f'{ode}_{suffa}'
                                if(args.add_medium and linktype=='d'):
                                    # get ode to get expression
                                    aode = get_parameters(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-medium'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to medium ODE
                                    medexp = add_diffusive_term(medexp, thisdiffconst,f'Values[{oaname}]',f'Values[{mediumode}]')
                                    # add term to target ODE
                                    odeaexpr = add_diffusive_term(aode.loc[oaname].at['expression'], thisdiffconst,f'Values[{mediumode}]', f'Values[{oaname}]')
                                    set_parameters(mediumode, exact=True, type='ode', expression=medexp, model=newmodel)
                                    set_parameters(oaname, exact=True, type='ode', expression=odeaexpr, model=newmodel)
                                if( c+1 < gridc ):
                                    suffb = f'{r+1},{c+2},{l+1}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_parameters(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to target ODE
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Values[{obname}]', f'Values[{oaname}]')
                                    set_parameters(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_parameters(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Values[{oaname}]', f'Values[{obname}]')
                                    set_parameters(obname, exact=True, expression=odebexpr, model=newmodel)
                                if( r+1 < gridr ):
                                    suffb = f'{r+2},{c+1},{l+1}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_parameters(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to target ODE
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Values[{obname}]', f'Values[{oaname}]')
                                    set_parameters(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_parameters(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Values[{oaname}]', f'Values[{obname}]')
                                    set_parameters(obname, exact=True, expression=odebexpr, model=newmodel)
                                if( l+1 < gridl ):
                                    suffb = f'{r+1},{c+1},{l+2}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_parameters(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to target ODE
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Values[{obname}]', f'Values[{oaname}]')
                                    set_parameters(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_parameters(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Values[{oaname}]', f'Values[{obname}]')
                                    set_parameters(obname, exact=True, expression=odebexpr, model=newmodel)
                else:
                    # insanity error
                    print('ERROR: this should not happen, dim is not 1,2 or 3. Contact developer!')
                    exit()
            # check if the ODE is a species
            elif( (mspecs is not None) and (ode in mspecs.index) ):
                if( mspecs.loc[ode].at['type'] != 'ode' ):
                    print(f'ERROR: {ode} is a species but it is not of type ODE')
                    continue
                # add medium unit if needed
                if(args.add_medium and linktype=='d'):
                    # add the new medium species
                    mediumode = f'{ode}_medium'
                    medexp = ''
                    add_species(mediumode, compartment_name=medium_name, type='ode', expression=medexp, initial_concentration=mspecs.loc[ode].at['initial_concentration'], model=newmodel)
                if( dim == 1):
                    # add coupling between unit and the medium
                    if(args.add_medium and linktype=='d'):
                        for r in range(gridr):
                            # name of ode in this unit
                            suffa = f'{r+1}'
                            oname = f'{ode}_{suffa}'
                            thisdiffconst = diffconst
                            if( args.cnoise ):
                                # if we add noise to couplings, then we need 1 parameter per reaction
                                (level, dist) = args.cnoise
                                thisdiffconst = f'k_{ode}_coupling_{suffa}-medium'
                                v = addnoise(coupleconst,float(level),dist)
                                add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                            # add term to medium ODE
                            medexp = add_diffusive_term( medexp, thisdiffconst, f'[{oname}]', f'[{mediumode}]')
                            set_species(mediumode, exact=True, expression=medexp, model=newmodel)
                            # get target ode to get expression
                            tode = get_species(oname, exact=True, model=newmodel)
                            # add term to target ODE
                            odexpr = add_diffusive_term( tode.loc[oname].at['expression'], thisdiffconst, f'[{mediumode}]', f'[{oname}]' )
                            set_species(oname, exact=True, expression=odexpr, model=newmodel)
                    if( args.network ):
                        for link in links:
                            if( link[0] == link[1]):
                                print(f' Warning: diffusive coupling onto the same unit not allowed, ignoring {link[0]} -> {link[1]}')
                                continue
                            suffa = f'{link[0]}'
                            oaname = f'{ode}_{suffa}'
                            suffb = f'{link[1]}'
                            obname = f'{ode}_{suffb}'
                            # get ode expressions from a and b
                            aode = get_species(oaname, exact=True, model=newmodel)
                            bode = get_species(obname, exact=True, model=newmodel)
                            if( linktype=='d' ):
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                if digraph:
                                    # add term to ODE a
                                    odeaexpr = add_irr_diffusive_term(aode.loc[oaname].at['expression'], thisdiffconst,'-',f'[{oaname}]')
                                    # add term to ODE b
                                    odebexpr = add_irr_diffusive_term(bode.loc[obname].at['expression'], thisdiffconst,'+',f'[{oaname}]')
                                else:
                                    # add term to ODE a
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{obname}]', f'[{oaname}]')
                                    # add term to ODE b
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'[{oaname}]', f'[{obname}]')
                                # we update oaname species here since it won't be modified in other types of connection
                                set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                            elif( linktype=='s' ):
                                thissyntaur = syntaur
                                thissyntaud = syntaud
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 tau parameter per reaction
                                    (level, dist) = args.cnoise
                                    thissyntaur = f'tau_r_{ode}_synapse_{suffa}-{suffb}'
                                    thissyntaud = f'tau_d_{ode}_synapse_{suffa}-{suffb}'
                                    r = addnoise(taurinit,float(level),dist)
                                    d = addnoise(taudinit,float(level),dist)
                                    add_parameter(name=thissyntaur, initial_value=r, model=newmodel)
                                    add_parameter(name=thissyntaud, initial_value=d, model=newmodel)
                                # add a new ODE to represent the proportion of bound post-synaptic receptor
                                brname = f'br_{ode}_{suffa},{suffb}'
                                brexp = f'( 1 / Values[{thissyntaur}] - 1 / Values[{thissyntaud}] ) * ( 1 - Values[{brname}] ) / ( 1 + exp( Values[{synv0}] - [{oaname}] ) ) -  Values[{brname}] / Values[{thissyntaud}]'
                                add_parameter(brname, type='ode', expression=brexp, initial_value=0.5, model=newmodel)
                                # add a synaptic maximum conductance parameter
                                syngc = f'g_c_{ode}_{suffa},{suffb}_synapse'
                                if( args.cnoise ):
                                    gc = addnoise(gcinit,float(level),dist)
                                else:
                                    gc = gcinit
                                add_parameter(syngc, type='fixed', initial_value=gc, model=newmodel)
                                # add a term to the postsynaptic ode corresponding to the voltage that is proportional to the bound receptor
                                odebexpr = bode.loc[obname].at['expression'] + f' + Values[{syngc}] * Values[{brname}] * ( Values[{synvsyn}] - [{oaname}] )'
                            # update species obname here as it was changed by all types
                            set_species(obname, exact=True, expression=odebexpr, model=newmodel)
                elif( dim == 2 ):
                    for r in range(gridr):
                        for c in range(gridc):
                            suffa = f'{r+1},{c+1}'
                            oaname = f'{ode}_{suffa}'
                            if(args.add_medium and linktype=='d'):
                                # add term to medium ODE
                                medexp = add_diffusive_term( medexp, diffconst, f'[{oaname}]', f'[{mediumode}]')
                                set_species(mediumode, exact=True, expression=medexp, model=newmodel)
                                # get ode to get expression
                                aode = get_species(oaname, exact=True, model=newmodel)
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-medium'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                # add term to medium ODE
                                medexp = add_diffusive_term( medexp, thisdiffconst, f'[{oname}]', f'[{mediumode}]')
                                set_species(mediumode, exact=True, expression=medexp, model=newmodel)
                                # add term to target ODE
                                odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{mediumode}]', f'[{oaname}]' )
                                set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                            if( c+1 < gridc ):
                                suffb = f'{r+1},{c+2}'
                                obname = f'{ode}_{suffb}'
                                # get ode to get expression
                                aode = get_species(oaname, exact=True, model=newmodel)
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                # add term to target ODE
                                odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{obname}]', f'[{oaname}]')
                                set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                # get ode to get expression
                                bode = get_species(obname, exact=True, model=newmodel)
                                # add term to target ODE
                                odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'[{oaname}]', f'[{obname}]')
                                set_species(obname, exact=True, expression=odebexpr, model=newmodel)
                            if( r+1 < gridr ):
                                suffb = f'{r+2},{c+1}'
                                obname = f'{ode}_{suffb}'
                                # get ode to get expression
                                aode = get_species(oaname, exact=True, model=newmodel)
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                # add term to target ODE
                                odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{obname}]', f'[{oaname}]')
                                set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                # get ode to get expression
                                bode = get_species(obname, exact=True, model=newmodel)
                                # add term to target ODE
                                odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'[{oaname}]', f'[{obname}]')
                                set_species(obname, exact=True, expression=odebexpr, model=newmodel)
                elif( dim == 3 ):
                    for r in range(gridr):
                        for c in range(gridc):
                            for l in range(gridl):
                                suffa = f'{r+1},{c+1},{l+1}'
                                oaname = f'{ode}_{suffa}'
                                if(args.add_medium and linktype=='d'):
                                    # add term to medium ODE
                                    medexp = add_diffusive_term( medexp, diffconst, f'[{oname}]', f'[{mediumode}]')
                                    set_species(mediumode, exact=True, expression=medexp, model=newmodel)
                                    # get ode to get expression
                                    aode = get_species(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-medium'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to medium ODE
                                    medexp = add_diffusive_term( medexp, thisdiffconst, f'[{oname}]', f'[{mediumode}]')
                                    set_species(mediumode, exact=True, expression=medexp, model=newmodel)
                                    # add term to target ODE
                                    odeaxpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{mediumode}]', f'[{oname}]' )
                                    set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                if( c+1 < gridc ):
                                    suffb = f'{r+1},{c+2},{l+1}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_species(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to target ODE
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{obname}]', f'[{oaname}]')
                                    set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_species(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'[{oaname}]', f'[{obname}]')
                                    set_species(obname, exact=True, expression=odebexpr, model=newmodel)
                                if( r+1 < gridr ):
                                    suffb = f'{r+2},{c+1},{l+1}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_species(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to target ODE
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{obname}]', f'[{oaname}]')
                                    set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_species(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'[{oaname}]', f'[{obname}]')
                                    set_species(obname, exact=True, expression=odebexpr, model=newmodel)
                                if( l+1 < gridl ):
                                    suffb = f'{r+1},{c+1},{l+2}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_species(oaname, exact=True, model=newmodel)
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    # add term to target ODE
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'[{obname}]', f'[{oaname}]')
                                    set_species(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_species(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'[{oaname}]', f'[{obname}]')
                                    set_species(obname, exact=True, expression=odebexpr, model=newmodel)
                else:
                    # insanity error
                    print('ERROR: this should not happen, dim is not 1,2 or 3. Contact developer!')
                    exit()

            # check if the ODE is a compartment
            elif( (mcomps is not None) and (ode in mcomps.index) ):
                if( ignc ):
                    print(f' Warning: {ode} is a compartment but ignore_compartments is set, nothing done')
                else:
                    if( mcomps.loc[ode].at['type'] != 'ode' ):
                        print(f'ERROR: {ode} is a compartment but it is not of type ODE')
                        exit()
                    if(linktype=='s'):
                        print(f'ERROR: {ode} is a compartment ODE, but compartments cannot have synaptic links')
                        exit()
                    if( dim == 1):
                        # for compartments we do nothing about the medium...
                        # let's just work on the network
                        if( args.network ):
                            for link in links:
                                if( link[0] == link[1]):
                                    print(f' Warning: diffusive coupling onto the same unit not allowed, ignoring {link[0]} -> {link[1]}')
                                    continue
                                suffa = f'{link[0]}'
                                oaname = f'{ode}_{suffa}'
                                suffb = f'{link[1]}'
                                obname = f'{ode}_{suffb}'
                                # get ode a to get expression
                                aode = get_compartments(oaname, exact=True, model=newmodel)
                                # get ode b to get expression
                                bode = get_compartments(obname, exact=True, model=newmodel)
                                thisdiffconst = diffconst
                                if( args.cnoise ):
                                    # if we add noise to couplings, then we need 1 parameter per reaction
                                    (level, dist) = args.cnoise
                                    thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                    v = addnoise(coupleconst,float(level),dist)
                                    add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                if digraph:
                                    # add term to ODE a
                                    odeaexpr = add_irr_diffusive_term(aode.loc[oaname].at['expression'], thisdiffconst,'-',f'Compartments[{oaname}].Volume')
                                    # add term to ODE b
                                    odebexpr = add_irr_diffusive_term(bode.loc[obname].at['expression'], thisdiffconst,'+',f'Compartments[{oaname}].Volume')
                                else:
                                    # add term to ODE a
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Compartments[{obname}].Volume', f'Compartments[{oaname}].Volume')
                                    # add term to ODE b
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Compartments[{oaname}].Volume', f'Compartments[{obname}].Volume')
                                set_compartment(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                set_compartment(obname, exact=True, expression=odebexpr, model=newmodel)
                    elif( dim == 2 ):
                        for r in range(gridr):
                            for c in range(gridc):
                                suffa = f'{r+1},{c+1}'
                                oaname = f'{ode}_{suffa}'
                                if( c+1 < gridc ):
                                    suffb = f'{r+1},{c+2}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_compartments(oaname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Compartments[{obname}].Volume', f'Compartments[{oaname}].Volume')
                                    set_compartment(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_compartments(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Compartments[{oaname}].Volume', f'Compartments[{obname}].Volume')
                                    set_compartment(obname, exact=True, expression=odebexpr, model=newmodel)
                                if( r+1 < gridr ):
                                    suffb = f'{r+2},{c+1}'
                                    obname = f'{ode}_{suffb}'
                                    # get ode to get expression
                                    aode = get_compartments(oaname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    thisdiffconst = diffconst
                                    if( args.cnoise ):
                                        # if we add noise to couplings, then we need 1 parameter per reaction
                                        (level, dist) = args.cnoise
                                        thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                        v = addnoise(coupleconst,float(level),dist)
                                        add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                    odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Compartments[{obname}].Volume', f'Compartments[{oaname}].Volume')
                                    set_compartment(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                    # get ode to get expression
                                    bode = get_compartments(obname, exact=True, model=newmodel)
                                    # add term to target ODE
                                    odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Compartments[{oaname}].Volume', f'Compartments[{obname}].Volume')
                                    set_compartment(obname, exact=True, expression=odebexpr, model=newmodel)
                    elif( dim == 3 ):
                        for r in range(gridr):
                            for c in range(gridc):
                                for l in range(gridl):
                                    suffa = f'{r+1},{c+1},{l+1}'
                                    oaname = f'{ode}_{suffa}'
                                    if( c+1 < gridc ):
                                        suffb = f'{r+1},{c+2},{l+1}'
                                        obname = f'{ode}_{suffb}'
                                        # get ode to get expression
                                        aode = get_compartments(oaname, exact=True, model=newmodel)
                                        # add term to target ODE
                                        thisdiffconst = diffconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 parameter per reaction
                                            (level, dist) = args.cnoise
                                            thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                            v = addnoise(coupleconst,float(level),dist)
                                            add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                        odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Compartments[{obname}].Volume', f'Compartments[{oaname}].Volume')
                                        set_compartment(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                        # get ode to get expression
                                        bode = get_compartments(obname, exact=True, model=newmodel)
                                        # add term to target ODE
                                        odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Compartments[{oaname}].Volume', f'Compartments[{obname}].Volume')
                                        set_compartment(obname, exact=True, expression=odebexpr, model=newmodel)
                                    if( r+1 < gridr ):
                                        suffb = f'{r+2},{c+1},{l+1}'
                                        obname = f'{ode}_{suffb}'
                                        # get ode to get expression
                                        aode = get_compartments(oaname, exact=True, model=newmodel)
                                        # add term to target ODE
                                        thisdiffconst = diffconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 parameter per reaction
                                            (level, dist) = args.cnoise
                                            thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                            v = addnoise(coupleconst,float(level),dist)
                                            add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                        odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Compartments[{obname}].Volume', f'Compartments[{oaname}].Volume')
                                        set_compartment(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                        # get ode to get expression
                                        bode = get_compartments(obname, exact=True, model=newmodel)
                                        # add term to target ODE
                                        odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Compartments[{oaname}].Volume', f'Compartments[{obname}].Volume')
                                        set_compartment(obname, exact=True, expression=odebexpr, model=newmodel)
                                    if( l+1 < gridl ):
                                        suffb = f'{r+1},{c+1},{l+2}'
                                        obname = f'{ode}_{suffb}'
                                        # get ode to get expression
                                        aode = get_compartments(oaname, exact=True, model=newmodel)
                                        # add term to target ODE
                                        thisdiffconst = diffconst
                                        if( args.cnoise ):
                                            # if we add noise to couplings, then we need 1 parameter per reaction
                                            (level, dist) = args.cnoise
                                            thisdiffconst = f'k_{ode}_coupling_{suffa}-{suffb}'
                                            v = addnoise(coupleconst,float(level),dist)
                                            add_parameter(name=thisdiffconst, initial_value=v, model=newmodel)
                                        odeaexpr = add_diffusive_term( aode.loc[oaname].at['expression'], thisdiffconst, f'Compartments[{obname}].Volume', f'Compartments[{oaname}].Volume')
                                        set_compartment(oaname, exact=True, expression=odeaexpr, model=newmodel)
                                        # get ode to get expression
                                        bode = get_compartments(obname, exact=True, model=newmodel)
                                        # add term to target ODE
                                        odebexpr = add_diffusive_term( bode.loc[obname].at['expression'], thisdiffconst, f'Compartments[{oaname}].Volume', f'Compartments[{obname}].Volume')
                                        set_compartment(obname, exact=True, expression=odebexpr, model=newmodel)
                    else:
                        # insanity error
                        print('ERROR: this should not happen, dim is not 1,2 or 3. Contact developer!')
                        exit()
            # not an entity in this model
            else:
                print(f'ERROR: {ode} is not a valid model entity')
                exit()

    # regulated synthesis (GRNs)
    if( args.grn ):
        # GRNs only work with dimension 1 and a network graph
        if (dim != 1) or (not args.network ):
            # error
            print( f'ERROR: regulatory synthesis reactions can only be created with dimension 1 and through a network (use option -n)' )
            exit()
        # if we are not ignoring compartments issue warning suggesting we should
        if not args.ignore_compartments:
            print(f' Warning: option --ignore-compartments is often desirable for building regulatory networks')
        # create all rate-laws for regulatory synthesis
        for k in indegrees:
            # add a rate law with this number of in-degrees
            add_regulatory_ratelaw(k)
        # iterate over all species with regulatory synthesis
        for sp in args.grn:
            # check that the species exists
            if( (seednspecs>0) and (sp in mspecs.index) ):
                # create global quantities for parameters
                GRN_V_const = f'V_synth_{sp}'
                GRN_a_const = f'a_synth_{sp}'
                GRN_h_const = f'h_synth_{sp}'
                if( not args.cnoise ):
                    add_parameter(name=GRN_V_const, initial_value=grnV, model=newmodel)
                    add_parameter(name=GRN_a_const, initial_value=grna, model=newmodel)
                    add_parameter(name=GRN_h_const, initial_value=grnh, model=newmodel)
                # create all regulatory synthesis reactions
                for target in reglinks:
                    rname = f'synthesis {sp}_{target}'
                    rscheme = f'  -> {sp}_{target};'
                    nregs = regulated[target]
                    thisV = GRN_V_const
                    thisa = GRN_a_const
                    thish = GRN_h_const
                    if( args.cnoise ):
                        # if we add noise to regulations, then we need 1 V per target
                        (level, dist) = args.cnoise
                        thisV = f'V_synth_{sp}_{target}'
                        val = -1
                        # make sure that val is positive (for V parameter)
                        while val < 0:
                            val = addnoise(grnV,float(level),dist)
                        add_parameter(name=thisV, initial_value=val, model=newmodel)
                    rmap = {'V': thisV}
                    i=0
                    for source in reglinks[target]:
                        i = i+1
                        rscheme = rscheme + f' {sp}_{source}'
                        rmap[f'M{i}'] = f'{sp}_{source}'
                        if( args.cnoise ):
                            # if we add noise to regulations, then we need 1 a and h per source
                            (level, dist) = args.cnoise
                            thisa = f'a_synth_{sp}_{source}-{target}'
                            val = addnoise(grna,float(level),dist)
                            # enforce boundaries: a should be in [-1,1]
                            if val < -1:
                                val = -1
                            if val > 1:
                                val = 1
                            add_parameter(name=thisa, initial_value=val, model=newmodel)
                            thish = f'h_synth_{sp}_{source}-{target}'
                            val = addnoise(grnh,float(level),dist)
                            val = float(int(val))
                            # enforce boundaries: h should be in [1,10]
                            if val < 1:
                                val = 1
                            if val > 10:
                                val = 10
                            add_parameter(name=thish, initial_value=val, model=newmodel)
                        rmap[f'h{i}'] = thish
                        rmap[f'a{i}'] = thisa
                        # add reaction from source  to link[1] with rate law regulated_by{regulated[link1]}
                    add_reaction(model=newmodel, name=rname, scheme=rscheme, mapping=rmap, function=f'regulated_by_{nregs}' )
            else:
                # error
                print( f'ERROR: Species {sp} does not exist in the model, no regulatory synthesis reactions added' )
                exit()

    #####
    # 13. set task parameters
    #####

    if( not args.ignore_tasks):
        # time course
        tc = get_task_settings('Time-Course', basic_only=False, model=seedmodel)
        # if report is not the default one, clear it
        if( tc['report']['report_definition'] != 'Time-Course' ):
            tc['report']['report_definition'] = 'Time-Course'
            tc['report']['filename'] = ""
        set_task_settings('Time-Course', {'scheduled': tc['scheduled'], 'problem': tc['problem'], 'method': tc['method'], 'report': tc['report']}, model=newmodel)
        # steady state
        ss = get_task_settings('Steady-State', basic_only=False, model=seedmodel)
        # if report is not the default one, clear it
        if( ss['report']['report_definition'] != 'Steady-State' ):
            ss['report']['report_definition'] = 'Steady-State'
            ss['report']['filename'] = ""
        set_task_settings('Steady-State', {'scheduled': ss['scheduled'], 'update_model': ss['update_model'], 'problem': ss['problem'], 'method': ss['method'], 'report': ss['report']},model=newmodel)

        # MCA
        mca = get_task_settings('Metabolic Control Analysis', basic_only=False, model=seedmodel)
        # if report is not the default one, clear it
        if( mca['report']['report_definition'] != 'Metabolic Control Analysis' ):
            mca['report']['report_definition'] = 'Metabolic Control Analysis'
            mca['report']['filename'] = ""
        set_task_settings('Metabolic Control Analysis', {'scheduled': mca['scheduled'], 'update_model': mca['update_model'], 'problem': mca['problem'], 'method': mca['method'], 'report': mca['report']},model=newmodel)

        # Lyapunov Exponents
        le = get_task_settings('Lyapunov Exponents', basic_only=False, model=seedmodel)
        # if report is not the default one, clear it
        if( le['report']['report_definition'] != 'Lyapunov Exponents' ):
            le['report']['report_definition'] = 'Lyapunov Exponents'
            le['report']['filename'] = ""
        set_task_settings('Lyapunov Exponents', {'scheduled': le['scheduled'], 'update_model': le['update_model'], 'problem': le['problem'], 'method': le['method'], 'report': le['report']},model=newmodel)

        # Time Scale Separation Analysis
        tsa = get_task_settings('Time Scale Separation Analysis', basic_only=False, model=seedmodel)
        # if report is not the default one, clear it
        if( tsa['report']['report_definition'] != 'Time Scale Separation Analysis' ):
            tsa['report']['report_definition'] = 'Time Scale Separation Analysis'
            tsa['report']['filename'] = ""
        set_task_settings('Time Scale Separation Analysis', {'scheduled': tsa['scheduled'], 'update_model': tsa['update_model'], 'problem': tsa['problem'], 'method': tsa['method'], 'report': tsa['report']},model=newmodel)

        # Cross section
        cs = get_task_settings('Cross Section', basic_only=False, model=seedmodel)
        # there is no standard report for cross section...
        if( cs['problem']['SingleVariable'] != ''):
            newv = fix_expression(cs['problem']['SingleVariable'], apdx1, ignc)
            print(f' Warning: the cross section task was updated to use {newv} as variable.')
            cs['problem']['SingleVariable'] = newv
            set_task_settings('Cross Section', {'scheduled': cs['scheduled'], 'update_model': cs['update_model'], 'problem': cs['problem'], 'method': cs['method']},model=newmodel)

        # Linear Noise Approximation
        lna = get_task_settings('Linear Noise Approximation', basic_only=False, model=seedmodel)
        # if report is not the default one, clear it
        if( lna['report']['report_definition'] != 'Linear Noise Approximation' ):
            lna['report']['report_definition'] = 'Linear Noise Approximation'
            lna['report']['filename'] = ""
        set_task_settings('Linear Noise Approximation', {'scheduled': lna['scheduled'], 'update_model': lna['update_model'], 'problem': lna['problem'], 'method': lna['method'], 'report': lna['report']},model=newmodel)

        # Sensitivities
        sen = get_sensitivity_settings(model=seedmodel)
        # if report is not the default one, clear it
        if( sen['report']['report_definition'] != 'Sensitivities' ):
            sen['report']['report_definition'] = 'Sensitivities'
            sen['report']['filename'] = ""
        seff = fix_expression(sen['effect'],apdx1, ignc)
        scau = fix_expression(sen['cause'],apdx1, ignc)
        ssec = fix_expression(sen['secondary_cause'],apdx1, ignc)
        if( (seff != sen['effect']) or (scau != sen['cause']) or (ssec != sen['secondary_cause']) ):
            print(f' Warning: sensitivies task is now using items of unit {apdx1}")')
            sen['effect'] = seff
            sen['cause'] = scau
            sen['secondary_cause'] = ssec
        set_sensitivity_settings(sen, model=newmodel)

        # Parameter scan
        ps = get_task_settings('Scan', basic_only=False, model=seedmodel)
        set_task_settings('Scan', {'scheduled': ps['scheduled'], 'update_model': ps['update_model'], 'problem': ps['problem'], 'method': ps['method']},model=newmodel)

        # we got the scanitems much further up... (due to a bug in COPASI/BasiCO)
        # when there are scan or random sampling items, we convert them to be those of the first unit
        srw = False
        for sit in scanitems:
            if( sit['type']=='parameter_set' ):
                print(f' Warning: a scan of parameter sets exists in the original model but was not included in the new model.')
            else:
                if( sit['type']=='scan' ):
                    newit = fix_expression(sit['item'], apdx1, ignc)
                    srw = True
                    add_scan_item(model=newmodel, type=sit['type'], num_steps=sit['num_steps'], item=newit, log=sit['log'], min=sit['min'], max=sit['max'], use_values=sit['use_values'], values=sit['values'] )
                else:
                    if( sit['type']=='random' ):
                        newit = fix_expression(sit['item'], apdx1, ignc)
                        srw = True
                        add_scan_item(model=newmodel, type=sit['type'], num_steps=sit['num_steps'], item=newit, log=sit['log'], min=sit['min'], max=sit['max'], distribution=sit['distribution'])
                    else:
                        if( sit['type']=='repeat' ):
                            add_scan_item(model=newmodel, type=sit['type'], num_steps=sit['num_steps'])
                        else:
                            tp = sit['type']
                            print(f' Warning: This scan task includes an unknonw type {tp}, likely from a new version of COPASI. Please file an issue on Github.')
        if( srw ): print(' Warning: in Parameter scan task the scanned or sampled items were converted to those of the first unit only.')

        # Optimization
        # we translate the objective function and parameters to the first unit
        nopt = get_opt_settings(model=seedmodel)
        # if report is not the default one, clear it
        if( nopt['report']['report_definition'] != 'Optimization' ):
            nopt['report']['report_definition'] = 'Optimization'
            nopt['report']['filename'] = ""
        set_opt_settings(nopt, model=newmodel)
        if( nopt['expression'] ):
            nopt['expression'] = fix_expression(nopt['expression'], apdx1, ignc)
            set_opt_settings(nopt,newmodel)
            ops = get_opt_parameters(model=seedmodel)
            if ops is not None:
                for p in ops.index:
                    # rename the CN
                    ops.loc[p, 'cn'] = fix_expression(ops.loc[p].at['cn'] ,apdx1, ignc)
                    # rename the index
                    newp = fix_expression(p,apdx1, ignc)
                    ops.rename(index={p: newp}, inplace=True)
                set_opt_parameters(ops, model=newmodel)
            cst = get_opt_constraints(model=seedmodel)
            if cst is not None:
                # deal with constraints
                for p in cst.index:
                    # rename the CN
                    cst.loc[p, 'cn'] = fix_expression(cst.loc[p].at['cn'] ,apdx1, ignc)
                    # rename the index
                    newp = fix_expression(p,apdx1, ignc)
                    cst.rename(index={p: newp}, inplace=True)
                set_opt_constraints(cst, model=newmodel)
            print(' Warning: in Optimization task the objective function and the search parameters were converted to those of the first unit only.')

        # We won't do Parameter estimation but need to issue a warning if it was set
        exps = get_experiment_filenames(model=seedmodel)
        if( len(exps)>0 ):
            print(' Warning: Parameter Estimation task settings were not copied to the new model.')

        # We won't do Time Course Sensitivities but need to issue a warning if it was set
        tcs = get_task_settings( 'Time-Course Sensitivities', basic_only=False, model=seedmodel)
        if( tcs['scheduled'] ):
            print(' Warning: Time Course Sensitivities task settings were not copied to the new model.')

    #TODO: what to do with plots?


    #####
    # 14. save model
    #####

    # get the base of the output model filename
    base,ext = os.path.splitext(newfilename)

    # save a COPASI file if the extension was .cps
    if( (ext.lower() == '.cps') and not args.sbml ):
        save_model(filename=newfilename, model=newmodel)
    # otherwise save an SBML file
    else:
        save_model(filename=newfilename, type='sbml', sbml_level=sbmll, sbml_version=sbmlv, model=newmodel)
    if( not args.quiet ):
        print(f"created new model {newfilename} with {desc} of {seedmodelfile}\n")

if __name__ == '__main__':
    main()
