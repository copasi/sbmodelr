# *sbmodelr*  User Manual

## Summary
*sbmodelr* takes as input an SBML or COPASI file describing a model. It then creates a new model that replicates the original in several copies ('units') where each one may interact with other ones. The topology of the connections can be arbitrary (described in a DOT network file), a 2D rectangular matrix, or a 3D cuboid array. Several options are available to randomize parameter values, supress creation of compartments, etc.

For a full set of options run `sbmodelr --help` on the command line.

### Examples
The rest of this document describes the many options that are available in *sbmodelr*, but specific examples of usage are provided in the `examples` folder, which illustrate real-world usage of this tool, mostly replicating existing modeling papers. The examples are:

 - [Cells in a medium](examples/Cells_in_medium#cells-in-medium): single cell organisms with species transported to the medium
 - [Row of cells with gap junctions](examples/Row_of_cells_gap_junctions#row-of-cells-with-gap-junctions): cells connected end-to-end by gap junctions
 - [Array of oscillating cells](examples/Array_of_oscillating_cells#array-of-cells): a square array of cells in a medium with diffusion
 - [Neuron networks I](examples/Neuron_networks_I#neuron-networks-i): small(ish) networks using a Hodgkin-Huxley-type model
 - [Neuron networks II](examples/Neuron_networks_II#neuron-networks-ii): large networks using the Izhikevich model
 - [Gene regulatory networks](examples/Gene_Regulatory_Networks#gene-regulatory-networks): easy creation of large (or small) gene networks from a simple model file

## Required options

*sbmodler* requires at least two command line arguments: 1) a base model file, and 2) the number of units to replicate. The simplest command that can be issued is: `sbmodelr mybasemodel.cps 2`; this would create a new file called `mybasemodel_2.cps` with two units that are exact copies of the model in `mybasemodel.cps`.

### File Input and Output
The **input** base model is either encoded in an [SBML](https:sbml.org) file (up to L3v2) or a [COPASI](https://copasi.org) file (extension `.cps`). The output of *sbmodelr* will be in the same format as the base model file unless forced to be in a specific format (explained below).

The **output** of *sbmodelr* may be either another SBML or COPASI file with a full exapanded model or, alternatively, the files required to run an equivalent multiscale simulation using the [Vivarium](https://github.com/vivarium-collective/process-bigraph) framework.

By default the output file will be named after the input file with an appendix to its name reflecting the number of replicate units. To specify a different output filename use the option `-o filename` or `--output filename`, for example the command `sbmodelr -o newmodel.cps basemodel.cps 2` would put the new model in the file `newmodel.cps` (without this option it would be `basemodel_2.cps`).

 - To force the output file to be written in SBML format add the option `--sbml`; this option can take an argument specifying the level and version of SBML required (one of `l1v2`,`l2v3`,`l2v4`,`l2v5`,`l3v1`,`l3v2`).
 - To force the output file to be written in COPASI format you will need to explicitly name the output file with option `-o filename.cps` or `--output filename.cps`, ensuring that the filename ends with `.cps` extension.
 - To create a *Vivarium* multiscale simulation, use the option `--vivarium`; that will then trigger the output of a JSON file with the topology of the new (compositional) simulation, and a Python file encapsulating the *Vivarium* class required to process the simulation. The two then need to be processed with *Vivarium* (version 2, or *Process-Bigraph*).

### Number of replicates and their connectivity

*sbmodelr* creates more than one copy of the base model, organized as a 2D rectangular matrix, a 3D cuboid array, or a set of units with arbitrary connections. The number of units created is specified with 1, 2 or 3 numbers after the base model filename.

 - *2D Rectangular matrices* are created with two numbers specifying the rows and columns of the matrix. The command `sbmodelr basemodel.cps 3 5` creates 15 units organized as a matrix with 3 rows and 5 columns. If interactions are specified (see below) then they will be between a unit's 4 neighboring units (left, right, top and bottom for the "bulk" units, while only 3 other units at the edges, and only 2 other units in the corners).  Note that this topology *is not* toroidal.
 - *3D Cuboid arrays* are created with three numbers specifying the rows, columns, and levels of the 3D array. The command `sbmodelr basemodel.cps 3 5 2` creates 30 units organized as a 3D array with 3 rows, 5 columns, and 2 levels. If interactions are specified (see below) then they will be between a unit's 6 neighboring units (left, right, top, bottom, front and back for the "bulk" units, while only 4 other units at the edges, and only 3 other units in the corners).  Note that this topology *is not* toroidal.
 - *Sets of units with arbitrary connections* are created with just one number after the base filename, specifying the total number of replicate units desired. The command `sbmodelr basemodel.cps 10` creates 10 units that are replicates of the base model, and as such they would have no connections between them. Any number of connections can be specified through a network file following the [DOT format](https://graphviz.org/doc/info/lang.html) of the [GraphViz package](https://graphviz.org/). Network files specifying connections are passed on with the `-n mynetfile.dot` or `--network mynetfile.dot` options. In addition to having to conform to the DOT format, *sbmodelr* also requires that the node names be integer numbers to specify each unit. Thus if we specify 10 units, the network file can only have nodes named `1` to `10`. The command `sbmodelr -n mynetfile.dot basemodel.cps 10` creates 10 replicate units that are connected through a network described in file `mynetfile.dot`. As an example, that file's content could be simply this: `graph { 1 -- 2 -- 3 -- 4 -- 5 -- 6 -- 7 -- 8 -- 9 -- 10}`, which defines a string of connected units with units 1 and 10 at the ends. The DOT language allows for two types of networks: undirected (specified by the keyword `graph` and edges as `--`) or directed (specified by the keyword `digraph` and edges as `->`). The outcome of using one or other depends on the types of connections required and is discussed below.  The examples further illustrate several network files. Network files do not need to have the `.dot` extension (or `.gv` that is also common), they can have any name and extension as long as they conform to the DOT syntax and the node names are integers.

While *sbmodelr* has special provisions to create 2D and 3D arrays, *any* kind of topology can be specified through the use of appropriate network files.

## Connecting the units
Currently there are four types of connections between units that can be added to the output model. They can be 1) transport of species, 2) diffusive connection of explicit ODEs, 3) regulatory interactions on the synthesis of species, and 4) chemical synapses through the method of Destexhe *et al.* (1994).

At present, *Vivarium* multiscale simulations can only use transport of species. This limitation is expected to be removed in future versions.

### Transport of species

If the base model has species that one wants to allow being transported between units, this can be specified with the option `-t species` or `--transport species`, where `species` is the species name. This will create transport reactions that are governed by mass action kinetics, where the rate constant is the same in both directions. A more general way of specifying transport is to use the option `--Hill-transport species` which will create transport steps following Hill kinetics. If connectivity is specified through a network file (option `-n netfile`) then the transport reactions will be reversible if the file specifies a `graph` with edges represented with `--`, or irreversible if the file specifies a `digraph` with edges represented with `->`. More specifically the rate laws used are described in the table below.

| transport option                    | network type | rate law                        | notes                     |
| ----------------------------------- | ------------ | ------------------------------- | ------------------------- |
|`-t species` or `--transport species`| graph        | $$k \cdot ( species_i - species_j )$$ | 2D and 3D matrices        |
|`-t species` or `--transport species`| digraph      | $$k \cdot species_i$$               | not for 2D or 3D matrices |
|`--Hill-transport species`           | graph        | $$\frac{Vmax \cdot (species_i^h - species_j^h )}{Km^h + species_i^h + species_j^h }$$ | 2D and 3D matrices        |
|`--Hill-transport species`           | digraph      | $$\frac{Vmax \cdot species_i^h}{Km^h + species_i^h}$$ | not for 2D or 3D matrices |

where *k* is a transport rate constant, *Vmax* is a maximal rate of transport, *Km* is the concentration of *species<sub>i</sub>* (and *species<sub>j</sub>*) when the rate is half of *Vmax*, *h* is a Hill exponent, where if it is 1 the rate is hyperbolic (essentially the Michaelis-Menten equation), or if larger than 1 the rate is sigmoidal; *i* and *j* are the indices of the two units.

When adding transport steps the parameters indicated in the rate laws above will get default values and will be the same for all transport reactions. To specify a value different from the default use the options in the table below. These parameters can also be randomized like the parameters of the base model, see section on *Randomizing parameter values* for more information.

| parameter | default | option to set value                       |
| --------- | ------- | ----------------------------------------- |
| *k*       | 1.0     | `-k value` or `--transport-k value`       |
| *Km*      | 1.0     | `--transport-Km value`                    |
| *Vmax*    | 1.0     | `--transport-Vmax value`                  |
| *h*       | 1.0     | `--transport-h value`                     |

### Diffusive connection of explicit ODEs

The option to indicate that a variable should be connected by a diffusive interaction is `--d variable` or `--ode-diffusive variable`. The `variable` must be defined as an explicit `ode` type in the base model, or this option will generate an error.

This type of connection allows connecting units by variables that are explicit ODEs, such as *species*, *global quantities*, or *compartments* that are defined as type `ode` (not `fixed`, `assignment`  or `reactions`).  The "diffusive" interaction is mathematically the same as a mass-action transport reaction. If units *i* and *j* are connected then the diffusive interaction adds the following terms to the right-hand side (rhs) of the respective ODEs of these variables:

| variable             | network type | new term on rhs of ODE            |
| -------------------- |------------- | --------------------------------- |
| variable<sub>i</sub> |  graph       | $$+ c·(variable_j - variable_i)$$ |
| variable<sub>j</sub> |  graph       | $$+ c·(variable_i - variable_j)$$ |
|                      |              |                                   |
| variable<sub>i</sub> |  digraph     | $$- c·variable_i$$                |
| variable<sub>j</sub> |  digraph     | $$+ c·variable_i$$                |

where *c* is a diffusive rate constant with a default value of 1.0. To use a different value for *c* use the option `-c value` or `--coupling-constant value`. This parameter can also be randomized like the parameters of the base model, see section on *Randomizing parameter values* for more information.

If a network was specified as a (bi-directional) graph, or we are creating 2D or 3D arrays, the interaction is symetric and acts like diffusion, where the unit with the highest value "flows" into the variable with the lowest value until they become equal. If a network is defined as a digraph (directed graph), then there is only "flow" from one unit to the other (whatever direction was defined in the network file). The table above describes the terms that are added to the right-hand side of the ODEs in each case.

Diffusive interactions can be used, for example, in connecting species that are transported, or diffuse between two compartments (*i.e.* formally the same as transport, but somehow these species were defined as explicit `ode` and so no reactions can be added). They also serve to connect two variables that represent electric potentials, where the constant *c* is interpreted as an conductivity between the two units.

This option is presently not supported if the output is a *Vivarium* multiscale simulation.

### Regulatory interactions on the synthesis of species

This type of connection is useful to create regulatory networks, particularly gene regulatory networks. It takes a species in the base unit and it will add a synthesis reaction for that species; that synthesis reaction is then modified (inhibited/activated) by the same species in other units.

The new synthesis reaction uses a general type of rate law that is composed of a product of regulatory terms, one for each modifier (*i.e.* the other units that affect this one). The general form of this equation is:

$$V \cdot \prod_i \frac{ 1 + ( 1 + a_i ) \cdot M_i^{h_i}}{1 + M_i^{h_i} }$$

where the subscript *i* represents all the units affecting this one, *M<sub>i</sub>* is the concentration of the *i*-th modifier species, parameter *h<sup>i</sup>* is a Hill coefficient, and parameter *a<sub>i</sub>* is an activation/inhibition strength. *V* is a basal synthesis rate (the rate when all *M<sub>i</sub>*=0). Parameter *a* can take values from -1 to +1, where negative values make the corresponding modifier be an inhibitor (repressor), and positive values an activator (inducer); a value of zero makes the corresponding modifier have no effect. The Hill coefficient *h* can take integer values between 1 and 10, where 1 makes the rate hyperbolic, and larger values make it an increasingly steep sigmoidal. The figure below shows the behavior of a regulatory term at different values of *M<sub>i</sub>*, *a<sub>i</sub>* and *h<sup>i</sup>*.

![regulatory term](diagrams/regfunction.png)

This type of connection cannot be used with 2D or 3D arrays, only with an explicit network file (see above, option `-n`), and it must be a `digraph` (directed graph, where the edges are unidirectional and specified with `->`).

In most uses of this type of connection, you want the resulting units to be contained inside the original compartment of the base unit (*e.g.* all genes in the same cell), this can be achieved by not replicating the compartments by using the option `--ignore-compartments`.

The three parameters of the rate law above will be assigned default values (whcih will be the same for all regulatory terms). To specify different values than the default use the options in the table below.  These parameters can also be randomized like the parameters of the base model, see section on *Randomizing parameter values* for more information.

| parameter | default | option to set value                       |
| --------- | ------- | ----------------------------------------- |
| *V*       | 1.0     | `--grn-V value`                           |
| *a*       | 1.0     | `--grn-a value`                           |
| *h*       | 2       | `--grn-h value`                           |


This option is presently not supported if the output is a *Vivarium* multiscale simulation.

### Chemical synapses

This type of connection is intended for electrophysiological models of neurons where the membrane potential is an explicit variable and the neurons are connected through chemical synapses. The base model should be of the Hodgkin-Huxley type where the variable to connect represents membrane potential. The connection is specified using `--s variable` or `--ode-synaptic variable`, where `variable` is the name of the variable representing the membrane potential.

This type of connection can only be used with variables that are of type `ode` (*i.e.* explicit ODEs, called "rate rules" in SBML), and these can be either *species* or *global quantities*. The connections must be defined by a network file of type `digraph` (directed graph, where the edges are unidirectional and specified with `->`). This cannot be used with the 2D and 3D array topologies.

Chemical synapses involve the release of a neurotransmitter by a presynaptic neuron caused by an action potential. The neurotransmitter diffuses across the synapse, eventually binding to a receptor in the postsynaptic neuron. As the receptor binds the neurotransmitter, it triggers an action potential in the postsynaptic neuron. Chemical synapses have the properties of being directional and depending on the diffusion, binding, and release or degradation of the neurotransmitter. [Destexhe *et al.* (1994)](https://doi.org/10.1162/neco.1994.6.1.14) published a simple kinetic model of chemical synapses that still provides a realistic reproduction of the phenomenon and is adopted in *sbmodelr*. Under this approach, each synapse requires adding one extra model variable: an ODE representing the proportion of bound postsynaptic neurotransmitter receptor. In a chemical synapse of neuron_i to neuron_j, where the membrane potential of neuron_i is *V<sub>i</sub>* and of neuron_j is *V<sub>j</sub>*, the proportion of bound receptor at neuron_j, *br<sub>i,j</sub>* is governed by the following differential equation:

$$\frac{d br_{i,j}}{dt} = \frac{ \( \frac{1}{tau\_r} - \frac{1}{tau\_d} \) \cdot (1 - br_{i,j})}{1 + e^{V_0 - V\_i}} -\frac{br_{i,j}}{tau\_d}$$

then the differential equation for the membrane potential at the postsynaptic neuron (neuron_j) is expanded by one term to represent the potential caused by the bound receptor:

$$\frac{d V_j}{dt} = ... + g \cdot br_{i,j} \cdot (V_{syn} - V_i)$$

The parameters involved are: *tau<sub>r</sub>* a characteristic time for release of neurotransmitter from the presynaptic neuron, *tau<sub>d</sub>* a characteristic time for clearance of the neurotransmitter from the bound postsynaptic receptor, *V<sub>0</sub>* a reversal potential, *V<sub>syn</sub>* the reversal potential of the synapse, and *g* is a synaptic weight (because *sbmodler* simply adds a term to the membrane potential equation, this does not get divided by the membrace capacitance thus,parameter *g* here is effectively the conductance of the synapse divided by the capacitance of the postsynaptic neuron). Note that a synapse becomes inhibitory if the parameter *V<sub>syn</sub>* becomes sufficiently negative (*e.g.* -40 mV), where the default value (20 mV) is for an excitatory synapse.

The five parameters of the two differential equations above will be assigned default values. To specify different values use the options described in the table below. These parameters can also be randomized like the parameters of the base model, see section on *Randomizing parameter values* for more information.

| parameter | default | option to set value                       |
| --------- | ------- | ----------------------------------------- |
| *g*       | 1.0     | `--synapse-g value`                       |
| *V0*      | -20.0   | `--synapse-V0 value`                      |
| *Vsyn*    | 20.0    | `--synapse-Vsyn value`                    |
| *tau_r*   | 0.5     | `--synapse-tau-r value`                   |
| *tau_d*   | 10      | `--synapse-tau-d value`                   |

By default *sbmodelr* will create one instance of the *g* parameter per synapse, and each one can later be changed independently within a modeling program like COPASI. However, sometimes it may be useful to have all the *g* parameters linked to a single master parameter. This can be achieved with the option `--synapse-link-g` which will then create an extra parameter called `g_c_{ode}_synapse` (where *{ode}* is the name of the variable that represents membrane potential) and all other *g* parameters will be linked to that one.

This option is presently not supported if the output is a *Vivarium* multiscale simulation.

## Randomizing parameter values

When replicating the base unit into several others, *sbmodelr* is able to allow specific parameters to take values different from the original. This is useful to generate populations of models that have parameter values distributed randomly. Currently *sbmodelr* allows parameter values to be sampled from uniform or normal distributions. In both cases the distributions are defined based on the value that the parameter takes in the base model.

The main option to specify a parameter to be randomized is `--pn parameter level distribution` where `parameter` is the name of the parameter, `distribution` is either `uni` or `norm` (for *uniform* and *normal* distributions). The `level` is a positive number that is the proportion of the parameter value that is added and subtracted from the parameter value, or the standard deviation, to characterize the distribution; see table below for more detail. Note that many parameters can be randomized with this option, simply add several `--pn parameter level distribution` entries, one for each parameter. Parameters can be any global quantities that are *fixed*, the initial concentrations of species , and initial sizes of compartments. For initial concentrations or initial sizes simply use the species or compartment name (in this context *sbmodelr* interprets a species or compartment name as their initial values).

The option above (`--pn`) does not cover the additional parameters that are created to establish the unit connections. To randomize *all* those parameters use the option `--cn level dist`, with `level` and `dist` having the same meaning as above.

The table below describes in more detail how these options work. In all cases *value* is the value that the parameter has in the base model.

| option                     | meaning |
| -------------------------- | ------- |
| `--pn parameter level uni` | parameter is sampled from a <ins>uniform</ins> distribution between (*value* - *level* * *value*) and (*value* + *level* * *value*) |
| `--cn level uni` | all connection parameters are sampled from a <ins>uniform</ins> distribution between (*value* - *level* * *value*) and (*value* + *level* * *value*) |
| `--pn parameter level norm` | parameter is sampled from a <ins>normal</ins> distribution with mean of *value* and standard deviation of (*level* * *value*) |
| `--pn parameter level norm` | all connection parameters are sampled from a <ins>normal</ins> distribution with mean of *value* and standard deviation of (*level* * *value*) |

As an example, the command `--pn Km 0.5 uni` applied to a model that has a parameter named *Km* with initial value 10 will result in the parameter *Km* in the replicas to have values that are sampled from a uniform distribution in the interval \[5,15\] (5=10-10\*0.5 and 15=10+10\*0.5). The command `--pn Km 0.5 norm` applied to the same model results in the parameter *Km* of the replicas to have values sampled from a normal distribution with mean of 10 and standard deviation of 5 (=10*0.5).

When randomizing connection parameters, *sbmodelr* will create one instance of that parameter for each replica. If no randomization is requested then there is only one parameter for all replicas (since the replicates would have the same exact value). In some cases it may be useful to force *sbmodelr* to create replicas of the connection parameters even though one may not wish to randomize them. This can be easily achieved with the option `--cn 0 uni` which will replicate the connection parameters but will assign them the same exact value (since the amplitude of the interval is zero). This is illustrated in the example [Neuron networks I case 3](https://github.com/copasi/sbmodelr/blob/main/examples/Neuron_networks_I/README.md#case-3).

This option is presently not supported if the output is a *Vivarium* multiscale simulation.

## Adding a medium compartment

In some situations one may want to create several replicas of a model and add an extra unit called *medium* that only contains the species that are transported, such that the species can be transported between the *medium* and all other units. Importantly, the *medium* compartment will not contain any reactions contained in the base unit. This can be achieved using the option `--add-medium`. The *medium* volume is set to a default of 1.0, but other values can be specified with option `--medium-volume value`.

A *medium* compartment mimics the situation in an experiment where cells are suspended in a medium and where some substances can move between the medium and the interior of the cell. A replica of a microbial culture where each microbe only interacts with the medium (*i.e.* not directly cell-cell) can be achieved by a command that does not connect the units, like this one:

        sbmodelr -t nutrient --add-medium mycell.cps 100

this would create a new model (called *mycell_100.cps*, as no name was specified) that contains 100 replicas of the base model (*mycell.cps*) where the species called *nutrient* can move between each cell and the medium (in both directions). The new model will include a compartment called *medium* that only contains species *nutrient*. The user can then load that model into COPASI and change the initial concentration of *nutrient* in the medium and set it to fixed, creating a gradient of nutrient concentration such that it gets transported into the cells.

This option is presently not supported if the output is a *Vivarium* multiscale simulation.

## Events

When models contain events these are also taken into consideration in the replication process. Events are dealt with in two different ways:

 - **Events that have a trigger only depending on a function of time** (without including any other model element in the trigger) are not replicated. However their targets are replicated appropriately. For example if there is an event when *Time* passes 10 (trigger) that changes the concentration of *S* (target), then the new model will continue to have only one event also triggered when *Time* passes 10, but now it has many targets, all the *S_1* ... *S_n* (*n* being the number of replicates) that get set to the same functions. If any other elements appear in the function, they will correspond to the same unit as the *S_i*.
 - **Events that depend on model elements** are entirely replicated, thus one event becomes *n* different events (for *n* replicate units). However in this case the targets will be changed to the respective unit (not any other unit). For example if there is an event that happens when species *Signal* becomes larger than 10, which then changes a quantity *tick* to be incremented, in the new model there will be *n* events, each one happening when *Signal_i* becomes larger than 10, which then changes *tick_i* to be incremented by one.

## Units, metadata and comments

All *units* used in the base model will also be used in the output model; model elements that are replicated inherit the units of the base model.

All *metadata* included in the base file are also copied to the new model. Each model element replicated will also get the same metadata as the original element in the base model. Additionally the model creation time will be a copy of the model creation time of the base model; the current date/time is added to a new modified time.

*Comments* (*i.e.* free text) attached to base model elements are copied to all the replicates in the new model. The model comments will include a statement that the model was created by *sbmodelr* and a copy of the full command line used.

## Task processing

When *sbmodelr* reads a base file that is in COPASI format (*.cps) it normally copies the settings of the tasks to the new file. To prevent this, and just create a file with default values for the tasks, add the option `--ignore-tasks`.

In the *Parameter Scan*, *Sensitivities*, *Cross Section* and *Optimization* tasks all the elements used will be translated to those of the first unit of the new model. In *Optimization* this includes both the objective function, the parameters and constraints. This is arbitrary and if the user requires that it reflects a different unit, it will have to be changed manually within COPASI.

## License

The software *sbmodelr* is Copyright © 2024-2026 Pedro Mendes, [Center for Cell Analysis and Modeling](https://health.uconn.edu/cell-analysis-modeling/), UConn Health. It is provided under the Artistic License 2.0, which is an OSI approved license. This license allows non-commercial and commercial use free of charge.

## Funding

This package was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number GM137787 as part of the [National Resource for Mechanistic Modeling of Cellular Systems](https://compcellbio.org/). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

