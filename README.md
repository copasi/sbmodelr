# *sbmodelr* - a tool to replicate one SBML/COPASI model into a set of replicas

## Summary
This is a python-based command line utility (*sbmodelr*) that reads an SBML or COPASI model file and then replicates the existing model into several units. These units may be organized as a network of arbitrarily connected units, a 2D rectangular grid, or a 3D cuboid array. Each unit contains a complete copy of the original model with all its species, reactions, compartments, events, and global quantities.

Connections between units in the new model can be:
 - species being transported between units
 - species acting as inhibitors/activators of the synthesis of other species (to make gene regulatory networks)
 - diffusive coupling of explicit ODEs ("rate rules" in SBML)
 - coupling of explicit ODEs through chemical synapse terms, appropriate for ODEs representing membrane potentials

An additional unit can be added -- called 'medium' -- which only contains the transported species, but is connected to all other units.

It is also possible to add randomness to parameter values, such that each unit becomes slightly different from each other.

Practical uses of *sbmodelr* include:
 - using a cell model to create a model of a tissue or organoid, 
 - use a gene transcription model to create a gene regulatory network
 - use a neuron model (e.g. the Hodgkin-Huxley) to create a network of neurons

The output of this program is a new model file with the more complex model. It is expected that the user may still have to tune parts of the resulting model in a regular modeling tool, such as [COPASI](https://copasi.org), [VCell](https://vcell.org), etc, where the model will be used for simulations. (*sbmodelr* only creates models, it does not carry out simulations.)

## Installation

The package works with python 3.8+ and requires the package *copasi-basico* (freely available on pypi).

You can install *sbmodelr* and its dependency by running:

        pip install git+https://github.com/copasi/sbmodelr.git

or optionally for development with:

        git clone https://github.com/copasi/sbmodelr
        pip install -e ./sbmodelr

<!--
You can install the package directly from pypi:

    pip install sbmodelr

-->

## Usage

See [User Manual](UserManual.md) for complete description of how to use *sbmodelr*. Detailed examples are provided in the `examples` folder.

## Credits

This program is inspired by [MEG](http://www.gepasi.org/meg.html) [1], a utility included in the old [Gepasi](http://www.gepasi.org) simulator. The COPASI GUI and the [BasiCO](https://github.com/copasi/basico) python API [2] both contain some functionality similar to that provided here, however they are limited to replicating compartments (with all their species and reactions) and connecting them by transport of species, but *do not* operate on global quantities, and can't add chemical synapse connections or regulatory interactions.

Thanks to Frank Bergmann for making [BasiCO](https://github.com/copasi/basico) and the whole [COPASI](https://copasi.org) team for that simulator, which is ultimately the backend that is working behind *sbmodelr*.

**References**
 1. [Mendes P, Kell DB (2001) MEG (Model Extender for Gepasi): a program for the modelling of complex, heterogeneous, cellular systems. Bioinformatics 17:288–289](https://doi.org/10.1093/bioinformatics/17.3.288)
 2. [Bergmann FT (2023) BASICO: A simplified Python interface to COPASI. Journal of Open Source Software 8:5553](https://doi.org/10.21105/joss.05553)

## Funding

This package was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number GM137787. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## License

This program is Copyright Pedro Mendes 2024 and is provided under the Artistic License 2.0, which is an OSI approved license. This license allows non-commercial and commercial use free of charge.
