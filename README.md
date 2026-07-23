# Infectious Disease Dynamics: a systems approach

This Github repository contains notebooks developed for use in the the University of Maryland course BSCI439C/BIOL708F for Winter 2025, led by Dr. Stephen Beckett and Dr. Gabi Steinbach. 
The interactive notebooks are designed to help users explore epidemiological models describing spread of infectious disease, and how different model structures and parameters can influence epidemiological population dynamics.
The notebooks have been developed in both Pluto.jl (using Julia) and marimo (using Python) frameworks, and both contain interactive elements.

### Pluto.jl notebooks  (file.jl)
Offline🖥️: Requires Julia: https://julialang.org/. Can be used after installation of the Pluto.jl package: https://github.com/fonsp/Pluto.jl

Online🌐:  Alternatively, could be loaded via binder at: https://binder.plutojl.org/  (note, I experienced binder timing out during the compilation process) 

### marimo notebooks (file.py)
Offline🖥️:  Requires Python: https://www.python.org/. Can be used after installation of the marimo package: https://docs.marimo.io/getting_started/installation/

Online🌐:   Alternatively, can be loaded via the interactive playground at https://marimo.io.

## Notebook summary
*Week 1:* file: **EpiNotebook.jl / EpiNotebook.py**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;These notebooks explore the SIR model, the SIRS model, and the stochastic SIR model.

*Week 2:* file: **ExtendedEpiModels.jl / ExtendedEpiModels.py**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;These notebooks explore the SEIR model, and the SEIaIsRD model (incorporating asymptomatic infections).

*Week 3:* file: **InterventionNotebook.jl / InterventionNotebook.py**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;These notebooks explore the SEIRD model, and using it to drive and compare outbreak intervention scenarios.

## Version notes
Software is continually revised and updated. For longevity, we note software and package versions at a stage where these notebooks have been run. This could be useful for setting up these files through containerization in the future.

.py files last checked with Python 3.14.2. Package versions are listed in python_packges.txt

.jl files last checked with Julia 1.11.5. Package versions are listed in julia_packages.txt
