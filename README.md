> [!WARNING]
> This is a work in progress branch.

# WORK IN PROGRESS
~~In this branch I am working on a generalized, modular file reader and its incorporation to the code.~~
This branch has transformed into dump of ideas that I am trying to work on...
Following things needs to be done:
- [x] the frame reader module
    + [x] trr 
    + [x] gro 
    + [x] xyz 
- [x] Some layer above the reader that will take care of the molecular structure
    + [x] in a separate structure file
    + [x] struct file format clearly defined
    + [x] struct format parser should be failproof
    + [x] The MOL should have a better name like GROUP
    + [x] User should be able to define multiple GROUP, each GROUP will have an optional name e.g. GROUP CH3
- [ ] Atleast the important files storing mid calculations should have their own modules with in and out procedures and shared internal structure...
    + [ ] binder
    + [ ] interface
- [ ] There could be a module for hydrogen bond analysis (3D histograms...) it needs fixing anyway
- [ ] A general density module (can be reused for both density and dipole)
- [ ] better BOXDATA module
    + [x] better PBC checking (by default origin in 0 0 0, possibility to move the origin...)
- [ ] Go through the modules and properly setup the interface (public, private, protected, implicit none in all functions and subroutines... hide the internal subroutines etc.)
    + [x] the functions should have implicit none everywhere
    + [ ] the interface needs to be properly set (private, public routines)
    + [ ] SFG_STRUCTURE - what a naming madness - needs to be fixed
- [ ] Incorporate into all programs
    + [ ] Interface_P (rename to Interface)
        - Since we want to calculate W-C interface, something must be done about it, so it is able to work with various mixtures e.g. water + organic ions
        - Since the trajectory can contain also a solid surface, there must be a way how to select the liquid atoms (the heavy ones), for the interface calculation we can keep the same coarse graining length as for water, but do some mass weighting on the contributions
        - another file?? NOBODY IS EVER GOING TO FILL IN THOSE CRAZY INDEX FILES... lol
    + [ ] Binder_P (rename to Binder)
        - there will be two binder files, one for water, second for OTHERS
        - the SFG_STRUCT must then be somehow tied to the binder... so we can catch wheather it has been changed...
        - the molecules or chromophores or whatever "units" of spectrum must be somehow associated to some distance from the interface (probably the COM)
    + [ ] SSP_CORR_P (rename to SFG_CORR)
        - Enable calculation of multiple correlation functions simultaineously e.g. water L0, water L2-3, other CH3, hydroxyls
- [ ] The programs should be able to somehow backup the results in case of IO errors
        - The reader modules should assure this - each function returns some kind of return code -> catch those, backup, abort
- [ ] Start making some tests (probably in form of a special program...)
- [ ] Work on the CMakeLists.txt
    + [ ] use the f18 standard flag
    + [ ] during the cleanup, the code should be more standardized, so possibly working even with gfortran??
- [ ] Update documentation

# SFG\_TOOLS - sum-frequency generation spectrum analysis toolkit

Welcome to the **SFG_TOOLS** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17517850.svg)](https://doi.org/10.5281/zenodo.17517850) documentation.

The code was released in version 1.0.2 with the JCTC article: [10.1021/acs.jctc.5c02160](https://doi.org/10.1021/acs.jctc.5c02160)

This set of programs processes trajectories of water molecules and surface hydroxyls to obtain vibrational sum-frequency generation O-H spectra $\chi^{(2)}(\omega)$.

# Table of Contents

1. [Getting Started](./documentation/user_guide/getting_started.md)
2. **Usage**
    - [Program switches](./documentation/user_guide/program_switches.md)
    - [Files description](./documentation/user_guide/files.md)
3. **Reference**
    - [Interface Program](./documentation/programs_description/Interface.md)
        + [Instantanous surfaces](./documentation/programs_description/Interface.md/#calculation-of-instantaneous-surfaces)
        + [Density function](./documentation/programs_description/Interface.md/#calculation-of-density-function)
    - [Binder Program](./documentation/programs_description/Binder.md)
        + [Dipole analysis](./documentation/programs_description/Binder.md/#dipole-analysis)
        + [Hydrogen bond distributions](./documentation/programs_description/Binder.md/#hydrogen-bond-distributions)
    - [SSP_CORR Program](./documentation/programs_description/SSP_CORR.md)
        + [Construction of D matrices and evaluation of vz](./documentation/programs_description/SSP_CORR.md/#construction-of-d-matrices-and-evaluation-of-vz)
        + [Calculation of the correlation function of mz and axx](./documentation/programs_description/SSP_CORR.md/#calculation-of-the-correlation-function-of-mz-and-axx)
        + [Filter function](./documentation/programs_description/SSP_CORR.md/#filter-function)
        + [Spectrum calculation](./documentation/programs_description/SSP_CORR.md/#spectrum-calculation)
    - [corr_to_spectrum Program](./documentation/programs_description/corr_to_spectrum.md)
4. **Appendix**
    - [Program call order](./documentation/user_guide/program_call_order.md)
