> [!WARNING]
> This is a work in progress branch.

# WORK IN PROGRESS
In this branch I am working on a generalized, modular file reader and its incorporation to the code.
Following things needs to be done:
- [x] the frame reader module
    + [x] trr 
    + [x] gro 
    + [x] xyz 
- [ ] Some layer above the reader that will take care of the molecular structure
    + [ ] probably some indexlike parameters in BOXDATA
    + [ ] give user responsibility for the indexes, later can be appended to automatically find some structures (it will be pain, since the atomnames are not unified etc...)
- [ ] Incorporate into all programs
    + [ ] Interface
    + [ ] Binder
    + [ ] SSP_CORR_P (rename to SFG_CORR)
- [ ] The programs should be able to somehow backup the results in case of IO errors
- [ ] Work on the CMakeLists.txt

# SFG\_TOOLS - sum-frequency generation spectrum analysis toolkit

Welcome to the **SFG_TOOLS** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17517850.svg)](https://doi.org/10.5281/zenodo.17517850) documentation.

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
