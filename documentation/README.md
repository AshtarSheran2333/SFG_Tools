# SFG\_TOOLS SFG spectrum analysis toolkit

Welcome to the **SFG_TOOLS** documentation.

This set of programs processes trajectories of water molecules and surface hydroxyls to obtain vibrational sum-frequency generation O-H spectra ( $ Im \chi^{(2)} $ ).

# Table of Contents

1. [Getting Started](./user_guide/getting_started.md)
2. **Usage**
    - [Program switches](./user_guide/program_switches.md)
    - [Files description](./user_guide/files.md)
3. **Reference**
    - [Interface Program](./programs_description/Interface.md)
        + [Instantanous surfaces](./programs_description/Interface.md/#calculation-of-instantaneous-surfaces)
        + [Density function](./programs_description/Interface.md/#calculation-of-density-function)
    - [Binder Program](./programs_description/Binder.md)
        + [Dipole analysis](./programs_description/Binder.md/#dipole-analysis)
        + [Hydrogen bond distributions](./programs_description/Binder.md/#hydrogen-bond-distributions)
    - [SSP_CORR Program](./programs_description/SSP_CORR.md)
        + [Construction of D matrices and evaluation of vz](./programs_description/SSP_CORR.md/#construction-of-d-matrices-and-evaluation-of-vz)
        + [Calculation of the correlation function of mz and axx](./programs_description/SSP_CORR.md/#calculation-of-the-correlation-function-of-mz-and-axx)
        + [Filter function](./programs_description/SSP_CORR.md/#filter-function)
        + [Spectrum calculation](./programs_description/SSP_CORR.md/#spectrum-calculation)
    - [corr_to_spectrum Program](./programs_description/corr_to_spectrum.md)
4. **Appendix**
    - [Program call order](./user_guide/program_call_order.md)
