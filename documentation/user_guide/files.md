# Input trajectory files
This section covers all the possible trajectory file formats supported by the tools.
> [!NOTE]
> All te input trajectory files should have its origin in the center of the simulation box.
		
## *.xyz files
Two separate .xyz files (one to store positions of atoms, second to store velocities) with standard .xyz structure. First line contains number of atoms, second line is a comment line followed by number of atoms lines corresponding to each atom position or velocities. Positions should be expressed in Angstroms, velocities should be expressed in Hartree atomic units $a_{0}E_{\text{h}}/\hbar$. 
				
## *.gro file
See [Gromacs manual](https://manual.gromacs.org/current/reference-manual/file-formats.html#gro).
		
## *.trr file
> [!IMPORTANT]
> Make sure that nstxout matches nstvout in your .mdp file. Note that if nstxout and nstvout is not 1, the [BOXDATA](./files.md#input-boxdata-file) [$DT](./files.md#table-of-all-boxdata-parameters) does not match .mdp file, and [$DT](./files.md#table-of-all-boxdata-parameters) should be corrected, the maximal [$DT](./files.md#table-of-all-boxdata-parameters) should be in accordance with Nyquist-Shannon sampling theorem counting with the [$FREQ](./files.md#table-of-all-boxdata-parameters).
>
> For example:			
>
> dt in .mdp is set to 1 fs, nstxout and nstvout is set to 1, then dt in .mdp matches [BOXDATA](./files.md#input-boxdata-file) [$DT](./files.md#table-of-all-boxdata-parameters).
>
> Yet another example:
>
> dt in .mdp is set to 1 fs, nstxout and nstvout is set to 2, that means one frame is skipped when writing to the .trr file, then .mdp dt does not match [BOXDATA](./files.md#input-boxdata-file) [$DT](./files.md#table-of-all-boxdata-parameters). [BOXDATA](./files.md#input-boxdata-file) [$DT](./files.md#table-of-all-boxdata-parameters) should be now set to 2 fs.
			
See Gromacs manual for [.trr files](https://manual.gromacs.org/current/reference-manual/file-formats.html#trr) and [.mdp file options](https://manual.gromacs.org/current/reference-manual/file-formats.html#mdp).
							
# Input BOXDATA file
File called BOXDATA contains all the user defined parameters of SFG spectra calculation. Note that if a parameter has more than one operand, each operand needs to be divided by new line character.
			
## Table of all BOXDATA file parameters

| parameter identifier           | n. operands | type    | unit                           | default value     | brief description 
| --- | --- | --- | --- | --- | --- |
| \$BOX\_DIMENSIONS              | 3           | REAL    | A                              | (/-1,-1,-1/)	  | dimensions of the simulation box |
| \$INTERFACE\_VOLUME\_ELEMENT   | 3           | REAL    | A                              | (/0.5,0.5,0.25/)  | size of volume element while calculating instantaneous surface |
| \$LAYERS\_LIMITS               | (up to) 7   | REAL    | A                              | (/0,0,0,.../)	  | height of water layers respective to the instantaneous surface |
| \$FREQ                         | 1           | REAL    | cm-1                           | 4000              | maximum frequency analysed |
| \$DFREQ                        | 1           | REAL    | cm-1                           | 1                 | frequency step |
| \$DT                           | 1           | REAL    | fs                             | -1                | timestep of the simulation |
| \$CORRLEN                      | 1           | REAL    | ps                             | 0                 | maximum lag of the correlation function |
| \$FILTER                       | 1           | REAL    | ps                             | 0                 | filter parameter |
| \$TEMPERATURE                  | 1           | REAL    | K                              | 300               | temperature of the simulation box |
| \$HBHIST\_ANGLE\_START         | 1           | REAL    | cos(angle)                     | -1                | hbhist y axis start value |
| \$HBHIST\_ANGLE\_END           | 1           | REAL    | cos(angle)                     | 1                 | hbhist y axis end value |
| \$HBHIST\_DIST\_START          | 1           | REAL    | A                              | 2.3               | hbhist x axis start value |
| \$HBHIST\_DIST\_END            | 1           | REAL    | A                              | 3.2               | hbhist x axis end value |
| \$DIPOLE\_R\_START             | 1           | REAL    | A                              | -5                | dipole\_vs\_r x axis start value |
| \$DIPOLE\_R\_END               | 1           | REAL    | A                              | 80                | dipole\_vs\_r x axis end value |
| \$DIPOLE\_R\_BINWIDTH          | 1           | REAL    | A                              | 0.025             | dipole\_vs\_r binwidth |
| \$DENSITY\_MIN\_R              | 1           | REAL    | A                              | -4                | density\_vs\_r x axis start |
| \$DENSITY\_MAX\_R              | 1           | REAL    | A                              | 25                | density\_vs\_r x axis end |
| \$DENSITY\_TOL                 | 1           | REAL    | A                              | 0.7               | density\_vs\_r bin overlap |
| \$INTERFACE\_DENSITY           | 1           | REAL    | -                              | 0.5               | the limit of coarse grained density to be marked as instantaneous surface |
| \$NATOM                        | 1           | INTEGER | -                              | -1                | number of atoms in simulation |
| \$NO                           | 1           | INTEGER | -                              | -1                | number of oxygens in simulation |
| \$NSTEP                        | 1           | INTEGER | -                              | -1                | number of simulation steps to be analysed |
| \$INTERFACE\_SKIP              | 1           | INTEGER | -                              | 0                 | number of skipped frames in order to calculate instantaneous surface |
| \$CROSS\_SKIP                  | 1           | INTEGER | -                              | 0                 | number of skipped time beginnings when calculating the cross-correlation function |
| \$SELF\_SKIP                   | 1           | INTEGER | -                              | 0                 | number of skipped time beginnings when calculating the self-correlation function |
| \$INTERFACE\_PUSHBACK          | 1           | INTEGER | Interface_volume_element(3)    | 100000000         | pushback when searching for new instantaneous surface |
| \$HBHIST\_ANGLE\_DIV           | 1           | INTEGER | -                              | 100               | hbhist number of y bins |
| \$HBHIST\_DIST\_DIV            | 1           | INTEGER | -                              | 45                | hbhist number of x bins |
| \$HBHIST_SKIP					 | 1		   | INTEGER | -							  | 0				  | number of skipped frames when evaluating hydrogen bond analyses |
| \$DENSITY\_BIN                 | 1           | INTEGER | -                              | 5                 | density\_vs\_r number of bins in 1 A |
| \$HYDROXYL_METAL				 | 1		   | CHAR(10) | -							  |	'X'				  | hydroxyl metal name |
| \$CANCEL_INTERFACE_CALCULATION | 0		   | NONE    | -                              |	-				  | [Interface program](../programs_description/Interface.md) will skip [calculation of instantaneous surface](../programs_description/Interface.md#calculation-of-instantaneous-surfaces), the instantaneous surfaces will be read from [interface.bin](./files.md#interfacebin) file |
| \$POLARIZATION				 | 1		   | CHAR(10) | -							  |	'SSP' (defaut) 'PPP' (allowed) | change between the chi xxz and chi zzz |

## Sample BOXDATA file

You can copy-paste content of this subsection, replace brackets by your values and save as BOXDATA to quickly prepare the input file. Note that each bracket represents an user input. Variable type, brief description and unit is contained in the brackets. Even an empty BOXDATA file can be provided, program should not start until errors in the BOXDATA file disappears, user will be informed about all the errors in BOXDATA file.
			
```
$BOX_DIMENSIONS
<REAL-X [Angstrom]>
<REAL-Y [Angstrom]>
<REAL-Z [Angstrom]>
$LAYERS_LIMITS
<REAL-L0/L1 [Angstrom]>
<REAL-L1/L2 [Angstrom]>
<REAL-L2/L3 [Angstrom]>
$NSTEP
<INTEGER-simulation steps [-]>
$DT
<REAL-simulation timestep [fs]>
$NATOM
<INTEGER-number of atoms in frame [-]>
$NO
<INTEGER-number of oxygens in frame [-]>
$INTERFACE_SKIP
<INTEGER-number of skipped frames in interface calculation [-]>
$INTERFACE_PUSHBACK
<INTEGER-pushback of testing interface in time t and t+1 [0.25 Angstrom]>
$CORRLEN
<REAL-max timelag of the correlation function [ps]>
$FREQ
<REAL-maximal analysed frequency [cm-1]>
$DFREQ
<REAL-frequency step [cm-1]>
$FILTER
<REAL-filter parameter [ps]>
```	

# Input -D (dipole_moment_and_polarizability_parameters) file (optional)
Structure of the file should be 9 lines, each of the line should contain only one number representing 
$\frac{dM_x}{dR_z}$, 
$\frac{dM_y}{dR_z}$, 
$\frac{dM_z}{dR_z}$, 
$\frac{dA_{xx}}{dR_z}$, 
$\frac{dA_{yy}}{dR_z}$, 
$\frac{dA_{zz}}{dR_z}$, 
$\frac{dA_{xy}}{dR_z} = \frac{dA_{yx}}{dR_z}$, 
$\frac{dA_{xz}}{dR_z} = \frac{dA_{zx}}{dR_z}$, and 
$\frac{dA_{yz}}{dR_z} = \frac{dA_{zy}}{dR_z}$ respectively.

These are the constants that will be used to calculate Axx and Mz see [SSP_CORR](../programs_description/SSP_CORR.md#ssp_corr). If the -D switch is not used, these values are set by default according to table below [^1]
					
| $\frac{dM_x}{dR_z}$ | $\frac{dM_y}{dR_z}$ | $\frac{dM_z}{dR_z}$ | $\frac{dA_{xx}}{dR_z}$ | $\frac{dA_{yy}}{dR_z}$ | $\frac{dA_{zz}}{dR_z}$ | $\frac{dA_{xy}}{dR_z}$ | $\frac{dA_{xz}}{dR_z}$ | $\frac{dA_{yz}}{dR_z}$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|0.005 | -0.047 | 1.028 | 1.168 | -0.247 | -1.331 | 0.319 | 0.281 | 3.444 |
			
# Program created files
This section covers all the files created by running the SFG calculation programs.
			
## interface.bin		
This file stores all the instantaneous surfaces calculated by program [Interface](../programs_description/Interface.md) in binary form (big endian).

The file contains:

1) header at the beginning:

| description | type | size (Bytes) |
| --- | --- | --- |
| number of points | integer | 8 |
| pair of x and y coordinate of a point | real * 8, real * 8 | number of points * 8 |

2) (0 and every [$NSKIP](./files.md#table-of-all-boxdata-parameters)th frame) number of points z coordinates:

| description | type | size (Bytes) |
| --- | --- | --- |
| upper and lower z coordinate for the point defined in header | real * 8, real * 8 | number of points * 8 |

> [!NOTE]
> the number of points should correspond to the [grid_interface](./files.md#grid_interface) file X*Y Npoint * 2
> e.g. number of frames 20, nskip is set to 9: the file has header, and two frames of z point coordinates.

## grid_interface
Plain text file containing information about division of the space while calculating [instantaneous surfaces](../programs_description/Interface.md#calculation-of-instantaneous-surfaces). Each row represents X,Y,Z axis respectively. Columns represent Number of points in given direction, minimal and maximal value in Angstrom respectively.
		
## density.dat
Data table containing a header marked by hashtag (\#) as a first character. The columns represent distance from the instantaneous surface in Angstrom, averaged density profile with respect to both upper and bottom interface, density profile with respect to upper interface and density profile with respect to bottom interface respectively. This file is generated by the [Interface program](../programs_description/Interface.md).
			
## interface.xyz (optional)
This file stores all the instantaneous surfaces calculated by program [Interface](../programs_description/Interface.md) in .xyz format, so the instantaneous surfaces can be visualised by standard tools. Creation of this file is subject of using -vmdout switch see [Interface program swicthes](./program_switches.md#interface).
	
## binder.bin
Binary file describing layer affiliation of each water molecule in each timestep of the simulation trajectory.
		
## hbonds.dat (optional)
Created by [-hbhist option](./program_switches.md#binder) calling [Binder program](../programs_description/Binder.md). Shows statistics of hydrogen bonds for each [layer](../programs_description/Binder.md):
			
- Average number of water molecules in specific [layer](../programs_description/Binder.md).
- Total average number of hydrogen bonds in specific [layer](../programs_description/Binder.md).
- Total average number of intralayer hydrogen bonds (donor and acceptor shares the same [layer](../programs_description/Binder.md)).
- Total average number of interlayer hydrogen bonds (donor and acceptor does not share the same [layer](../programs_description/Binder.md)).
- Average number of hydrogen bond interactions with other [layer](../programs_description/Binder.md).  
		
## dipole.dat (optional)
Data table containing a header marked by hashtag (\#) as a first character. The columns represent distance from the [instantaneous surface](../programs_description/Interface.md#calculation-of-instantaneous-surfaces) in Angstrom, average water dipole orientation (cosine of angle between dipole and Z axis of the box) as a fuction of distance from the closest point of the [instantaneous surface](../programs_description/Interface.md#calculation-of-instantaneous-surfaces) (bottom and top, respectively). This file is generated by the [Binder program](../programs_description/Binder.md).
			
## Histogram_*.3dhist (optional)
Data table containing a header marked by hashtag (\#) as a first character. The columns represent O-O distance (x-axis), cosine of O-O vector to the [interface normal](../programs_description/Binder.md#hydrogen-bond-distributions), and occurence probabilities for specific orientation (see table below) respectively. This file is generated by the [Binder program](../programs_description/Binder.md).
			
| column | data |
| --- | --- |
| 1 | *X* |
| 2 | *Y* |
| 3 | inter *Z* |
| 4 | intra *Z* |
| 5 | total *Z* |

## *-selfterms.dat
Data table containing only [self correlation terms T_{self}](../programs_description/SSP_CORR.md#calculation-of-the-correlation-function-of-mz-and-axx) of [calculated correlation function](../programs_description/SSP_CORR.md#calculation-of-the-correlation-function-of-mz-and-axx) of Axx and Mz.
		
## *-crossterms.dat (optional)
Data table containing only [cross correlation terms T_{cross}](../programs_description/SSP_CORR.md#calculation-of-the-correlation-function-of-mz-and-axx) of [calculated correlation function]() of Axx and Mz.

Three datacolumns included (different switch function $\zeta_L$ )
1) original switch function mentioned in [the program dexription](../programs_description/SSP_CORR.md#calculation-of-the-correlation-function-of-mz-and-axx)
2) switch function centered around dipole moment
```math
\zeta_L \left( \vec{r}_{m}(t), \vec{r}_{m'}(0)\right) =
    \begin{cases}
        1; &\mathbf{r}_{m'}(0) \text{ belongs to analysed layer }L  \\
        0; &\text{otherwise} 
    \end{cases}
```
3) switch function centered around polarizability
```math
\zeta_L \left( \vec{r}_{m}(t), \vec{r}_{m'}(0)\right) =
    \begin{cases}
        1; &\mathbf{r}_{m}(t) \text{ belongs to analysed layer }L  \\
        0; &\text{otherwise} 
    \end{cases}
```
## *-hydroxylterms.dat (optional)
Data table containing only hydroxyl [self correlation terms](../programs_description/SSP_CORR.md#calculation-of-the-correlation-function-of-mz-and-axx) of [calculated correlation function](../programs_description/SSP_CORR.md#calculation-of-the-correlation-function-of-mz-and-axx) of Axx and Mz.
				
## *-spectrum.dat
Data table containing a header marked by hashtag (\#) as a first character followed by the final SFG spectrum calculated by [SSP\_CORR program](../programs_description/SSP_CORR.md). Columns represent wavenumber, real part of the spectrum, imaginary part of the spectrum, amplitude of the spectrum and phase shift of the spectrum respectively.
			
## *-shiftedspectrum.dat
Data table containing a header marked by hashtag (\#) as a first character followed by the final SFG spectrum calculated by [SSP\_CORR program](../programs_description/SSP_CORR.md). Columns represent wavenumber, real part of the spectrum, imaginary part of the spectrum, amplitude of the spectrum and phase shift of the spectrum respectively. Data contained in this file are the same as in [*-spectrum.dat](./files.md#-spectrumdat), but the whole spectrum is shifted so the spectrum at frequency [\$FREQ](./files.md/#table-of-all-boxdata-file-parameters) is equal to zero.

[^1]: Khatib, R.; Backus, E. H. G.; Bonn, M.; Perez-Haro, M.-J.; Gaigeot, M.-P.; Sulpizi, M. Water Orientation and Hydrogen-Bond Structure at the Fluorite/Water Interface. Sci. Rep. *2016*, 6, 24287.
