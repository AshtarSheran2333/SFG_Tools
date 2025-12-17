# Program switches
This section describes required and optional options for the programs. 

## Interface
- -I(input) &lt;**arg**&gt; (optional)&lt;**arg**&gt;
	- Takes up to 2 arguments. Pay attention to the order of input files.
	- In case of [.xyz files](../user_guide/files.md/#xyz-files): position file + matching velocity file.
	- In case of [.gro file](../user_guide/files.md/#gro-file): just multiple frame [.gro file](../user_guide/files.md/#gro-file) with velocities.
	- In case of [.trr file](../user_guide/files.md/#trr-file): [.trr file](../user_guide/files.md/#trr-file) + atleast one frame of matching [.gro file](../user_guide/files.md/#gro-file).

- (optional) "-G(good)"
	- skip checking of the input file contents (if you are sure that the input files have proper format, number of frames...)

- (optional) "-V(vmdout)"
	- Program will also create [interface.xyz](../user_guide/files.md/#interfacexyz-optional) file in case of need to visualize the interface using visualization software. 

- (optional) "-H(help)"
	- Help dialog appears when selected.

## Binder
- -I(input) &lt;**arg**&gt; (optional)&lt;**arg**&gt;
	- Takes up to 2 arguments. Pay attention to the order of input files.
	- In case of [.xyz files](../user_guide/files.md/#xyz-files): position file + matching velocity file.
	- In case of [.gro file](../user_guide/files.md/#gro-file): just multiple frame [.gro file](../user_guide/files.md/#gro-file) with velocities.
	- In case of [.trr file](../user_guide/files.md/#trr-file): [.trr file](../user_guide/files.md/#trr-file) + atleast one frame of matching [.gro file](../user_guide/files.md/#gro-file).

- (optional) -D(dipole)
	- Program will also calculate [average water dipole orientation](../programs_description/Binder.md/#dipole-analysis) based on distance from the instantaneous surface the results are then stored in [dipolevsr.dat](../user_guide/files.md/#dipolevsrdat-optional).

- (optional) -B(bonds)
	- Program will also [analyse hydrogen bond networks](../programs_description/Binder.md/#hydrogen-bond-distributions) in all the [water layers](../programs_description/Binder.md). Data will be stored in the [Histogram*_.3dhist files](../user_guide/files.md/#histogram_3dhist-files-optional). File [hbonds.dat](../user_guide/files.md/#hbondsdat-optional) will be generated.

- (optional) -G(good)
	- skip checking of the input file contents (if you are sure that the input files have proper format, number of frames...)

- (optional) -H(help)
	- Help dialog appears when selected.

## SSP_CORR
- -I(input) &lt;**arg**&gt; (optional)&lt;**arg**&gt;
	- Takes up to 2 arguments. Pay attention to the order of input files.
	- In case of [.xyz files](../user_guide/files.md/#xyz-files): position file + matching velocity file.
	- In case of [.gro file](../user_guide/files.md/#gro-file): just multiple frame [.gro file](../user_guide/files.md/#gro-file) with velocities.
	- In case of [.trr file](../user_guide/files.md/#trr-file): [.trr file](../user_guide/files.md/#trr-file) + atleast one frame of matching [.gro file](../user_guide/files.md/#gro-file).
		
- "-L(layers) &lt;**arg1**&gt; &lt;**arg2**&gt;
	- Selection of layers to be analysed (see [Binder](../programs_description/Binder.md)).
		- **arg1** - from Layer number: 0 - 15
		- **arg2** - to Layer number 0-15

- (optional) -O(output) &lt;**arg**&gt;
	- Takes output file tag. By default the output file tag is empty string.

- (optional) -S(selfterms)
	- program will calculate the [self correlation function](../user_guide/files.md/#-selftermsdat)

- (optional) -C(crossterms) &lt;**arg**&gt;
	- Program will also calculate [the crossmolecular terms](../user_guide/files.md/#-crosstermsdat-optional) with cutoff &lt; **arg** &gt; (in Angstroms)

- (optional) -N(neighborlist)
	- neigborlist will be used when calculating [the crossmolecular terms](../user_guide/files.md/#-crosstermsdat-optional), this reduces calculation time for larger systems

- (optional) -A(AMSUM)
	- program will sum [dipole moment derivatives](../programs_description/SSP_CORR.md) and [polarizability derivatives](../programs_description/SSP_CORR.md) of the selected layer and calculate correlation function of the properties - This is mostly an experimental feature that is not that useful. The main idea was to substract [self correlation terms](../user_guide/files.md/#-selftermsdat), and obtain the [cross correlation terms](../user_guide/files.md/#-crosstermsdat-optional) for a fraction of computation time. One cannot tune cross correlation cutoff, so all the cross correlation terms from the selected layer are included, but the cross correlation terms between the selected layer and everything else are missing. This does not seem to produce relevant results.

- (optional) -D(Dipole_moment_and_polarizability_derivatives) &lt;**arg**&gt;
	- Takes [derivatives file](../user_guide/files.md/#input--deriv--file-optional) as an argument.
	- Allows to change $\frac{d\mu_x}{dr_z}$, $\frac{d\mu_y}{dr_z}$, $\frac{d\mu_z}{dr_z}$, $\frac{d\alpha_{xx}}{dr_z}$, $\frac{d\alpha_{yy}}{dr_z}$, $\frac{d\alpha_{zz}}{dr_z}$, $\frac{d\alpha_{xy}}{dr_z} = \frac{d\alpha_{yx}}{dr_z}$, $\frac{d\alpha_{xz}}{dr_z} = \frac{d\alpha_{zx}}{dr_z}$, and $\frac{d\alpha_{yz}}{dr_z} = \frac{d\alpha_{zy}}{dr_z}$ parameters.

- (optional) -F(Filter) &lt;**arg**&gt;
	- Takes real number [ps] as parameter to the [filter function](../programs_description/SSP_CORR.md/#filter-function) used on the correlation function; $\text{filter}(\tau) = \exp{\left(\frac{tau}{\lt arg \gt}\right)^2}$.

- (optional) -G(Good)
	- The contents of input files will not be checked. This can save some time **if the user is sure** that all the input files contains all the frames needed for the calculation.

- (optional) -H(Help)
	- Help dialog appears when selected.

## corr_to_spectrum
- -i &lt;**arg**&gt;
	- takes name of correlation function file

- (optional) -filter &lt;**arg**&gt;
	- Takes real number [ps] as parameter to the [filter function](../programs_description/SSP_CORR.md/#filter-function) used on the correlation function; $\text{filter}(\tau) = \exp{\left(\frac{tau}{\lt arg \gt}\right)^2}$.

- (optional) -o &lt;**arg**&gt;
	- if not selected the output will be named spektrum.xyz
