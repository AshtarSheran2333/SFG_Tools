# corr_to_spectrum

Once the correlation function is prepared, this program can be used to recalculate the spectrum with different [$TEMPERATURE](../user_guide/files.md/#table-of-all-boxdata-file-parameters) (affects only magnitude of the spectrum according to the [equation](./SSP_CORR.md/#spectrum-calculation)) or [$FILTER](../user_guide/files.md/#table-of-all-boxdata-file-parameters) (affects shape of the spectrum).


## Input and output file table

| input files | output files |
| --- | --- |
| [BOXDATA](../user_guide/files.md/#table-of-all-boxdata-file-parameters) | [*-spectrum](../user_guide/files.md/#-spectrumdat) |
| [*terms.dat](../user_guide/files.md/#-selftermsdat) | [*-shiftedspectrum](../user_guide/files.md/#-shiftedspectrumdat) |
