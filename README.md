# SEDmaker
A quick python script to produce SEDs from BPASS models for a range of SFHs.

## Still under work - only simple Star Formation Histories and BPASS models for now. Output is still formatted for HyperZ.

## To work, please clone the git repo into the same directory as your BPASS downloaded files (or change paths in the code accordingly @ line 15) 

Check:
- test_SFHs() to get a flavor of the range of SFH available
- produce_grid_SED_HyperZ() gives a rough example of grid generation of SEDs. Modify the ranges of age, e-folding, slope (linear models) and metallicities as you wish. The output is only for HyperZ and gives you all necessary files (comprising the .param files).