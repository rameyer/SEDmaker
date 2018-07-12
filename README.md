SEDmaker
======
A quick python script to produce SEDs from BPASS models for a range of SFHs.

## Still under work - only simple Star Formation Histories and BPASS models for now. Output is still formatted for HyperZ.


### Installing SEDmaker 

Provided you have the right Python (3.6) with the distutils build, and updated numpy and matplotlib packages, this should be as easy as:

```
$ git clone https://github.com/rameyer/SEDmaker
$ cd SEDmaker
$ python3 setup.py install
$ cd ..
$ python3
>>> import SEDmaker
```

If this fails, you may need more detailed installation instruction. Please contact me!

### Downloading the BPASS models

You need first to download the BPASS models (Stanway, Eldridge et al., 2017; http://bpass.auckland.ac.nz/) before SEDmaker works. Please unpack each IMF folder in the BPASS folder provided in the git repo. 


### Getting started

SEDmaker is a small package to help produce Spectral Energy Distributions based on the BPASS models. In order to do so, SEDmaker has a range of Star Formation Histories that can be chosen from and adjusted for different age, e-folding time, slope of the Star Formation Rate (SFR) decrease etc... This is then supplemented by the different BPASS models for various Initial Mass Functions (IMF), metallicities and inclusion of binary star populations.

To see the range of possible SFHs:

```
$ python3
>>> import SEDmaker.SEDmaker as sed
>>> sed.test_SHFs()
```

To produce a grid of HyperZ SEDs, call the following function:
```
def produce_grid_SED_HyperZ(type = 'constant', age_array = [10], metallicity_array = ['001'],
                            params = [1], min_wavelength = 0, max_wavelength = 6000,
                            sampling = 1, imf= 'imf135_300', sed_save_path = PACKAGE_PATH+'/SED/',
                            BPASSv = '2.1'):
	''' 
	An example of a grid production of SEDs in HyperZ specific format.
	Still under work to be fully determined from inputs.
	Inputs:
		type: either constant, linear, exponential
		age_array: the array of ages of the galaxies you want to produce the SEDs of
		metallicity_array: metallicities of the SEDs, in the BPASS format
		params: array of slope(e-folding time) for the linear(exponential) SFH case.
		min_wavelength: minimum wavelength of the output SEDs (integer)
		max_wavelength: maximum wavelength of the output SEDs (integer)
		sampling: parameter to undersample the BPASS SEDs by 1 (no sampling), 2,3,4,.. Integer.
		imf: The imf to use. Must be one of those used by the BPASS team
		sed_save_path: The directory to save the SEDs, by default the SED
		BPASSv: BPASS version, in case you have many.
	Outputs:
		-- (save the SEDs as text files in given directory + hyperz .param files)
	'''
```
