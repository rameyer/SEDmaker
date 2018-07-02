SEDmaker
======
A quick python script to produce SEDs from BPASS models for a range of SFHs.

## Still under work - only simple Star Formation Histories and BPASS models for now. Output is still formatted for HyperZ.


### Installing SEDmaker 

Provided you have the right Python (3.6) with distutils build, and the numpy and matplotlib, this should be as easy as:

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


#### Getting started

SEDmaker is a small package to help produce Spectral Energy Distributions based on the BPASS models (see  Eldridge, Stanway et al.,2017 and http://bpass.auckland.ac.nz/). In order to do so, SEDmaker has a range of Star Formation Histories that can be chosen from and adjusted for different age, e-folding time, slope of the Star Formation Rate (SFR) decrease etc... This is then supplemented by the different BPASS models for various Initial Mass Functions (IMF), metallicities and inclusion of binary star populations.

To see the range of possible SFHs:

```
$ python3
>>> import SEDmaker as sed
>>> sed.SEDmaker.test_SHFs()
```

To produce a grid of HyperZ SEDs: 
```
$ python3
>>> import SEDmaker as sed
>>> sed.SEDmaker.produce_grid_SED_HyperZ()
```

ToDo: These will be more parametrized to allow the user to pick from the range of SFH listed [here](https://github.com/rameyer/SEDmaker/blob/master/docs/functions.md). 
