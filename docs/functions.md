SEDmaker : SFH available
======

- Constant SFH

```
SEDmaker.SEDmaker.SFH_constant(age,spectrafile,SFR =1):
	'''
	Compute the SED from a constant Star Formation History and a BPASS SED file

	Inputs:
		SFR: Star Formation Rate [M_sun yr^{-1}] (normalized to 1, absolute amplitude 
											does not impact SED fitting)
		age: Time of the constant SF [Myr]
		spectrafile: Nx52 array containing the BPASS SEDs
	Return:
		SED: 1-100000 Angstrom flux in a 1-N array.
		SFR_array: SFR for each timestep, effectivly the SFH
	'''
```

- Linear SFH



```
SEDmaker.SEDmaker.SFH_linear(age,slope,spectrafile,SFR = 1):
	'''
	Compute the SED from a linear SFH and a BPASS SED file. 

	Inputs:
		age: duration of the SF [Myr]
		SFR: Initial SFR (normalized to 1 by default)
		slope: slope of the linear SFH [M_sun yr^{-1} Myr^{-1}]
		spectrafile:  Nx52 array containing the BPASS SEDs
	Return 
		SED: 1-100000 Angstrom flux in a 1-N array.
		SFR_array: SFR for each timestep, effectivly the SFH
	'''
```


- Exponential SFH

```
SEDmaker.SEDmaker.SFH_exponential(age,tau,spectrafile,SFR=1):
	'''
	Compute the SED from a linear SFH and a BPASS SED file. 

	Inputs:
		age: duration of the SF [Myr]
		SFR: Initial SFR (normalized to 1 by default)
		tau: e-folding time of the exponential SFH [Myr] 
			(please specify the sign for decreasing (-) or increasing (+) SFR)
		spectrafile:  Nx52 array containing the BPASS SEDs
	Return 
		SED: 1-100000 Angstrom flux in a 1-N array.
		SFR_array: SFR for each timestep, effectivly the SFH
	'''
```