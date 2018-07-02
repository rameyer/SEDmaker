import numpy as np
import matplotlib.pyplot as plt
import pkg_resources

### Kroupa IMF (imf135_300) slope 0.1-0.5 : -1.3 , 0.5-300: -2.30 . Extends to 300 M_sol

PACKAGE_PATH = str(np.genfromtxt(pkg_resources.resource_filename('SEDmaker', 'data/index_data.txt'),dtype='str'))

IMF = 'imf135_300'

binaries = True

def _agebins_BPASS():
	'''
	Private function returning the logarithmic agebins of the BPASS models

	Inputs:
		-- 
	Return:
		dt_array: 1x52 array containing the agebins dt in [yr]
	'''

	dt_array = np.zeros(52)

	for i  in range(52):
		if i == 0:
			dt_array[i] = 10**6.05
		else:
			dt_array[i] = 10**(6.15+0.1*i) - 10**(6.05+0.1*i)

	return dt_array

def _time_midsteps_BPASS():
	return 10**np.linspace(6,11.1,52)/1e6


def _compute_SED(SFR_array,spectrafile):
	'''
	A private function computing the SED from the BPASS spectra and the SFR given
	as an array of values for each 52 timestep corresponding to the logscale of the 
	BPASS models.

	Inputs:
		SFR_array: 1x52-array containing the SFR at each epoch.
		spectrafile: Nx52 array containing the BPASS SEDs
	Return:
		SED: 1-100000 Angstrom flux in a 1-N array.
	'''
	assert len(SFR_array)==52,'Please ensure the SFR is given for every BPASS agebins'
	assert len(spectrafile[0])==52,'Please ensure the SED are given for every BPASS agebins'
	
	agebins_yr = _agebins_BPASS()*1e6

	return np.sum(spectrafile*agebins_yr*SFR_array, axis = 1)

	'''
	for i  in range(52):
		if i == 0:
			dt = 10**6.05
			SED = spectrafile[:,i+1]*dt*SFR_array[i] 
		else:
			dt = 10**(6.15+0.1*i) - 10**(6.05+0.1*i)
			SED += spectrafile[:,i+1]*dt*SFR_array[i]

	return SED
	'''

def SFH_constant(age,spectrafile,SFR =1):
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
	i_max = int(np.ceil(np.log10(age) / 0.1))
	assert i_max < 52, 'Please use an SFH with duration less than ~100 Gyr.'
	assert SFR > 0, 'Please use a positive SFR' 
	assert age > 0, 'Please use a positive age' 


	SFR_array = np.ones(52)*SFR
	SFR_array = SFR_array * (np.arange(52)<i_max)

	return _compute_SED(SFR_array=SFR_array,spectrafile=spectrafile), SFR_array

def SFH_linear(age,slope,spectrafile,SFR = 1):
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

	assert int(np.ceil(np.log10(age) / 0.1)), 'Please use an SFH with duration less than ~100 Gyr.'
	assert SFR > 0, 'Please use a positive SFR' 
	assert age > 0, 'Please use a positive age' 

	t =  age - _time_midsteps_BPASS()	 # Myr

	SFR_array = SFR + slope * t 
	SFR_array = SFR_array  * (t>0) * (SFR_array>0)

	return _compute_SED(SFR_array=SFR_array, spectrafile=spectrafile), SFR_array

def SFH_exponential(age,tau,spectrafile,SFR=1):
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

	assert int(np.ceil(np.log10(age) / 0.1)), 'Please use an SFH with duration less than ~100 Gyr.'
	assert SFR > 0, 'Please use a positive SFR' 
	assert age > 0, 'Please use a positive age' 	
	
	t =  age - _time_midsteps_BPASS()	 # Myr

	SFR_array = SFR * np.exp( t / tau )
	SFR_array = SFR_array * (t>0) * (SFR_array>0)

	return _compute_SED(SFR_array=SFR_array, spectrafile=spectrafile), SFR_array

def test_SFHs():
	'''
	A private test function to show the range of SFHs available with this code. 
	Inputs:
		--
	Outputs:
		--
	'''

	z_met = '001'
	
	if binaries: 
		spectra_single = np.loadtxt(PACKAGE_PATH+'/BPASS/BPASSv2.1_'+ IMF +'/spectra-bin.z' + z_met 
										+ '.dat')[0:100000:20]
	else:
		spectra_single = np.loadtxt(PACKAGE_PATH+'/BPASS/BPASSv2.1_'+ IMF +'/spectra.z' + z_met 
										+ '.dat')[0:100000:20]

	_ , SFH_1 = SFH_constant(30, spectrafile= spectra_single)
	_ , SFH_2 = SFH_constant(10, spectrafile= spectra_single, SFR = 10)
	_ , SFH_3 = SFH_linear(age = 50,slope = 2, spectrafile= spectra_single,SFR = 1)
	_ , SFH_4 = SFH_linear(age = 10,slope = -1, spectrafile= spectra_single,SFR = 10)
	_ , SFH_5 = SFH_exponential(age = 20,tau = 10, spectrafile= spectra_single,SFR=1)
	_ , SFH_6 = SFH_exponential(age = 20,tau = -5, spectrafile= spectra_single, SFR = 1)

	age_bins = _time_midsteps_BPASS()	# Myr

	plt.plot(age_bins,SFH_1,label = 'Constant 30 Myr, SFR = 1')
	plt.plot(age_bins,SFH_2,label = 'Constant 10 Myr, SFR = 10')
	plt.plot(age_bins,SFH_3,label = 'Linear 50 Myr, Slope = 2')
	plt.plot(age_bins,SFH_4,label = 'Linear 10 Myr, Slope = -1')
	plt.plot(age_bins,SFH_5,label = r'Exponential 20 Myr, $\tau$ = 10')
	plt.plot(age_bins,SFH_6,label = r'Exponential 20 Myr, $\tau$ = -5')
	plt.legend()
	plt.xlim(1,100)
	plt.ylim(0,20)
	if binaries:
		plt.title('BPASS SFH, binaries included')
	else:
		plt.title('BPASS SFH, NO binaries')
	plt.ylabel(r'SFR [$M_\odot$ yr $^{-1}$]')
	plt.xlabel(r'Star populations ages [Myr]')
	plt.show()


def produce_grid_SED_HyperZ():
	''' 
	An example of a grid production of SEDs in HyperZ specific format. 
	Still under work to be fully determined from inputs.
	Inputs:
		--
	Outputs:
		-- (save the SEDs as text files in given directory)
	'''
	f = open(PACKAGE_PATH+'/SED/BPASS_imf135_300_expi.param','w+')
	f_nobin = open(PACKAGE_PATH+'/SED/BPASS_imf135_300_expi_nobin.param','w+')
	f_bin = open(PACKAGE_PATH+'/SED/BPASS_imf135_300_expi_bin.param','w+')
	for z_met in ['em5','em4','001','002','003','004','006','010','014','020','030','040']:
		for age  in [10,20,50,100,200,400]:
			for tau in [0.03,0.06,0.1,0.2,0.4,0.8,1.5,3,5,10]: 

				spectra_single = np.loadtxt(PACKAGE_PATH+'/BPASS/BPASSv2.1_'+ IMF +'/spectra.z' + z_met 
											+ '.dat')[0:6000]
				spectra_bin = np.loadtxt(PACKAGE_PATH+'/BPASS/BPASSv2.1_'+ IMF +'/spectra-bin.z'+ z_met 
											+ '.dat')[0:6000]

				SED_single,_ = SFH_exponential(age=age,tau = tau*1e3,spectrafile=spectra_single)
				SED_bin,_ = SFH_exponential(age=age,tau = tau*1e3,spectrafile=spectra_bin)

				wavelength = np.linspace(1,1e5,1e5)[0:6000]

				# Normalization (if you wish)
				F_lambda_single = SED_single #/ (SED_single[2000])
				F_lambda_bin = SED_bin  #/ (SED_bin[2000])
								
				np.savetxt(PACKAGE_PATH+'/SED/imf135_300/SFH_exp_a' + str(age) + '_t' + str(tau)+ '_z' + z_met + '.sed', 
						   np.transpose(np.vstack((wavelength,F_lambda_single))))
				
				np.savetxt(PACKAGE_PATH+'/SED/imf135_300/SFH_exp_a' + str(age) + '_t' + str(tau)+ '_z' + z_met + '_bin.sed', 
						   np.transpose(np.vstack((wavelength,F_lambda_bin))))

				f.write(PACKAGE_PATH + '/SED/imf135_300/SFH_exp_a' + str(age) + '_t' + str(tau) +  '_z' 
					    + z_met + '.sed   AS  \n')
				f_bin.write(PACKAGE_PATH+'/SED/imf135_300/SFH_exp_a' + str(age) + '_t' + str(tau)+ '_z' 
					    + z_met + '_bin.sed   AS  \n')
				f_nobin.write(PACKAGE_PATH+'/SED/imf135_300/SFH_exp_a' + str(age) + '_t' + str(tau)+ '_z' 
					    + z_met + '.sed   AS  \n')
				f.write(PACKAGE_PATH+'/SED/imf135_300/SFH_exp_a' + str(age)+ '_t' + str(tau) + '_z' 
					    + z_met + '_bin.sed   AS \n')


