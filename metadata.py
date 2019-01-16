#Comments - Check if spin data available till junk -75M for precessing cases and if not, then do not build the waveform. 


from init_data import initial_data
import numpy as np
import os
from metadata_functions import *
from eccentricity import *
import time


def metadata(dirpath, jkrad_time, Error_Series = False):

	print("\n \n *(metadata) >> Extracting Metadata - Analysing Data \n")
	# Extract Filename	
	filename = dirpath.split("/")[-1]

	#Find junk time in COM frame
	retarted_junktime = jkrad_time - 75.
	
	
	# Required Information from Parameter file
	spin1 = np.empty(3)
	spin2 = np.empty(3)
	p1 = np.empty(3)
	p2 = np.empty(3)
	r1 = np.empty(3)
	r2 = np.empty(3)
	v1 = np.empty(3)
	v2 = np.empty(3)

	#Extract Initial parameters from parfile
	initdata = initial_data(dirpath)
	spin1 = initdata['spin_BH1']
	spin2 = initdata['spin_BH2']
	p1 = initdata['momentum_BH1']
	p2 = initdata['momentum_BH2']

	#Check if initial system is quasicircular
	qc_sys=False
	if (p1[0]==-1*p2[0]) and (p1[1]==-1.*p2[1]):
	        qc_sys = True
	

	try:
		q = float((filename.split("q")[-1]).split("_")[0])
	except ValueError:
		m1 = initdata['puncture_mass1']
		m2 = initdata['puncture_mass2']
		q = np.around(m1/m2, decimal=1)

	m_plus = q/(1.+q) 
	m_minus = 1./(1.+q)
	#eta = m_plus*m_minus/(m_plus+m_minus)**2.


	# Find simulation type - nonspinning, align spins or precessing
	simtype = simulation_type(spin1,spin2) 

	
	# Update Position - Always
	[r1,r2,v1,v2] = updatepos(dirpath, retarted_junktime)

	
	#Update Spins - Required for precession
	[s1,s2, spin_cutofftime] = updatespins(dirpath, retarted_junktime, spin1, spin2)
	
	spin1 = s1
	spin2 = s2


	# Update Mass - 
	m1, m2, mass_cutofftime = updatemass(dirpath, retarted_junktime)
	m_plus = m1
	m_minus = m2
	
	
	#Orbital Frequency - 
	orbital_freq = update_orbital_frequency(dirpath, retarted_junktime)

	orbital_freq_phasecomp = find_omega22_st(dirpath, retarted_junktime)

	print("*(metadata): Orbital Frequency using rxv formula = %g and using orbital phase = %g \n"%(orbital_freq, orbital_freq_phasecomp))
	# Compute relevant quantities
	# Separation Unit Vector
	delta_r = r2-r1
	init_sep = np.linalg.norm(delta_r)
	nhat = 	delta_r/init_sep
	
	assert(mag(nhat)==1.), "Error: Norm of nhat vector = %g (not normalized)"%mag(nhat)	

	#Mass Ratio
	q = m_plus/m_minus
	eta = m_plus*m_minus/(m_plus + m_minus)**2.
	
	#Newtonian Orbital Angular Momentum Unit Vector
	p1 = m_plus*v1		
	p2 = m_minus*v2
	
	orb_angmom = np.cross(r1,p1) + np.cross(r2,p2)
	Lhat = orb_angmom/np.linalg.norm(orb_angmom)
	assert(mag(Lhat)==1.), "Error: Norm of Lhat vector = %g (not normalized)"%mag(Lhat)	
	
	#Horizon Mass
	s1_mag = np.sqrt(s1[0]**2 + s1[1]**2 + s1[2]**2)
	s2_mag = np.sqrt(s2[0]**2 + s2[1]**2 + s2[2]**2)

	horizonmass_plus = np.sqrt(m_plus**2 + 0.25*s1_mag**2/m_plus**2)
	horizonmass_minus = np.sqrt(m_minus**2 + 0.25*s1_mag**2/m_minus**2)


	#Dimensionless Spins
	a1x = s1[0]/horizonmass_plus**2
	a1y = s1[1]/horizonmass_plus**2
	a1z = s1[2]/horizonmass_plus**2
	a1 = np.array((a1x, a1y, a1z))

	a2x = s2[0]/horizonmass_minus**2
	a2y = s2[1]/horizonmass_minus**2
	a2z = s2[2]/horizonmass_minus**2
	a2 = np.array((a2x, a2y, a2z))

	assert(mag(a1)<=1), "Error: |a1| = %g. Dimensionless spin magnitude cannot be greater than 1!"%mag(a1)
	assert(mag(a2)<=1), "Error: |a2| = %g. Dimensionless spin magnitude cannot be greater than 1!"%mag(a2)
	# Computing eccentricity and mean anomaly
	[mean_anomaly, eccentricity] = ecc_and_anomaly(dirpath,  jkrad_time)
	
	if qc_sys and ecc>1.e-3:
		print("*(Metadata) >> Warning: Eccentricity of system too high (%g) for being quasicircular \n"%eccentricity)

	if np.isnan(mean_anomaly) or np.isinf(mean_anomaly):
		raise NameError("Mean Anomaly Value not correct(%g). Please check the computation method"%mean_anomaly)

	print("\n \n *(metadata) >> Final Results \n")
	print("*(metadata) >> Mass Ratio of system, Symmetric Mass Ratio = {} \n".format(q))
	print("*(metadata) >> Initial separation  = {} and nhat = {} \n".format(init_sep, nhat))
	print("*(metadata) >> Momentum of BH1 = {}, and Momentum of BH2 = {} \n".format( p1, p2))
	
	print("*(metadata) >> Spin of BH1 = {} and Spin of BH2={} \n".format( spin1, spin2))
	print("*(metadata) >> Masses of BH: mplus = {}, m_minus = {} and q = {} \n". format(m_plus,m_minus, q ))
	print("*(metadata) >> Orbital Angular Momentum vector is = {} \n".format( Lhat))
	print("*(metadata) >> Mean Anomaly =  {} and Eccentricity = {} \n".format( mean_anomaly, eccentricity))
	
	if simtype=='precessing': comp_sim = 'GT0560'
	else: comp_sim = 'GT0582'

	nr_metadata = {}
	nr_metadata['Format'] = 1
	nr_metadata['NR-group'] = 'Georgia Tech'
	nr_metadata['NR-code'] = 'MAYA'
	nr_metadata['type'] = 'NRinjection'		
	nr_metadata['simulation-type'] = simtype
	nr_metadata['name'] = simulation_name(dirpath)
	nr_metadata['alternative-names'] = filename
	nr_metadata['modification-date'] = time.strftime("%Y-%m-%d")
	nr_metadata['Contact Person'] = 'Deirdre Shoemaker'
	nr_metadata['point-of-contact-email'] = 'deirdre.shoemaker@physics.gatech.edu'
	nr_metadata['license'] = 'public'
	nr_metadata['INSPIRE-bibtex-keys'] = 'Jani:2016wkt'	
	nr_metadata['NR-techniques'] = 'Puncture-ID, BSSN, Psi4-integrated, Extrapolated-Waveform, ApproxKillingVector-Spin, Christodoulou-Mass'	

	if Error_Series==True:
		nr_metadata['files-in-error-series'] = '\'\''	#need to check for each simulation
		nr_metadata['comaparable-simulations'] = '\'\''
		print("*(Metadata)>> Please update the code with list of error series for this simulation. ")
	else:
		nr_metadata['files-in-error-series'] = ''
		nr_metadata['comparable-simulation'] = comp_sim	#need to check for each simulation

	nr_metadata['production-run'] = 1		#set to 0 for lower resolution runs
	nr_metadata['object1'] = 'BH'
	nr_metadata['object2'] = 'BH'
	nr_metadata['init_sep'] = init_sep
	nr_metadata['mass1'] = round(m_plus, 8)
	nr_metadata['mass2'] = round(m_minus, 8)
	nr_metadata['eta'] = round(eta, 8)
	nr_metadata['spin1x'] = round(a1x, 8)
	nr_metadata['spin1y'] = round(a1y, 8)
	nr_metadata['spin1z'] = round(a1z, 8)
	nr_metadata['spin2x'] = round(a2x, 8)
	nr_metadata['spin2y'] = round(a2y, 8)
	nr_metadata['spin2z'] = round(a2z, 8)
	nr_metadata['nhatx'] = round(nhat[0], 8)
	nr_metadata['nhaty'] = round(nhat[1], 8)
	nr_metadata['nhatz'] = round(nhat[2], 8)
	nr_metadata['LNhatx'] = round(Lhat[0], 8)
	nr_metadata['LNhaty'] = round(Lhat[1], 8)
	nr_metadata['LNhatz'] = round(Lhat[2], 8)
	nr_metadata['PN_approximant'] = 'None'
        nr_metadata['eccentricity'] = eccentricity
	nr_metadata['mean_anomaly'] = mean_anomaly
	nr_metadata['Omega'] = orbital_freq
	nr_metadata['Warning1'] = "Please do not use m=0 modes. For caution, they are set to zero!"


	return nr_metadata


