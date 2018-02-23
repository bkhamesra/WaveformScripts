from init_data import initial_data
import numpy as np
import os
from metadata_functions import *
from eccentricity import *
import time


def metadata(dirpath, jkrad_time, Error_Series = False):

# Extract Filename	
	filename = dirpath.split("/")[-1]

#Find junk time in COM frame
	jkrad_time = jkrad_time - 75.
	
# Required Information from Parameter file
	spin1 = np.empty(3)
	spin2 = np.empty(3)
	p1 = np.empty(3)
	p2 = np.empty(3)
	r1 = np.empty(3)
	r2 = np.empty(3)
	v1 = np.empty(3)
	v2 = np.empty(3)
	warning1 = ''
	warning2 = ''

	initdata = initial_data(dirpath)
	spin1 = initdata['spin_BH1']
	spin2 = initdata['spin_BH2']
	p1 = initdata['momentum_BH1']
	p2 = initdata['momentum_BH2']

	try:
		q = float((filename.split("q")[-1]).split("_")[0])
	except ValueError:
		q=1.
	m_plus = q/(1.+q) 
	m_minus = 1./(1.+q)
	eta = m_plus*m_minus/(m_plus+m_minus)**2.

# Find simulation type - nonspinning, align spins or precessing
	simtype = simulation_type(spin1,spin2) 

# Update Position - Always
	[r1,r2,v1,v2] = updatepos(dirpath, jkrad_time)

#Update Spins - Required for precession
	[s1,s2, cutofftime] = updatespins(dirpath, jkrad_time, spin1, spin2)
	
	if np.array_equal(s1,spin1) or np.array_equal(s2,spin2):
		if simtype=='precessing':		
			warning1 = 'Spins Values not found near junk radiation. Precession effects cannot be neglected'
			print("*(Metadata) >> Warning:%s \n "%warning1)
		elif simtype=='aligned-spins':
			warning1 = 'Spin Information not available at t_begin. Will be provided at t=0.'
			print("*(Metadata) >> Warning:%s \n "%warning1)
	else:
		spin1 = s1
		spin2 = s2
		#print("*(Metadata) >> Updated spin1 = {} and spin2 = {} \n".format(spin1, spin2))

# Update Mass - Required for precession
	[m1,m2] = updatemass(dirpath, jkrad_time, spin1, spin2)
	
	if m1==-1. or m2==-1.:
		if simtype=='precessing':
			warning2 = 'Mass Information not found near junk radiation. Will be provided from t=0'
			print ("*(Metadata) >> Warning: %s \n"%warning2)
		elif simtype=='aligned-spins':
			warning2 = 'Mass Information not available at t_begin. Will be provided at t=0.'
			print ("*(Metadata) >> Warning: %s \n"%warning2)
	else:	
		m_plus = m1
		m_minus = m2

# Compute relevant quantities
	delta_r = r2-r1
	init_sep = np.linalg.norm(delta_r)
	nhat = 	delta_r/init_sep
	
	p1 = m_plus*v1		
	p2 = m_minus*v2
	q = m_plus/m_minus
	eta = m_plus*m_minus/(m_plus + m_minus)**2.
	
	orb_angmom = np.cross(r1,p1) + np.cross(r2,p2)
	Lhat = orb_angmom/np.linalg.norm(orb_angmom)

# Computing eccentricity and mean anomaly
	[mean_anomaly, eccentricity] = ecc_and_anomaly(dirpath,  jkrad_time)

	print("*(metadata) >> Final Information \n")
	print("*(metadata) >> Mass Ratio of system = {} \n".format(q))
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
	nr_metadata['NR-techniques'] = 'Puncture-ID, BSSN, Psi4-integrated, Finite-Radius-Waveform, ApproxKillingVector-Spin, Christodoulou-Mass'	

	if Error_Series==True:
		nr_metadata['files-in-error-series'] = '\'\''	#need to check for each simulation
		nr_metadata['comaparable-simulations'] = '\'\''
		print("*(Metadata)>> Please enter the name of simulation and error series manually to metadata.py")
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
	nr_metadata['spin1x'] = round(spin1[0]/m_plus**2, 8)
	nr_metadata['spin1y'] = round(spin1[1]/m_plus**2, 8)
	nr_metadata['spin1z'] = round(spin1[2]/m_plus**2, 8)
	nr_metadata['spin2x'] = round(spin2[0]/m_minus**2, 8)
	nr_metadata['spin2y'] = round(spin2[1]/m_minus**2, 8)
	nr_metadata['spin2z'] = round(spin2[2]/m_minus**2, 8)
	nr_metadata['nhatx'] = round(nhat[0], 8)
	nr_metadata['nhaty'] = round(nhat[1], 8)
	nr_metadata['nhatz'] = round(nhat[2], 8)
	nr_metadata['LNhatx'] = round(Lhat[0], 8)
	nr_metadata['LNhaty'] = round(Lhat[1], 8)
	nr_metadata['LNhatz'] = round(Lhat[2], 8)
	nr_metadata['PN_approximant'] = 'None'
        nr_metadata['eccentricity'] = eccentricity
	nr_metadata['mean_anomaly'] = mean_anomaly
	nr_metadata['Warning1'] = "Please do not use m=0 modes. For caution, they are set to zero!"
	nr_metadata['Warning2'] = warning1
	nr_metadata['Warning3'] = warning2
	return nr_metadata


