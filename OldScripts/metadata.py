'''This script creates metadata dictionary required to generate h5 strain data. This nr data needs to be feeded to the bundlers.py code (by Lionel) with Strain data to generate h5 strain data.
Author - BHAVESH KHAMESRA'''


import numpy as np
import os
import warnings
import time
from eccentricity import ecc_and_anomaly
 
# Outputs numerical values corresponding to specific parameter in file	
def output_data(parfile, data):
	
	datafile = file(parfile)	
	datafile.seek(0)
	for line in datafile:
		if data in line:
			break
	line = line.split()
	data_value = float(line[-1])
	datafile.close()
	return data_value

# Finds the GTID correspondingto specific waveform
def simulation_name(dirpath):
	
	file_name = dirpath.split('/')[-1]
	wf_junkrad = dirpath.split(file_name)[0] + 'wf_junkrad.txt'
    
	wfdata = np.genfromtxt(wf_junkrad, dtype=None, comments='#', usecols=(0,1), skip_header=1, delimiter = '\t', names = ('GTID', 'simname'))
    	GTname, wfname = wfdata['GTID'], wfdata['simname']

	if np.array(np.where(wfname==file_name)).size==0:
		print('*(metadata) >> GT simulation name incorrectly listed in wf_junkrad.txt. Please check the file.')
 	elif np.array(np.where(wfname==file_name)).size>1:
		print('*(metadata) >> Multiple occurences of GT simulation - {} in wf_junkrad.txt - Please check the file.'.format(file_name))
	else: 
		idx = (np.where(wfname == file_name))[0][0]
		return GTname[idx]


def mag(vector):
	magnitude = np.sqrt(vector[0]**2. + vector[1]**2. + vector[2]**2.)
	return magnitude

def metadata(dirpath,  jkrad_time, PNEvo = False, Error_Series = False):	

#Check if the necessary files exists: parfile, shifttracker, ihspin, hnmass
	
	filename = dirpath.split("/")[-1]
    	parfile = (dirpath.split('/')[-1]) + ('.par')
	shifttracker0 = dirpath+'/ShiftTracker0.asc'
	shifttracker1 = dirpath+'/ShiftTracker1.asc'
#	hn_mass0 = dirpath+'/hn_mass_spin_0.asc'
#	hn_mass1 = dirpath+'/hn_mass_spin_1.asc'
	bh_diag1 = dirpath+'/BH_diagnostics.ah1.gp'
	bh_diag2 = dirpath+'/BH_diagnostics.ah2.gp'
	
	if not os.path.isfile(parfile) :
		raise ValueError(" Parfile missing in the directory. Please check again. \n")
	if not os.path.isfile(shifttracker0) or not os.path.isfile(shifttracker1):
		raise ValueError(" ShiftTracker files not found. Please check again. \n")
	
#	hnmass_file = 0
#	if not os.path.isfile(hn_mass0) or not os.path.isfile(hn_mass1):
#		print("*(metadata) >> Warning: hn_mass_spin files not found. \n")
#		hnmass_file = -1

	bhdiag_file=1	
	if not os.path.isfile(bh_diag1) or not os.path.isfile(bh_diag2):
		print("*(metadata) >> Warning: BH_Diagnostic files not found. Mass information will be provided using initial data \n")
		bhdiag_file = -1
		

	ihspin0 = dirpath+'/ihspin_hn_0.asc'
	ihspin1 = dirpath+'/ihspin_hn_1.asc'
	ihspin3 = dirpath+'/ihspin_hn_3.asc'
	ihspin4 = dirpath+'/ihspin_hn_4.asc'
	spinvar = 0

	if not os.path.isfile(ihspin0):
		if os.path.isfile(ihspin3):
			ihspin0 = ihspin3
			ihspin1 = ihspin4
			spinvar = 1
		else:
			print("*(metadata) >> Warning: Spin information not available. Spins taken from initial data in parfile \n")
	else:
		spinvar = 1

# Define masses, mass ratio q and eta - Mass ratio extracted from filename - so ensure that filename has massratio. 

	q = float((filename.split("q")[-1]).split("_")[0])
	m_plus = q/(1.+q) 
	m_minus = 1./(1.+q)
	eta = m_plus*m_minus/(m_plus+m_minus)**2.
	print("*(metadata) >> Mass Ratio of system = {} \n".format(q))

	if m_plus<m_minus:
		print("*(metadata) >> Warning: Mass 1 < Mass 2. Not as per the LIGO convention.")
	else:
	
# Define spins, position, momentum and angular momentum
		r1 = np.empty(3)
		r2 = np.empty(3)
		v1 = np.empty(3)
		v2 = np.empty(3)
		p1 = np.empty(3)
		p2 = np.empty(3)
		spin1 = np.empty(3)
		spin2 = np.empty(3)
		warning = None
			
#Shift Junk radiation in COM frame
		jkrad_time = jkrad_time - 75.
	
# If ihspin files are present, then check if they extend till junk radiation. If they do, then read spin information from ihspin_extrapolated data. If not then use ihspin files
# If ihspin files are absent read data from parameter file
# Comments - changes required - interpolation only useful if ihspin data available for whole simulation but not as frequent as shifttracker. 
		if spinvar==1:	
			time1_spin = np.loadtxt(ihspin0, usecols = (0,), comments = '#')
			time2_spin = np.loadtxt(ihspin1, usecols = (0,), comments = '#')
			minlen = min(len(time1_spin), len(time2_spin))
			time_spin = time1_spin[:minlen]

			if time_spin[-1]>=jkrad_time:
				if os.path.isfile(dirpath+'/Interpolated_Spindata.txt'):  
 					print("*(metadata) >> Spin information available from interpolated data. \n")
					ihspin = dirpath+'/Interpolated_Spindata.txt'			 
					time_intpld = np.loadtxt(ihspin, usecols = (0,), comments = '#')
					temp_arr = np.where(time_intpld>=jkrad_time)
					cutoff_idx = np.amin(temp_arr)     
					time_cutoff = time_intpld[cutoff_idx]
					(spin1[0],spin1[1], spin1[2], spin2[0], spin2[1], spin2[2]) = np.loadtxt(ihspin, usecols= (1,2,3,4,5,6), comments = '#')[cutoff_idx]
					print("*(metadata) >> At time = {}, Spin1 = ({}, {}, {}) and Spin2 = ({}, {}, {}) using interpolation. \n".format(time_cutoff, spin1[0], spin1[1], spin1[2], spin2[0], spin2[1], spin2[2]))
				else :

					print("*(metadata) >> Warning: Interpolated spins absent. Spin information available in ihspin files after junk radiation \n")
					temp_arr = np.where(time_spin>=jkrad_time)
					cutoff_idx = np.amin(temp_arr)    
					time_cutoff = time_spin[cutoff_idx] 
					(spin1[0],spin1[1], spin1[2]) = np.loadtxt(ihspin0, usecols= (1,2,3), comments = '#')[cutoff_idx]
					(spin2[0],spin2[1], spin2[2]) = np.loadtxt(ihspin1, usecols= (1,2,3), comments = '#')[cutoff_idx]
					print("*(metadata) >> At time = {}, Spin1 = ({}, {}, {}) and Spin2 = ({}, {}, {}) using ihspins. \n".format(time_cutoff, spin1[0], spin1[1], spin1[2], spin2[0], spin2[1], spin2[2]))
					
									
			else:

		# If spin information not available in ihspin, then use parameter file. 
				print("*(metadata) >> Warning: Spin information not available in ihspin files. Hence, it will be provided from parameter file. \n")
				spinwarn = 1
				spin1[0] = output_data(parfile, 'par_s_plus[0]')
				spin1[1] = output_data(parfile, 'par_s_plus[1]')
				spin1[2] = output_data(parfile, 'par_s_plus[2]')
				spin2[0] = output_data(parfile, 'par_s_minus[0]')
				spin2[1] = output_data(parfile, 'par_s_minus[1]')
				spin2[2] = output_data(parfile, 'par_s_minus[2]')
				spinvar = 0
				warning = "Spin and Mass information not available at t_begin, will be provided at t=0"


		# Compute Mass Information	
		if bhdiag_file==1 and spinvar==1:
			time_hn1, irr_m1 = np.loadtxt(bh_diag0, usecols = (1,26), comments = '#')
			time_hn2, irr_m2 = np.loadtxt(bh_diag1, usecols = (1,26), comments = '#')
			minlen = min(len(time_hn1), len(time_hn2))-1
		
			if (time_hn0[minlen]>=time_cutoff):
				temparr_hn = np.where(time_hn>=time_cutoff)
				cutoff_idx_mass = np.amin(temparr_hn)     
			else: 		
				#cutoff_idx_mass = minlen
				raise ValueError("*(metadata) >>  Horizon Mass data unavailable at junk time, t = {} \n".format(time_cutoff))

				
			#if time_hn0[cutoff_idx_mass]!=time_cutoff:
			#	print("*(metadata) >>  Warning: Time of cropping (after junk radiation) = {} while mass is computed at t={} \n".format( time_cutoff, time_hn0[cutoff_idx_mass]))
			
			spin1_mag = mag(spin1)		 #np.loadtxt(hn_mass0, usecols = (1,), comments='#')[cutoff_idx_mass]
			spin2_mag = mag(spin2)		#np.loadtxt(hn_mass1, usecols = (1,), comments='#')[cutoff_idx_mass]
				

			m1 = np.sqrt(irr_m1[cutoff_idx_mass]**2. + spin1_mag**2./(4.*irr_m1[cutoff_idx_mass]**2.))
			m1_err = (m1-m_plus)*100./m_plus
				
			m2 = np.sqrt(irr_m2[cutoff_idx_mass]**2. + spin2_mag**2./(4.*irr_m2[cutoff_idx_mass]**2.))
			m2_err = (m2-m_minus)*100./m_minus
			print("*(metadata) >> Comparing with initial data, Change in mass 1 = {}% and Change in mass 2 = {}% \n".format(m1_err, m2_err))
	
			m_plus = m1
			m_minus = m2
			masswarn = 0
		else: 
			print("*(metadata) >> Mass Information not available in BH_diagnostic files. Information will be taken from parameter files. \n")
			masswarn =1

		# Find the spin configuration
		if (np.count_nonzero(spin1) ==0 and np.count_nonzero(spin2)==0): 
			simtype = 'non-spinning'
			spinvar = 1
			time_cutoff = jkrad_time
			warning = None
		elif (np.count_nonzero(spin1[0:2])>0 or np.count_nonzero(spin2[0:2])>0) :
			simtype = 'precessing'
		
		else:							
			simtype = 'aligned-spins'
			spinvar = 1
			time_cutoff = jkrad_time
		
		# Error if precessing Waveforms lack spin information 
		if simtype=='precessing' and (spinvar==0 or masswarn==1):
			raise ValueError('Spin Information not Available for this precessing waveform. Metadata information cannot be computed.')
		
		# Warning:
		if simtype=='aligned-spins':
			if spinwarn==1 and masswarn==1:
				warning = 'Spin and Mass Information unavailable at t_begin, will be provided at t=0'	
			elif  spinwarn==1:
				warning = 'Spin Information unavailable at t_begin, will be provided at t=0'
				
			elif  masswarn==0:
				warning = 'Mass Information unavailable at t_begin, will be provided at t=0.'

		# Obtain the position and momentum data
		if spinvar==1:
			time_shtr= np.loadtxt(shifttracker0, usecols = (1,), comments = '#')
			temparr_shtr = np.where(time_shtr>=time_cutoff)
			cutoff_idx_shtr = np.amin(temparr_shtr)     
			if time_shtr[cutoff_idx_shtr]!=time_cutoff:
				print("*(metadata) >>  Warning: Time of cropping (after junk radiation) t = {} while position and momentum computed at t={} \n".format( time_cutoff, time_shtr[cutoff_idx_shtr]))
			
			( r1[0],r1[1],r1[2], v1[0], v1[1], v1[2]) = np.loadtxt(shifttracker0, usecols=( 2, 3, 4, 5, 6, 7), comments = '#')[cutoff_idx_shtr]
			( r2[0],r2[1],r2[2], v2[0], v2[1], v2[2]) = np.loadtxt(shifttracker1, usecols=( 2, 3, 4, 5, 6, 7), comments = '#')[cutoff_idx_shtr]
			
			print("*(metadata) >> At time = {}, r1 = ({}, {}, {}) and r2 = ({}, {}, {}) \n".format(time_shtr[cutoff_idx_shtr], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]))
			
			
			p1 = m_plus*v1		#check with Karan
			p2 = m_minus*v2
			q = m_plus/m_minus
			eta = m_plus*m_minus/(m_plus + m_minus)**2.
# If using Post Newtonian, modify this part also for position and momentum.
		else:
			print ("*(metadata) >> Position and Momentum data provided from parameter file")
			r1[0] = output_data(parfile, 'x0[0]')
			r1[1] = output_data(parfile, 'y0[0]')
			r1[2] = output_data(parfile, 'z0[0]')
			r2[0] = output_data(parfile, 'x0[1]')
			r2[1] = output_data(parfile, 'y0[1]')
			r2[2] = output_data(parfile, 'z0[1]')
		
			p1[0] = output_data(parfile, 'par_P_plus[0]')
			p1[1] = output_data(parfile, 'par_P_plus[1]')
			p1[2] = output_data(parfile, 'par_P_plus[2]')
			p2[0] = output_data(parfile, 'par_P_minus[0]')
			p2[1] = output_data(parfile, 'par_P_minus[1]')
			p2[2] = output_data(parfile, 'par_P_minus[2]')
		

		delta_r = r2-r1
		init_sep = np.linalg.norm(delta_r)
		nhat = 	delta_r/init_sep	

		
		print("*(metadata) >> Initial separation  = {} and nhat = {} \n".format(init_sep, nhat))
		print("*(metadata) >> Momentum of BH1 = {}, and Momentum of BH2 = {} \n".format( p1, p2))
		
		print("*(metadata) >> Spin of BH1 = {} and Spin of BH2={} \n".format( spin1, spin2))
		print("*(metadata) >> Masses of BH: mplus = {}, m_minus = {} and q = {} \n". format(m_plus,m_minus, q ))

		orb_angmom = np.cross(r1,p1) + np.cross(r2,p2)
		Lhat = orb_angmom/np.linalg.norm(orb_angmom)#Check with Juan
		print("*(metadata) >> Orbital Angular Momentum vector is = {} \n".format( Lhat))

 		#Computing eccentricity and mean anomaly
		[mean_anomaly, eccentricity] = ecc_and_anomaly(dirpath,  jkrad_time)
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
		nr_metadata['NR-techniques'] = 'Puncture-ID, BSSN, Psi4-integrated, Finite-Radius-Waveform, ApproxKillingVector-Spin, Christodoulou-Mass'	#check with Pablo/Deirdre
		
		if Error_Series==False:
			nr_metadata['files-in-error-series'] = '\'\''	#need to check for each simulation
		else: nr_metadata['files-in-error-series'] = -1
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

		nr_metadata['Warning1'] = warning
		return nr_metadata


