import numpy as np
import os 

def error(msg):
	print("*(Metadata) >> Error:%s" %msg)
	raise ValueError(msg)

def mag(vector):

	magnitude = np.sqrt(vector[0]**2. + vector[1]**2. + vector[2]**2.)
	return magnitude


def importfile(dirpath, filename, Importance = 'Optional'):
	
	filepath = os.path.join(dirpath, filename)
	if Importance=='Optional' and not os.path.isfile(filepath):
		print("*(metadata) >> Warning: %s file missing, Corresponding data will be taken from parameter file" %filename)
		return filepath
	elif not os.path.isfile(filepath):
		error("{} File not found. Please Check again.".format(filename))
	else:
		return filepath

	
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
	wf_junkrad = "wf_junkrad.txt"		#'/localdata/bkhamesra3/research_localdata/UsefulScripts/LIGO/LIGO_Scripts/LIGO_Metadata/LIGO_Scripts/wf_junkrad.txt'
    
	wfdata = np.genfromtxt(wf_junkrad, dtype=None, comments='#', usecols=(0,1), skip_header=1, delimiter = '\t', names = ('GTID', 'simname'))
    	GTname, wfname = wfdata['GTID'], wfdata['simname']

	if np.array(np.where(wfname==file_name)).size==0:
		error('*(metadata) >> GT simulation name incorrectly listed in wf_junkrad.txt. Please check the file.')
 	elif np.array(np.where(wfname==file_name)).size>1:
		error('*(metadata) >> Multiple occurences of GT simulation - {} in wf_junkrad.txt - Please check the file.'.format(file_name))
	else: 
		idx = (np.where(wfname == file_name))[0][0]
		return GTname[idx]



def simulation_type(spin1, spin2):

	if (np.count_nonzero(spin1) ==0 and np.count_nonzero(spin2)==0): 
		simtype = 'non-spinning'
	elif (np.count_nonzero(spin1[0:2])>0 or np.count_nonzero(spin2[0:2])>0) :
		simtype = 'precessing'
	else:							
		simtype = 'aligned-spins'
	
	return simtype



def updatespins(dirpath, jkrad_time, spin1, spin2):
	
	datadir = os.path.join(dirpath, "data")
	ihspin3 = importfile(datadir,'ihspin_hn_3.asc', Importance = 'Optional')
	ihspin4 = importfile(datadir,'ihspin_hn_4.asc', Importance = 'Optional')

	time_cutoff = jkrad_time
	s1 = np.copy(spin1)
	s2 = np.copy(spin2)

	if not os.path.isfile(ihspin3):
		ihspin0 = importfile(datadir, 'ihspin_hn_0.asc', 'Optional')
		ihspin1 = importfile(datadir, 'ihspin_hn_1.asc', 'Optional')
	else:
		ihspin0 = ihspin3
		ihspin1 = ihspin4

	if os.path.isfile(ihspin0):

		time1_spin = np.loadtxt(ihspin0, usecols = (0,), comments = '#')
		time2_spin = np.loadtxt(ihspin1, usecols = (0,), comments = '#')
		minlen = min(len(time1_spin), len(time2_spin))
		time_spin = time1_spin[:minlen]

		if time_spin[-1]>=jkrad_time:

			temp_arr = np.where(time_spin>=jkrad_time)
			cutoff_idx = np.amin(temp_arr)    
			time_cutoff = time_spin[cutoff_idx] 
			(s1[0],s1[1], s1[2]) = np.loadtxt(ihspin0, usecols= (1,2,3), comments = '#')[cutoff_idx]
			(s2[0],s2[1], s2[2]) = np.loadtxt(ihspin1, usecols= (1,2,3), comments = '#')[cutoff_idx]

		#	print("*(metadata) >> Warning: Interpolated spins absent. Spin information available in ihspin files after junk radiation \n")
			print("*(metadata) >> At time = {}, Spin1 = ({}, {}, {}) and Spin2 = ({}, {}, {}) using ihspins. \n".format(time_cutoff, s1[0], s1[1], s1[2], s2[0], s2[1], s2[2]))
		#		return [s1,s2, time_cutoff]
	else:	
		print("*(Metadata) >> Warning: Spin Information not available in IHSpin files. \n")
	return [s1,s2,time_cutoff]

def updatemass(dirpath, jkrad_time, spin1, spin2):

	datadir = os.path.join(dirpath, "data")
	spin1_mag = mag(spin1)	
	spin2_mag = mag(spin2)		

	bh_diag0 = importfile(datadir,'BH_diagnostics.ah1.gp')
	bh_diag1 = importfile(datadir,'BH_diagnostics.ah2.gp')

	if os.path.exists(bh_diag1):
		time_bh1, irr_m1 = np.loadtxt(bh_diag0, usecols = (1,26), unpack=True, comments = '#')
		time_bh2, irr_m2 = np.loadtxt(bh_diag1, usecols = (1,26), unpack =True, comments = '#')
		minlen = min(len(time_bh1), len(time_bh2))-1
		
		if (time_bh1[minlen]>=jkrad_time):
			cutoff_idx_mass = np.amin(np.where(time_bh1>=jkrad_time))
			m1 = irr_m1[cutoff_idx_mass]	
			m2 = irr_m2[cutoff_idx_mass]
			print("*(Metadata) >> Mass Information provided from BH_Diagnostic files - m1 = {}, m2 = {} at t = {} \n".format(m1, m2, time_bh1[cutoff_idx_mass]))
			return [m1,m2]
		else:
			return [-1.,-1.]


	else: 		
		print("*(metadata) >>  Horizon Mass data unavailable at junk time, t = {} \n".format(jkrad_time))
		return [-1.,-1.]

def updatepos(dirpath, time_cutoff):

	datadir = os.path.join(dirpath, "data")	
	r1, r2, v1, v2 = np.zeros((4,3))	

	shifttracker0 = importfile(datadir, 'ShiftTracker0.asc', Importance='Required')
	shifttracker1 = importfile(datadir, 'ShiftTracker1.asc', Importance='Required')

	time_shtr= np.loadtxt(shifttracker0, usecols = (1,), comments = '#')
	temparr_shtr = np.where(time_shtr>=time_cutoff)
	cutoff_idx_shtr = np.amin(temparr_shtr)     
	
	( r1[0],r1[1],r1[2], v1[0], v1[1], v1[2]) = np.loadtxt(shifttracker0, usecols=( 2, 3, 4, 5, 6, 7), comments = '#')[cutoff_idx_shtr]
	( r2[0],r2[1],r2[2], v2[0], v2[1], v2[2]) = np.loadtxt(shifttracker1, usecols=( 2, 3, 4, 5, 6, 7), comments = '#')[cutoff_idx_shtr]
	
	print("*(metadata) >> At time = {}, r1 = ({}, {}, {}) and r2 = ({}, {}, {}) \n".format(time_shtr[cutoff_idx_shtr], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]))
	return [r1,r2, v1,v2]

