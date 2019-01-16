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



def updatespins(dirpath, retarted_junktime, spin1, spin2):
	
	datadir = os.path.join(dirpath, "data")
	ihspin_bh1 = importfile(datadir,'ihspin_hn_0.asc', Importance = 'Optional')
	ihspin_bh2 = importfile(datadir,'ihspin_hn_1.asc', Importance = 'Optional')

	simtype = simulation_type(spin1, spin2)
	s1 = np.copy(spin1)
	s2 = np.copy(spin2)


	if not os.path.isfile(ihspin_bh1) or not os.path.isfile(ihspin_bh2):
		ihspin0 = importfile(datadir, 'ihspin_hn_3.asc', 'Optional')
		ihspin1 = importfile(datadir, 'ihspin_hn_4.asc', 'Optional')
	else:
		ihspin0 = ihspin_bh1
		ihspin1 = ihspin_bh2

	if os.path.isfile(ihspin0) and os.path.isfile(ihspin1):

		time1_spin = np.loadtxt(ihspin0, usecols = (0,), comments = '#')
		time2_spin = np.loadtxt(ihspin1, usecols = (0,), comments = '#')
		#minlen = min(len(time1_spin), len(time2_spin))
		time_spin = min(time1_spin[-1], time2_spin[-1])	#time1_spin[:minlen]

		if time_spin[-1]>=retarted_junktime:

			t1_idx = np.amin(np.where(time1_spin>=retarted_junktime))
			t2_idx = np.amin(np.where(time2_spin>=retarted_junktime))

			time_cutoff = time1_spin[t1_idx] 
			(s1[0],s1[1], s1[2]) = np.loadtxt(ihspin0, usecols= (1,2,3), comments = '#')[t1_idx]
			(s2[0],s2[1], s2[2]) = np.loadtxt(ihspin1, usecols= (1,2,3), comments = '#')[t2_idx]

			print("*(Metadata) >> At time = {}, Spin1 = ({}, {}, {}) and Spin2 = ({}, {}, {}) using ihspins. \n".format(time_cutoff, s1[0], s1[1], s1[2], s2[0], s2[1], s2[2]))

		elif not(simtype=='precessing')

			time_cutoff = time1_spin[-1]			

			(s1[0],s1[1], s1[2]) = np.loadtxt(ihspin0, usecols= (1,2,3), comments = '#')[-1]
			(s2[0],s2[1], s2[2]) = np.loadtxt(ihspin1, usecols= (1,2,3), comments = '#')[-1]
			
			print("*(Metadata) >> Warning: Spin Information not available in after retarted junktime. Last available information will be taken from ihspin files\n")
			print("*(Metadata) >> At time = {}, Spin1 = ({}, {}, {}) and Spin2 = ({}, {}, {}) using ihspins. \n".format(time_cutoff, s1[0], s1[1], s1[2], s2[0], s2[1], s2[2]))
		else:
			raise ValueError('Spin Information not available for %s simulation. Metadata cannot be generated \n')

		return [s1,s2,time_cutoff]
	else:	
		if not (simtype=='precessing'):
			print("*(Metadata) >> Warning: Spin Information not available in IHSpin files. Spin Information will be taken from parfile\n")
		    	return [spin1, spin2, 0]
		else:
			raise ValueError('*(Metadata) >> Spin Function: Spin Information not found for precessing simulation. Metadata cannot be created \n')


def updatemass(dirpath, retarted_junktime):

	datadir = os.path.join(dirpath, "data")

	bh_diag0 = importfile(datadir,'BH_diagnostics.ah1.gp')
	bh_diag1 = importfile(datadir,'BH_diagnostics.ah2.gp')

	if os.path.exists(bh_diag1):
		time_bh1, irr_m1 = np.loadtxt(bh_diag0, usecols = (1,26), unpack=True, comments = '#')
		time_bh2, irr_m2 = np.loadtxt(bh_diag1, usecols = (1,26), unpack =True, comments = '#')
		#minlen = min(len(time_bh1), len(time_bh2))-1
		time_mass = min(time_bh1[-1], time_bh2[-1])	#time_bh1[:minlen]
		
		if (time_mass[-1]>=retarted_junktime):

			t1_idx = np.amin(np.where(time_bh1>=retarted_junktime))
			t2_idx = np.amin(np.where(time_bh2>=retarted_junktime))
			time_cutoff = time_bh1[t1_idx]
			m1 = irr_m1[cutoff_idx_mass]	
			m2 = irr_m2[cutoff_idx_mass]
		else:
			m1 = irr_m1[-1]
			m2 = irr_m2[-1]
			time_cutoff = min(time_bh1[-1], time_bh2[-1])
			print("*(Metadata) >> Warning: Mass information not available after junk radiation. Mass will be provided at t = {} \n".format(time_cutoff))

		print("*(Metadata) >> Mass Information provided from BH_Diagnostic files - m1 = {}, m2 = {} at t = {} \n".format(m1, m2, time_bh1[cutoff_idx_mass]))
		return [m1,m2, time_cutoff]

	else: 		
		raise ValuesError("*(Metadata) >> Mass Function - BH_diagnostics files not found. Metadata cannot be created \n")
		

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


def func_phase(varphase):

	varphi = np.copy(varphase)
	for i in range(len(varphase)):
		if abs(varphase[i-1]-varphase[i]-2.*np.pi)<0.1:
			varphi[i:] = varphi[i:] + 2.*np.pi
	return varphi


def write_sep_data(filename,hdr, outdir, data):
	output_traj = open(os.path.join(outdir, filename),'w')
	np.savetxt(output_traj, data, header=hdr, delimiter='\t', newline='\n')
	output_traj.close()


def find_omega22_st(dirpath, time_cutoff):
	
	datadir = os.path.join(dirpath, "data")	

	shifttracker0 = importfile(datadir, 'ShiftTracker0.asc', Importance='Required')
	shifttracker1 = importfile(datadir, 'ShiftTracker1.asc', Importance='Required')

	time_shtr = np.loadtxt(shifttracker0, usecols = (1), comments = '#')
	temparr_shtr = np.intersect1d(np.where(time_shtr>=time_cutoff),  np.where(time_shtr<=150))

	time_shtr1, r1x, r1y, r1z = np.loadtxt(shifttracker0, usecols = (1,2,3,4), comments = '#', unpack=True)
	time_shtr2, r2x, r2y, r2z = np.loadtxt(shifttracker1, usecols = (1,2,3,4), comments = '#', unpack=True)
	
	x = r1x - r2x
	y = r1y - r2y

	phase = np.arctan2(y, x)
	
	phi = func_phase(phase[temparr_shtr])
        t = time_shtr1[temparr_shtr]
	omega = np.gradient(phi)/np.gradient(t)	

	return omega
   	

def update_orbital_frequency(dirpath, time_cutoff):
	
	datadir = os.path.join(dirpath, "data")	

	shifttracker0 = importfile(datadir, 'ShiftTracker0.asc', Importance='Required')
	shifttracker1 = importfile(datadir, 'ShiftTracker1.asc', Importance='Required')

	time_shtr1, r1x, r1y, r1z, v1x, v1y, v1z = np.loadtxt(shifttracker0, usecols = (1,2,3,4,5,6,7), comments = '#', unpack=True)
	time_shtr2, r2x, r2y, r2z, v2x, v2y, v2z = np.loadtxt(shifttracker1, usecols = (1,2,3,4,5,6,7), comments = '#', unpack=True)
	
	sepx = r1x - r2x
	sepy = r1y - r2y
	sepz = r1z - r2z

	sep = np.array((sepx, sepy, sepz))
	nhat = sep/mag(sep)

	dv = np.array(v1x-v2x, v1y-v2y, v1z-v2z))
	dnhat_dt = dv/mag(sep)
	
	omega_orbfreq = np.cross(nhat, dnhat_dt)
		
	cutoff_idx_shtr = np.amin(np.where(time_shtr>=time_cutoff))
	omega_orbfreq = omega_obr[cutoff_idx_shtr]
	return omega_orbfreq
   	
