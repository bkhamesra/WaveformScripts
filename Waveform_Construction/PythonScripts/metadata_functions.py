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
	line = line.split('=')[1]
	data_value = float(line.split()[0])
	datafile.close()
	return data_value


# Finds the GTID correspondingto specific waveform
def simulation_name(dirpath):
	
	file_name = dirpath.split('/')[-1]
	wf_junkrad = "wf_junkrad.txt"		#'/localdata/bkhamesra3/research_localdata/UsefulScripts/LIGO/LIGO_Scripts/LIGO_Metadata/LIGO_Scripts/wf_junkrad.txt'
	wfdata = np.genfromtxt(wf_junkrad, dtype=None, comments='#', usecols=(0,1), delimiter = '\t', names = ('GTID', 'simname'))
    	
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



def updatespins(dirpath, retarted_junktime, spin1, spin2, verbose=True):
	
	datadir = os.path.join(dirpath, "data")
	ihspin0_bh1 = os.path.join(datadir,'ihspin_hn_0.asc')
	ihspin1_bh2 = os.path.join(datadir,'ihspin_hn_1.asc')

	ihspin3_bh1 = os.path.join(datadir,'ihspin_hn_3.asc')
	ihspin4_bh2 = os.path.join(datadir,'ihspin_hn_4.asc')

	
	simtype = simulation_type(spin1, spin2)
	s1 = np.copy(spin1)
	s2 = np.copy(spin2)

	spindata_from_hznfndr = True
	spindata_from_sphrad = True
	
	if not os.path.isfile(ihspin3_bh1) or not os.path.isfile(ihspin4_bh2):
		spindata_from_sphrad = False
	
	if not os.path.isfile(ihspin0_bh1) or not os.path.isfile(ihspin1_bh2):
		spindata_from_hznfndr = False
	
	if spindata_from_hznfndr and spindata_from_sphrad:
		time_ih0 =  np.loadtxt(ihspin0_bh1, usecols = (0,), comments = '#')
		time_ih1 =  np.loadtxt(ihspin1_bh2, usecols = (0,), comments = '#')
		time_ih3 =  np.loadtxt(ihspin3_bh1, usecols = (0,), comments = '#')
		time_ih4 =  np.loadtxt(ihspin4_bh2, usecols = (0,), comments = '#')
	
		tlast_ih01 = min(time_ih0[-1], time_ih1[-1])
		tlast_ih34 = min(time_ih3[-1], time_ih4[-1])
	
		choose_ih01 = (tlast_ih01 >= retarted_junktime)
		choose_ih34 = (tlast_ih34 >= retarted_junktime)
	
		if choose_ih01: 
			ihspin0 = importfile(datadir, 'ihspin_hn_0.asc', 'Optional')
			ihspin1 = importfile(datadir, 'ihspin_hn_1.asc', 'Optional')
		elif tlast_ih01>=tlast_ih34:
			ihspin0 = importfile(datadir, 'ihspin_hn_0.asc', 'Optional')
			ihspin1 = importfile(datadir, 'ihspin_hn_1.asc', 'Optional')
	    	elif tlast_ih34>tlast_ih01:
			ihspin0 = importfile(datadir, 'ihspin_hn_3.asc', 'Optional')
			ihspin1 = importfile(datadir, 'ihspin_hn_4.asc', 'Optional')
		
		 		    
	elif not spindata_from_hznfndr and spindata_from_sphrad:
		ihspin0 = importfile(datadir, 'ihspin_hn_3.asc', 'Optional')
		ihspin1 = importfile(datadir, 'ihspin_hn_4.asc', 'Optional')
	else:
		ihspin0 = ihspin0_bh1
		ihspin1 = ihspin1_bh2


	if os.path.isfile(ihspin0) and os.path.isfile(ihspin1):

		time1_spin = np.loadtxt(ihspin0, usecols = (0,), comments = '#')
		time2_spin = np.loadtxt(ihspin1, usecols = (0,), comments = '#')
		#minlen = min(len(time1_spin), len(time2_spin))
		time_spin = min(time1_spin[-1], time2_spin[-1])	#time1_spin[:minlen]
		
			
		if time_spin>=retarted_junktime:

			t1_idx = np.amin(np.where(time1_spin>=retarted_junktime))
			t2_idx = np.amin(np.where(time2_spin>=retarted_junktime))

			time_cutoff = time1_spin[t1_idx] 
			(s1[0],s1[1], s1[2]) = np.loadtxt(ihspin0, usecols= (1,2,3), comments = '#')[t1_idx]
			(s2[0],s2[1], s2[2]) = np.loadtxt(ihspin1, usecols= (1,2,3), comments = '#')[t2_idx]

			if verbose: print("*(Metadata) >> At time = {}, Spin1 = ({}, {}, {}) and Spin2 = ({}, {}, {}) using ihspins. \n".format(time_cutoff, s1[0], s1[1], s1[2], s2[0], s2[1], s2[2]))


		elif not(simtype=='precessing'):


			time_cutoff = time1_spin[-1]			

			(s1[0],s1[1], s1[2]) = np.loadtxt(ihspin0, usecols= (1,2,3), comments = '#')[-1]
			(s2[0],s2[1], s2[2]) = np.loadtxt(ihspin1, usecols= (1,2,3), comments = '#')[-1]
			
			if verbose: print("*(Metadata) >> Warning: Spin Information not available in after retarted junktime. Last available information will be taken from ihspin files\n")
			if verbose: print("*(Metadata) >> At time = {}, Spin1 = ({}, {}, {}) and Spin2 = ({}, {}, {}) using ihspins. \n".format(time_cutoff, s1[0], s1[1], s1[2], s2[0], s2[1], s2[2]))
		else:
			raise ValueError('Spin Information not available for %s simulation. Metadata cannot be generated \n'%simulation_name(dirpath))

		return [s1,s2,time_cutoff]
	else:	
		if not (simtype=='precessing'):
			if verbose: print("*(Metadata) >> Warning: Spin Information not available in IHSpin files. Spin Inf rmation will be taken from parfile\n")
		    	return [spin1, spin2, 0]
		else:
			raise ValueError('*(Metadata) >> Spin Function: Spin Information not found for precessing simulation. Metadata cannot be created \n')


def updatemass(dirpath, retarted_junktime, verbose):

	datadir = os.path.join(dirpath, "data")

	bh_diag0 = importfile(datadir,'BH_diagnostics.ah1.gp')
	bh_diag1 = importfile(datadir,'BH_diagnostics.ah2.gp')

	if os.path.exists(bh_diag1):
		time_bh1, irr_m1 = np.loadtxt(bh_diag0, usecols = (1,26), unpack=True, comments = '#')
		time_bh2, irr_m2 = np.loadtxt(bh_diag1, usecols = (1,26), unpack =True, comments = '#')
		#minlen = min(len(time_bh1), len(time_bh2))-1
		time_mass = min(time_bh1[-1], time_bh2[-1])	#time_bh1[:minlen]
		
		if (time_mass>=retarted_junktime):

			t1_idx = np.amin(np.where(time_bh1>=retarted_junktime))
			t2_idx = np.amin(np.where(time_bh2>=retarted_junktime))
			time_cutoff = time_bh1[t1_idx]
			m1 = irr_m1[t1_idx]	
			m2 = irr_m2[t1_idx]
		else:
			m1 = irr_m1[-1]
			m2 = irr_m2[-1]
			time_cutoff = min(time_bh1[-1], time_bh2[-1])
			if verbose:print("*(Metadata) >> Warning: Mass information not available after junk radiation. Mass will be provided at t = {} \n".format(time_cutoff))

		if verbose: print("*(Metadata) >> Mass Information provided from BH_Diagnostic files - m1 = {}, m2 = {} at t = {} \n".format(m1, m2, time_cutoff))
		return [m1,m2, time_cutoff]

	else: 		
		raise ValueError("*(Metadata) >> Mass Function - BH_diagnostics files not found. Metadata cannot be created \n")
		

def updatepos(dirpath, time_cutoff, verbose=False):

	datadir = os.path.join(dirpath, "data")	
	r1, r2, v1, v2 = np.zeros((4,3))	

	shifttracker0 = importfile(datadir, 'ShiftTracker0.asc', Importance='Required')
	shifttracker1 = importfile(datadir, 'ShiftTracker1.asc', Importance='Required')

	time_shtr= np.loadtxt(shifttracker0, usecols = (1,), comments = '#')
	temparr_shtr = np.where(time_shtr>=time_cutoff)
	cutoff_idx_shtr = np.amin(temparr_shtr)     
	
	( r1[0],r1[1],r1[2], v1[0], v1[1], v1[2]) = np.loadtxt(shifttracker0, usecols=( 2, 3, 4, 5, 6, 7), comments = '#')[cutoff_idx_shtr]
	( r2[0],r2[1],r2[2], v2[0], v2[1], v2[2]) = np.loadtxt(shifttracker1, usecols=( 2, 3, 4, 5, 6, 7), comments = '#')[cutoff_idx_shtr]
	
	if verbose: print("*(metadata) >> At time = {}, r1 = ({}, {}, {}) and r2 = ({}, {}, {}) \n".format(time_shtr[cutoff_idx_shtr], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]))
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

	data_shtr1 = np.loadtxt(shifttracker0, usecols = (1,2,3,4), comments = '#', unpack=True).T
	data_shtr2 = np.loadtxt(shifttracker1, usecols = (1,2,3,4), comments = '#', unpack=True).T

	if (len(data_shtr1)>len(data_shtr2)):
	    data_shtr1 = data_shtr1[:len(data_shtr2)]
	elif (len(data_shtr2)>len(data_shtr1)):
	    data_shtr2 = data_shtr2[:len(data_shtr1)]

	time_shtr1, r1x, r1y, r1z = data_shtr1.T
	time_shtr2, r2x, r2y, r2z = data_shtr2.T
	
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

	time_shtr = np.loadtxt(shifttracker0, usecols = (1), comments = '#')
	temparr_shtr = np.intersect1d(np.where(time_shtr>=time_cutoff),  np.where(time_shtr<=150))

	data_shtr1 = np.loadtxt(shifttracker0, usecols = (1,2,3,4,5,6,7), comments = '#', unpack=True).T
	data_shtr2 = np.loadtxt(shifttracker1, usecols = (1,2,3,4,5,6,7), comments = '#', unpack=True).T

	if (len(data_shtr1)>len(data_shtr2)):
	    data_shtr1 = data_shtr1[:len(data_shtr2)]
	elif (len(data_shtr2)>len(data_shtr1)):
	    data_shtr2 = data_shtr2[:len(data_shtr1)]

	time_shtr1, r1x, r1y, r1z, v1x, v1y, v1z = data_shtr1.T 
	time_shtr2, r2x, r2y, r2z, v2x, v2y, v2z = data_shtr2.T
	
	sepx = r1x - r2x
	sepy = r1y - r2y
	sepz = r1z - r2z

	sep = np.array((sepx, sepy, sepz))
	nhat = sep/mag(sep)

	dv = np.array((v1x-v2x, v1y-v2y, v1z-v2z))
	dnhat_dt = dv/mag(sep)
	
	omega_orbfreq = np.cross(nhat, dnhat_dt, axis=0)
	omega_orbfreq_mag = mag(omega_orbfreq) 

	cutoff_idx_shtr = np.amin(np.where(time_shtr1>=time_cutoff))
	omega_orbfreq_cutoff = omega_orbfreq_mag[cutoff_idx_shtr]
	return omega_orbfreq_cutoff
   

def determine_production_run(dirpath):

	produc_run = 1
	
	#Do not use runs with separation <8 as production runs:
	datadir = os.path.join(dirpath, "data")
	filename = dirpath.split("/")[-1]
	if filename[-3:]=='all':
    	    parfile = os.path.join(datadir, filename[:-4] + '.par')
	else:
    	    parfile = os.path.join(datadir, filename + '.par')
	
	
	if not os.path.isfile(parfile) :
	    raise ValueError(" Parfile missing in the directory. Please check again. \n")
	
	par_b = output_data(parfile, 'par_b')
	separation = 2.*par_b
	
	if separation<7: produc_run=0.
       
	#Do not use runs with resolution less than M/120

	regrid_boxes = False
	carpet_regrid2 = False
	datafile = file(parfile)
        datafile.seek(0)
        for line in datafile:
            if 'CarpetRegrid2' in line:
		carpet_regrid2 = True
                break
	    elif 'RegridBoxes' in line:
		regrid_boxes = True
		break
	datafile.close()	
 
	if carpet_regrid2:	
	    num_lvl1 = output_data(parfile,'num_levels_1')
	    num_lvl2 = output_data(parfile,'num_levels_2')
	    num_lvl_max = max(num_lvl1, num_lvl2)
	elif regrid_boxes:
	    datafile = file(parfile)
            datafile.seek(0)

    	    num_lvls_arr = []
            num_lvl_occ = 0

            for line in datafile:

                if 'refinement_levels' in line:
    		    num_lvl_occ +=1
    		    num_lvls = int((line.split('=')[1]).split()[0])
                    num_lvls_arr.append( num_lvls)
    		    if num_lvl_occ==2:
    		        break
    	
            datafile.close()	
	    num_lvl_max = np.amax(num_lvls_arr)
	

	dx = output_data(parfile, 'CoordBase::dx')
	dx_fine = dx/(2.**(num_lvl_max-1.))

	resolution = 1./dx_fine
	if resolution<120: produc_run=0
	
	print('*(Metadata) >> %s: Separation of this run = %g and resolution =M/%g, hence the production run status = %g'%(simulation_name(dirpath), separation, resolution, produc_run))

	return produc_run

def get_mass_from_stdout(dirpath):
        filename = dirpath+"/data/stdout"
        if not os.path.isfile(filename):
                filename = dirpath+"/data/stdout-1"
        if not os.path.isfile(filename):
                raise
        adm_mass_p = 0
        adm_mass_m = 0
        with open(filename) as f:
                line = f.readline()
                while line:
                        if('ADM mass from r= 1000000' in line):
                                adm_mass_p = float(f.readline().split('=')[-1])
                                adm_mass_m = float(f.readline().split('=')[-1])
                                break
                        line  = f.readline()
        if(adm_mass_p == 0 or adm_mass_m == 0):
                raise
        else:
                return adm_mass_p, adm_mass_m
