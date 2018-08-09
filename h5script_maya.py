'''
File to create hdf5 file for LAL injection from NS NR data.
		
The hdf5 files are to be used in the "LIGO NR-Injection infrastructure":https://arxiv.org/pdf/1703.01076.pdf

Notes - Keep the junkradiation file in the Waveform directory or modify its path

'''

from numpy import loadtxt,sqrt,angle,unwrap,diff,pi
import matplotlib
matplotlib.use('Agg')

from matplotlib.pyplot import *
import glob, os
import romspline
from bundlers import nr2h5
from metadata import metadata
from metadata_functions import output_data, error, simulation_name, find_omega22_st
from init_data import initial_data
import shutil as sh

def crop_data(dirpath, time, hp, hx):
 
    #Open the text file which has the data about the time duration of junk radiation    
    file_name = dirpath.split('/')[-1]
    wf_junkrad = 'wf_junkrad.txt'
    
    wfdata = np.genfromtxt(wf_junkrad, dtype=None, comments='#', usecols=(1,4), skip_header=1, delimiter = '\t', names = ('simname', 'time'))
    
    wfname, jkrad_time = wfdata['simname'], wfdata['time']
    
    
    #Find the index where to crop off the data
    
    if np.array(np.where(wfname==file_name)).size==0:
	msg = 'GT simulation not in wf_junkrad.txt. Please check the file.'
	raise ValueError(msg)
    elif np.array(np.where(wfname==file_name)).size>1:
	msg = 'Multiple occurences of GT simulation - {} in wf_junkrad.txt - Please check the file.'.format(file_name)
	raise NameError(msg)
    else: 
	idx = (np.where(wfname == file_name))[0][0]
	

    if np.array(np.where(time==jkrad_time[idx])).size==0:
	temp_arr = np.where(time>jkrad_time[idx])
	cutoff_idx = np.amin(temp_arr)
    else:
	cutoff_idx = (np.where(time==jkrad_time[idx]))[0][0]  
    

    #Crop the data to remove junk radiation
    time = time[cutoff_idx:]
    hp = hp[cutoff_idx:]
    hx = hx[cutoff_idx:]
    
    return (time, hp, hx, cutoff_idx, jkrad_time[idx])



def create_plot(varx, vary, label_x, label_y, figpath, figname):

	plot(varx, vary )
	xlabel(label_x)
	ylabel(label_y)
	savefig(figpath+'/'+figname)
	close()
  

def create_single_h5 (dirpath, movepath_h5, movepath_wf, verbose=False):			 #dirpath - path to directory which contains the strain data of the waveform
			

    file_name = dirpath.split('/')[-1]
    print("*(Create_Single_h5) >> Constructing h5 waveform for simulation - {} \n ".format(file_name)) 
  
    #Add path of directories to read data and save figures- 
    figdir = os.path.join(dirpath, "figures")
    datadir = os.path.join(dirpath, "data")
    strdir = os.path.join(datadir, "Strain")

    #Check if the h5 file already exists
    simname = simulation_name(dirpath)
    h5check_path = glob.glob(os.path.join(movepath_h5, "*/%s.h5"%simname))

    if len(h5check_path)>0:
        if os.path.exists(h5check_path[0]):
            raise NameError("*(Create_Single_h5) >> h5file %s already exists. Please delete the old file before constructing the new one. Waveform getting moved to Failed Directory"%simname)


    #Create a dictionary to store all the modes 
    nr_strain_data = {}

    lmax=ltemp = 0


    # Safety Check
    if not os.path.isfile(strdir+'/Strain_l2_m2.txt'): 
	error('*(Create_Single_h5) >> Strain Files Missing \n')
    #else:
 	#t22, hp_22, hx_22 = np.loadtxt(strdir + "/Strain_l2_m2.txt", comments="#", usecols=(0,1,2), unpack=True)
	#amp_22 = abs(hp_22 - 1j*hx_22)
	#t_maxamp = t22[amp_22 == np.amax(amp_22)]
   	#print("*(Create_Single_h5) >> Time at max amplitude = %f \n"%t_maxamp)


    #The for loop ensures all the modes are included in this dictionary
    for strfile in sorted(glob.glob(strdir+'/Strain*.txt')):
    
	# NOTE: here we determine l and m from the filename. This assumes that each
	# strain file has data for only one l and one m.
	l = int((os.path.basename(strfile).split('_')[1])[1])
	m = int((os.path.basename(strfile).split('_')[-1]).split('.')[0][1:])  

	# Find Lmax for metadata
    	ltemp = l
	if ltemp>lmax:
	    lmax = ltemp
   
	# Load the data
	data_array = loadtxt(strfile, comments='#')#, skiprows=2)

    	# Unpack the data array (get time, plus and cross )
	time_uncrop, hp_uncrop, hx_uncrop  = data_array[:,0], data_array[:,1], data_array[:,2]

    	# Crop the data to remove the junk radiation
	time_crop, hp, hx, crop_idx, junkrad_time = crop_data(dirpath, time_uncrop, hp_uncrop, hx_uncrop)

	time = time_crop 


	# Calculate the waveform amplitude and phase  (Unwrapping phase  - Request from NR-injection Infrastructure
	A = abs(hp -1j*hx)
	phase = -1.0*unwrap(angle(hp - 1j*hx))  #Check this with Juan

    	# NOTE: The NR-Injection infrastructure requests that phase decreases if m is positive. Here we multiply by -1 to enforce this.
	phase *= -1

	# Set all m=0 mode to 0:
	# This is just temporary measure to fix the issue with m=0 modes. Better method would be to find the right tuning frequency by comparison with Quasi normal modes
	if m==0: #l==2 and m==0:
	    A *= 0.
	    phase *=0.	


    	# Create a dictionary that contains the NR data. This will be used as an input to the function that actually makes the h5 file.
        nr_strain_data[ (l,m) ] = { 'amp':A, 'phase':phase, 't':time }
   

    # Compute the index of maximum of strain amplitude - all modes included

    amp_sq = np.zeros(len(nr_strain_data[(2,2)]['t']))

    for key in nr_strain_data:
	 amp_sq = amp_sq + nr_strain_data[key]['amp']**2 
	
    amp_max_idx = np.where(amp_sq==np.amax(amp_sq))
    amp22_max_idx = np.where(nr_strain_data[(2,2)]['amp']==np.amax(nr_strain_data[(2,2)]['amp']))

	
    # Check for difference between peak of amplitude-all modes,  and peak of 2,2 mode amplitude
    
    diff = (nr_strain_data[key]['t'][amp_max_idx] - nr_strain_data[key]['t'][amp22_max_idx])
    if diff>5:
	raise NameError('*(Create_single_h5) >>  Unexpected difference of %gM observed between time of peak of 22 mode and time of peak of strain amplitude (all modes). Some mode is behaving incorrectly'%diff)
    elif verbose:
	print("*(Create_single_h5) >>  Time difference between peak of strain amplitude and (2,2) mode amplitude is %gM \n"%(diff))


    # Mean Center the data about max amplitude

    for key in nr_strain_data:

	time = nr_strain_data[key]['t'] - nr_strain_data[key]['t'][amp_max_idx]
	nr_strain_data[key]['t'] = time	

	# Create plots of amplitude and phase for mean centered time 
	if key==(2,2):
	    create_plot(time, nr_strain_data[key]['amp'], 'Time', 'h+', figdir,  'Amp_22_plot.png')
            create_plot(time, nr_strain_data[key]['phase'], 'Time', 'Phase', figdir,'Phase_22_plot.png')
    	
	amp_lm = nr_strain_data[key]['amp']	#t22 - t22[abs(h22)==np.amax(abs(h22))]
   	phase_lm = nr_strain_data[key]['phase']				#unwrap(angle(hp_22 - 1j*hx_22))

	
	#Testing Purpose Only - Save the file and compare with spline interpolant
	data = np.column_stack((time, amp_lm, phase_lm))
        hdr = "Time \t Amplitude \t Phase \n"
	
	timeshift_dir = os.path.join(strdir, "Strain_TimeShifted")
	if not os.path.exists(timeshift_dir): os.makedirs(timeshift_dir)
        l,m = key[0], key[1]
	strfile_lm = os.path.join(timeshift_dir,"Strain_{}{}_timeshift.asc".format(l,m))
        np.savetxt(strfile_lm, data, header=hdr, delimiter='\t')


  
    # Compute the gravitational wave frequency

    time_22 = nr_strain_data[(2,2)]['t'] 
    amp_22 = nr_strain_data[(2,2)]['amp']	#t22 - t22[abs(h22)==np.amax(abs(h22))]
    phase_22 = nr_strain_data[(2,2)]['phase']				#unwrap(angle(hp_22 - 1j*hx_22))

    omega_22 = -1.*np.gradient(phase_22)/np.gradient(time_22)#(phase_uncrop[1:] - phase_uncrop[:-1])/(time_uncrop[1:]-time_uncrop[:-1])

    create_plot(time_22, omega_22, 'Time', 'MOmega', figdir, 'Omega_22_plot.png')
	
   
    # Compute Momega(2,2)
    time_200_idx = np.amin(np.where(time_22>=time_22[0]+200))
    Mf_initial = omega_22[0]#np.amin(omega_22[0:time_200_idx])
		 
   
    #Check with frequency of ShiftTracker
    omega_st = find_omega22_st(dirpath, junkrad_time-75)
    if verbose:
	print('*(Create_single_h5) >> Gravitational Wave frequency from strain = %g/M, and from orbital frequency = %g/M \n'%(Mf_initial, 2.*omega_st[0]))

    #if verbose:
    print("*(Create_single_h5) >>  Frequency of (2,2) mode after junk radiation -  Momega = {} \n".format(Mf_initial))

    # Below we estimate the starting frquency of the waveform scaled at 1Msolar
    G = 6.67428e-11 # m^3/(kg s^2)
    mass_sun = 1.98892e30 # kg
    c = 2.99792458e8 # m/s
    mass_sec = G * mass_sun / ( c*c*c )
    simulation_mass_sec = ( 1.0 ) * mass_sec # NOTE that the 1.0 here is for 1solar mass
    physical_freq_initial_at_1Msol = (Mf_initial/(2.*np.pi)) / simulation_mass_sec
   
    if verbose: 
	print('*(Create_single_h5) >> Gravitational Wave frequency from strain = %gHz, and from orbital frequency = %gHz \n'%(physical_freq_initial_at_1Msol, (2.*omega_st[0]/(2.*np.pi)/simulation_mass_sec)))
   
    # Compute initial data
    nr_initdata = initial_data(dirpath) 

    # Compute the nr_metadata
    nr_meta_data = metadata(dirpath, junkrad_time)
    nr_meta_data['Lmax'] = lmax
    nr_meta_data['Omega'] = Mf_initial
    nr_meta_data['f_lower_at_1MSUN'] = physical_freq_initial_at_1Msol

    # Create the h5 file and move to corresponding directory
    outpath = movepath_h5+('/%s.h5'%nr_meta_data['name'])

    print('*(Create_single_h5) >> Creating HDF5 files:{}'.format(outpath))
    nr2h5( nr_strain_data, nr_meta_data, nr_initdata, output_path=outpath, wfdir_path=dirpath, verbose=True, testsuite=False )

    if nr_meta_data['simulation-type']=='non-spinning':
	movepath_h5 = os.path.join(movepath_h5, 'NonSpinning')
	movepath_wf = os.path.join(movepath_wf, 'NonSpinning')
    elif nr_meta_data['simulation-type'] == 'aligned-spins':
	movepath_h5 = os.path.join(movepath_h5, 'AlignedSpin')
	movepath_wf = os.path.join(movepath_wf, 'AlignedSpin')
    elif nr_meta_data['simulation-type'] == 'precessing':
	movepath_h5 = os.path.join(movepath_h5, 'Precessing')
	movepath_wf = os.path.join(movepath_wf, 'Precessing')

    sh.move(dirpath, movepath_wf)
    sh.move(outpath, movepath_h5)


#wf_direc = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/Precessing/Lq_D6.2_q2.50_a0.6_th015_m140"
h5output_path = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/H5Files" 
wfoutput_path = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Completed"
failed_path = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Failed/NonSpinning"
wf_direc = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/NonSpinning"


#create_single_h5(wf_direc, h5output_path, wfoutput_path, verbose=True)

for direc in os.listdir(wf_direc):
 	direc_path = os.path.join(wf_direc,direc)
	if os.path.isdir(direc_path):
		print("(h5script_Maya)* >> Starting the python script to create h5 files for waveform - {} \n".format(direc))
		try:
  			create_single_h5(direc_path, h5output_path, wfoutput_path, verbose=True)
		except (ValueError, NameError) as error :	#Remove except condition if testing Failed waveforms
			print("An Error occured. The file is being moved to Failed Directory")
			sh.move(direc_path, failed_path )
		except (sh.Error) as error :	#Remove except condition if testing Failed waveforms
			print("An Error occured. The file is being moved to Failed Directory")
			sh.move(direc_path, failed_path )
