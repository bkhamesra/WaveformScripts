'''
File to create hdf5 file for LAL injection from NS NR data.
		
The hdf5 files are to be used in the "LIGO NR-Injection infrastructure":https://dcc.ligo.org/DocDB/0123/T1500606/002/NRInjectionInfrastructure.pdf

Notes - Keep the junkradiation file in the Waveform directory or modify its path

'''

from numpy import loadtxt,sqrt,angle,unwrap,diff,pi
import matplotlib
matplotlib.use('Agg')

from matplotlib.pyplot import *
import glob, os, shutil
import romspline
from bundlers import nr2h5
from metadata import metadata
from metadata_functions import output_data, error
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
	msg = 'GT simulation name incorrectly listed in wf_junkrad.txt. Please check the file.'
	error(msg)
    elif np.array(np.where(wfname==file_name)).size>1:
	msg = 'Multiple occurences of GT simulation - {} in wf_junkrad.txt - Please check the file.'.format(file_name)
	error(msg)
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
			

    #Change to this directory and find the file name
    #os.chdir(dirpath)
    file_name = dirpath.split('/')[-1]
    print(" \n Constructing h5 waveform for simulation - {} \n ".format(file_name)) 
  
    #Construct directories - 
    figdir = os.path.join(dirpath, "figures")
    datadir = os.path.join(dirpath, "data")
    strdir = os.path.join(datadir, "Strain")


    #Create a dictionary to store all the modes 
    nr_strain_data = {}

    lmax=ltemp = 0


    # Safety Check
    if not os.path.isfile(strdir+'/Strain_l2_m2.txt'): 
	error('Strain Files Missing \n')
    else:
 	t22, hp_22, hx_22 = np.loadtxt(strdir+"/Strain_l2_m2.txt", comments="#", usecols=(0,1,2), unpack=True)
	amp_22 = abs(hp_22+1j*hx_22)
	t_maxamp = t22[amp_22==np.amax(amp_22)]
   	print("h5_maya >> Time at max amplitude = %f"%t_maxamp)

    #The for loop ensures all the modes are included in this dictionary
    for strfile in glob.glob(strdir+'/Strain*.txt'):
    
	# NOTE: here we determine l and m from the filename. This assumes that each
	# strain file has data for only one l and one m.
	l =int((os.path.basename(strfile).split('_')[1])[1])
	m =int((os.path.basename(strfile).split('_')[-1]).split('.')[0][1:])  

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

	#Mean Center the data about max amplitude
	time = time_crop - t_maxamp


	# Calculate the waveform amplitude and phase  '''Comment - Why unwrapping'''
	A = abs(hp +1j*hx)
	phase = -unwrap(angle(hp + 1j*hx))  #Check this with Juan

    	# NOTE: The NR-Injection infrastructure requests that phase decreases if m is positive. Here we multiply by -1 to enforce this.
	phase *= -1

	# Set all m=0 mode to 0	
	if m==0: #l==2 and m==0:
		A *= 0.
		phase *=0.	

    	# Create plots of strain, amplitude, phase and omega for (2,2) mode as a check. 	
	if l==2 and m==2:
		
		phase_uncrop = unwrap(angle(hp_uncrop+1j*hx_uncrop))
		omega = (phase_uncrop[1:] - phase_uncrop[:-1])/(time_uncrop[1:]-time_uncrop[:-1])

		create_plot(time, hp, 'Time', 'h+', figdir,  'hplus_22_plot.png')
		create_plot(time, hx, 'Time', 'hx', figdir, 'hx_22_plot.png')

	      	create_plot(time, A, 'Time', 'Amplitude', figdir, 'Amp_22_plot.png' )

		create_plot(time, phase, 'Time', 'Phase', figdir,'Phase_22_plot.png')

		create_plot(time_uncrop[1:], omega, 'Time', 'MOmega', figdir, 'Omega_22_plot.png')
	
	# Compute Momega(2,2)
	if l==2 and m==2:
	
		Mf_initial = omega[crop_idx]
		 
		if verbose:
			print("\n *(Create_single_h5)>>  Frequency of (2,2) mode -  Momega = {} \n".format(Mf_initial))

    	# Create a dictionary that contains the NR data. This will be used as an input to the function that actually makes the h5 file.
        nr_strain_data[ (l,m) ] = { 'amp':A, 'phase':phase, 't':time }
    
	
	
    # Below we estimate the starting frquency of the waveform scaled at 1Msolar
    G = 6.67428e-11 # m^3/(kg s^2)
    mass_sun = 1.98892e30 # kg
    c = 2.99792458e8 # m/s
    mass_sec = G * mass_sun / ( c*c*c )
    simulation_mass_sec = ( 1.0 ) * mass_sec # NOTE that the 1.0 here is for 1solar mass
    physical_freq_initial_at_1Msol = (Mf_initial/(2.*np.pi)) / simulation_mass_sec
 
   
    # Compute initial data
    nr_initdata = initial_data(dirpath) 

    # Compute the nr_metadata
    nr_meta_data = metadata(dirpath, junkrad_time)
    nr_meta_data['Lmax'] = lmax
    nr_meta_data['Omega'] = Mf_initial
    nr_meta_data['f_lower_at_1MSUN'] = physical_freq_initial_at_1Msol

    outpath = movepath_h5+('/%s.h5'%nr_meta_data['name'])
    print(' \n Creating HDF5 files:{}'.format(outpath))
    nr2h5( nr_strain_data, nr_meta_data, nr_initdata, output_path=outpath, verbose=True )

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


#wf_direc = "/nethome/numrel/datafiles/Finalized_Waveforms/Waveform_files/Remaining/Precessing/fr_b5_a0.6_oth.315_M77"
h5output_path = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms" 
wfoutput_path = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Completed"
failed_path = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Failed/NonSpinning"
wf_direc = "/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/NonSpinning"


#create_single_h5(wf_direc, h5output_path, wfoutput_path, verbose=True)

for direc in os.listdir(wf_direc):
 	direc_path = os.path.join(wf_direc,direc)
	print direc_path
	if os.path.isdir(direc_path):
		print("Starting the python script to create h5 files for waveform - {}".format(direc))
		try:
  		    create_single_h5(direc_path, h5output_path, wfoutput_path, verbose=True)
		except (ValueError, shutil.Error):	#Remove except condition if testing Failed waveforms
			print("An Error occured. The file is being moved to Failed Directory")
			sh.move(direc_path, failed_path )
