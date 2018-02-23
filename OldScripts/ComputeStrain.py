from nrutils import scsearch, gwylm, alert
import numpy as np
import os, shutil
import matplotlib.pyplot as plt
import multiprocessing as mp
#get_ipython().magic(u'matplotlib inline')

def write_strdata(data, outfile):
    strain = open(outfile,'w')
    hdr = "#Time \t hplus \t hcross \t amp \t phase \n"
    np.savetxt(strain, data, header=hdr, newline="\n", delimiter="\t" )
    strain.close()

    
def simulation_type(spin1, spin2):

    if (np.count_nonzero(spin1) ==0 and np.count_nonzero(spin2)==0):
        simtype = 'non-spinning'
    elif (np.count_nonzero(spin1[0:2])>0 or np.count_nonzero(spin2[0:2])>0) :
        simtype = 'precessing'
    else:
        simtype = 'aligned-spins'

    return simtype



def simulation_name(waveform):
	
	waveform_name = waveform.__dict__['simname']
	wf_junkrad = 'wf_junkrad.txt'
    
	wfdata = np.genfromtxt(wf_junkrad, dtype=None, comments='#', usecols=(0,1), skip_header=1, delimiter = '\t', names = ('GTID', 'simname'))
    	GTname, wfname = wfdata['GTID'], wfdata['simname']

	if np.array(np.where(wfname==waveform_name)).size==0:
		error('*(metadata) >> GT simulation name incorrectly listed in wf_junkrad.txt. Please check the file.')
 	elif np.array(np.where(wfname==waveform_name)).size>1:
		error('*(metadata) >> Multiple occurences of GT simulation - {} in wf_junkrad.txt - Please check the file.'.format(file_name))
	else: 
		idx = (np.where(wfname == waveform_name))[0][0]
		return GTname[idx]

#Path of output directories - 

def ComputeStrain(waveforms, outdir, lmax):


    #Search the waveform
    tempdir = "/localdata/bkhamesra3/LIGO_Waveforms/Waveform_Tests"
    wf = waveforms#scsearch(keyword=waveforms, institute='gt', verbose=True)
    wfdir = os.path.join(tempdir, wf[0].__dict__['simname'])
    numreldir = '/nethome/numrel/datafiles/Waveforms'
   
    #Create the necessary directories
    strain_dir = os.path.join(wfdir,"data/Strain")
    strfigdir = os.path.join(wfdir,"figures/strain_figures")
    datadir = os.path.join(wfdir,"data")
    figdir = os.path.join(wfdir,"figures")
    
    if not os.path.exists(strain_dir):
        os.makedirs(strain_dir)
    if not os.path.exists(strfigdir):
        os.makedirs(strfigdir)
    
    

    #Compute the metadata
    ylm_66 = os.path.join(os.path.join(numreldir, wf[0].__dict__['relative_simdir']), 'Ylm_WEYLSCAL4::Psi4r_l6_m6_r75.00.asc')
    ylm_88 = os.path.join(os.path.join(numreldir, wf[0].__dict__['relative_simdir']), 'Ylm_WEYLSCAL4::Psi4r_l8_m8_r75.00.asc')

    if not os.path.exists(ylm_88):
	lmax=6

    if not os.path.exists(ylm_66):
	lmax=2

    lm = []
    for l in range(2,lmax+1):
        for m in range(-l,l+1):
            lm.append((l,m))
    
    
    y = gwylm(wf[0], lm=lm, clean=True);                  
    #y.__dict__
    
    #Copy the relevant files
    simdir = y.__dict__['simdir']
    simname = y.__dict__['simname']
    

    def copyfile(filename):
    	filepath = os.path.join(simdir, filename)
	if os.path.exists(filepath):
	    shutil.copy(filepath, datadir)
	else:
	   print("%s >> Warning - %s file was not found"%(simname,filename))

    

    files = ["BH_diagnostics.ah1.gp","BH_diagnostics.ah2.gp", "BH_diagnostics.ah3.gp", (simname+".par"),"hn_mass_spin_0.asc", "hn_mass_spin_1.asc", "ihspin_hn_0.asc", "ihspin_hn_1.asc", "ihspin_hn_2.asc", 
    "ihspin_hn_3.asc", "ihspin_hn_4.asc", "ShiftTracker0.asc", "ShiftTracker1.asc"]
    
    for f in files:
	copyfile(f)

    s1 = y.__dict__["S1"]
    s2 = y.__dict__["S1"]

    simtype=simulation_type(s1, s2)

    #Update the output directory based on spins
    if simtype=="non-spinning":
        outdir = os.path.join(outdir, "NonSpinning")
    elif simtype=="aligned-spins":
        outdir = os.path.join(outdir, "AlignedSpin")
    else:
        outdir = os.path.join(outdir, "Precessing")
    
    
    #Compute the strain plots
    for l,m in lm:
        [time, hp, hx] = y.lm[l,m]['strain'].wfarr.T
        data = np.column_stack((time, hp, hx))
        write_strdata(data, os.path.join(strain_dir,"Strain_l%d_m%d.asc"%(l,m)))
        
        plt.plot(time, hp, 'b', label='h+')
        plt.plot(time, hx, color='darkorange', label='hx')
        plt.ylabel('Strain')
        plt.xlabel('Time')
        plt.title("Strain - (%d,%d) Mode"%(l,m))
        plt.legend()
        plt.savefig(os.path.join(strfigdir, "Strain_%d%d.png"%(l,m)), dpi=500)
        plt.close()

    print("Strain plots created for %s waveform"%simname)
    shutil.move(wfdir, outdir)
    



lmax = 8

outdir = '/nethome/numrel/datafiles/Finalized_Waveforms/Waveform_files/New_BK'
if not os.path.exists(outdir):
    os.makedirs(outdir)

waveform = scsearch(verbose=True,institute='gt')
pool = mp.Pool(processes=4)	
for wf in waveform:
    ComputeStrain(waveform, outdir, lmax)
