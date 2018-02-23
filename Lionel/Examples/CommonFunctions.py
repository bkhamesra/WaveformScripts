
# coding: utf-8

# In[2]:

from nrutils import scsearch, gwylm, alert
import numpy as np
import os, shutil, sys, h5py

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


import multiprocessing as mp


# In[4]:

#Record initial stdout to switch on and off the verbose output of functions
old_stdout = sys.stdout
py3 = sys.version_info[0] > 2 #creates boolean value for test that Python major version > 2

os.chdir('/localdata/bkhamesra3/LIGO_Waveforms/Waveform_Comparison/Mode_20/')


# In[5]:

# Disable Print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore Print
def enablePrint():
    sys.stdout = old_stdout


# # Psi4 plots 

# In[7]:

def amp(re,im):
    return np.sqrt(re**2. + im**2.)

def gatech_psi4plots(y):
    
    psitime_20, re_20, im_20 = y.lm[2,0]['psi4'].wfarr.T
    psitime_21, re_21, im_21 = y.lm[2,1]['psi4'].wfarr.T
    psitime_22, re_22, im_22 = y.lm[2,2]['psi4'].wfarr.T
    
    tmax = psitime_22[amp(re_22, im_22)==np.amax(amp(re_22, im_22))]
    
    ax = []
    plt.figure(figsize=(10,8))
    
    ax.append(plt.subplot2grid((2,2), (0,0), colspan=2))
    ax[0].plot(psitime_20, amp(re_20, im_20),'k--', lw=1, label="Amplitude")
    ax[0].plot(psitime_20, re_20, 'b', lw=1)
    
    ax.append(plt.subplot2grid((2,2), (1,0)))
    ax[1].plot(psitime_21, amp(re_21, im_21),'k--', lw=1)
    ax[1].plot(psitime_21, re_21,'r')
    
    ax.append(plt.subplot2grid((2,2), (1,1)))
    ax[2].plot(psitime_22, amp(re_22, im_22),'k--', lw=2)
    ax[2].plot(psitime_22, re_22, 'g')
    #ax[2].plot(t20, im20, 'g')
    
    for i in range(3):
        ax[i].set_ylabel('Psi4 (2,%d)'%i)
        ax[i].set_xlabel('Time')
        ax[i].set_xlim(tmax-200, tmax+200)
        
    plt.tight_layout()
    plt.show()
    plt.close()


def gatech_strainplots(y):
    
    time_20, hp_20, hx_20 = y.lm[2,0]['strain'].wfarr.T
    time_21, hp_21, hx_21 = y.lm[2,1]['strain'].wfarr.T
    time_22, hp_22, hx_22 = y.lm[2,2]['strain'].wfarr.T
    
    tmax = time_22[amp(hp_22, hx_22)==np.amax(amp(hp_22, hx_22))]
    
    ax = []
    plt.figure(figsize=(10,8))
    ax.append(plt.subplot2grid((2,2), (0,0), colspan=2))
    ax[0].plot(time_20, amp(hp_20, hx_20),'k--', lw=1, label="Amplitude")
    ax[0].plot(time_20, hp_20, 'b', lw=1)
    
    ax.append(plt.subplot2grid((2,2), (1,0)))
    ax[1].plot(time_21, amp(hp_21, hx_21),'k--', lw=1)
    ax[1].plot(time_21, hp_21,'r')
    
    ax.append(plt.subplot2grid((2,2), (1,1)))
    ax[2].plot(time_22, amp(hp_22, hx_22),'k--', lw=2)
    ax[2].plot(time_22, hp_22, 'g')
    #ax[2].plot(t20, im20, 'g')
    
    for i in range(3):
        ax[i].set_ylabel('Strain (2,%d)'%i)
        ax[i].set_xlabel('Time')
        ax[i].set_xlim(tmax-200, tmax+200)
        
    plt.tight_layout()
    plt.show()
    plt.close()



# In[8]:

#Strain-Psi4 Comparison
def gatech_strain_psi4_compare(y):
   
    time_20, hp_20, hx_20 = y.lm[2,0]['strain'].wfarr.T
    time_21, hp_21, hx_21 = y.lm[2,1]['strain'].wfarr.T
    time_22, hp_22, hx_22 = y.lm[2,2]['strain'].wfarr.T
      
    psitime_20, re_20, im_20 = y.lm[2,0]['psi4'].wfarr.T
    psitime_21, re_21, im_21 = y.lm[2,1]['psi4'].wfarr.T
    psitime_22, re_22, im_22 = y.lm[2,2]['psi4'].wfarr.T
    
       
    tmax = time_22[amp(hp_22,hx_22)==np.amax(amp(hp_22, hx_22))]
    ax = []
    time = [time_20, psitime_20, time_21, psitime_21, time_22, psitime_22]
    re = [hp_20, re_20, hp_21, re_21, hp_22, re_22]
    im = [hx_20, im_20, hx_21, im_21, hx_22, im_22]
    
    plt.figure(figsize=(15,10))
    
    for i in range(6):
        ax.append( plt.subplot2grid((3,2), (i/2,i%2)))
        ax[i].plot(time[i], amp(re[i], im[i]),'k--', lw=1)
        ax[i].plot(time[i], re[i], 'b', lw=1)
        ax[i].set_xlabel("Time")
        ax[i].set_xlim(tmax-200, tmax+200)
        
        if i%2==0: ax[i].set_ylabel('Strain (2,%d)'%(i/2))
        else: ax[i].set_ylabel('Psi4 (2,%d)'%(i/2))
    
    
    
    plt.tight_layout()
    plt.show()
    plt.close()


# In[9]:

#SXS - Psi4
def sxs_psi4plots(wf_sxs):
    
    data_22 = wf_sxs['Y_l2_m2.dat']
    data_21 = wf_sxs['Y_l2_m1.dat']
    data_20 = wf_sxs['Y_l2_m0.dat']

    t22, re22, im22 = np.transpose(data_22)
    t21, re21, im21 = np.transpose(data_21)
    t20, re20, im20 = np.transpose(data_20)
    
    tmax = t22[amp(re22,im22)==np.amax(amp(re22,im22))]
    ax = []
    
    plt.figure(figsize=(10,8))
    ax.append( plt.subplot2grid((2,2), (0,0), colspan=2))
    ax[0].plot(t20, amp(re20, im20),'k--')
    ax[0].plot(t20, re20, 'b', lw=1)
    
    ax.append(plt.subplot2grid((2,2), (1,0)))
    ax[1].plot(t21, amp(re21, im21),'k--')
    ax[1].plot(t21, re21,'r', lw=1)
    
    ax.append(plt.subplot2grid((2,2), (1,1)))
    ax[2].plot(t22, amp(re22, im22),'k--')
    ax[2].plot(t22, re22, 'b', lw=1)
    #ax3.plot(t20, im20, 'g')
    
    for i in range(3):
        ax[i].set_ylabel('Psi4 (2,%d)'%i)
        ax[i].set_xlabel('Time')
        ax[i].set_xlim(tmax-200, tmax+200)
    
    plt.tight_layout()
    plt.show()
    plt.close()


#SXS Strain Plots
def sxs_strainplots(wf_sxs):
    
    data_20 = np.loadtxt(wf_sxs['20'])
    data_21 = np.loadtxt(wf_sxs['21'])
    data_22 = np.loadtxt(wf_sxs['22'])

    t22, re22, im22 = data_22.T
    t21, re21, im21 = data_21.T
    t20, re20, im20 = data_20.T
    
    tmax = t22[amp(re22,im22)==np.amax(amp(re22,im22))]
    ax = []
    
    plt.figure(figsize=(10,8))
    ax.append( plt.subplot2grid((2,2), (0,0), colspan=2))
    ax[0].plot(t20, amp(re20, im20),'k--', lw=2)
    ax[0].plot(t20, re20, 'b')
    
    ax.append(plt.subplot2grid((2,2), (1,0)))
    ax[1].plot(t21, amp(re21, im21),'k--', lw=2)
    ax[1].plot(t21, re21,'r')
    
    ax.append(plt.subplot2grid((2,2), (1,1)))
    ax[2].plot(t22, amp(re22, im22),'k--', lw=2)
    ax[2].plot(t22, re22, 'g')
    #ax3.plot(t20, im20, 'g')
    
    for i in range(3):
        ax[i].set_ylabel('Strain (2,%d)'%i)
        ax[i].set_xlabel('Time')
        ax[i].set_xlim(tmax-200, tmax+200)
    
    plt.tight_layout()
    plt.show()
    plt.close()


# In[10]:

def sxs_strain_psi4_compare(str_sxs, wf_sxs, wftype):
   
    time_20, hp_20, hx_20 = np.loadtxt(str_sxs['20']).T
    time_21, hp_21, hx_21 = np.loadtxt(str_sxs['21']).T
    time_22, hp_22, hx_22 = np.loadtxt(str_sxs['22']).T
      
    psitime_20, re_20, im_20 = np.transpose(wf_sxs['Y_l2_m0.dat'])
    psitime_21, re_21, im_21 = np.transpose(wf_sxs['Y_l2_m1.dat'])
    psitime_22, re_22, im_22 = np.transpose(wf_sxs['Y_l2_m2.dat'])
    
       
    tmax = time_22[amp(hp_22,hx_22)==np.amax(amp(hp_22, hx_22))]
    
    psitmax = psitime_22[amp(re_22,im_22)==np.amax(amp(re_22, im_22))]
    psitime_20 = psitime_20 - psitmax
    psitime_21 = psitime_21 - psitmax
    psitime_22 = psitime_22 - psitmax
    
    ax = []
    time = [time_20, psitime_20, time_21, psitime_21, time_22, psitime_22]
    re = [hp_20, re_20, hp_21, re_21, hp_22, re_22]
    im = [hx_20, im_20, hx_21, im_21, hx_22, im_22]
    
    plt.figure(figsize=(15,10))
    
    for i in range(6):
        ax.append( plt.subplot2grid((3,2), (i/2,i%2)))
        ax[i].plot(time[i], amp(re[i], im[i]),'k--', lw=1)
        ax[i].plot(time[i], re[i], 'b', lw=1)
        ax[i].set_xlabel("Time")
        ax[i].set_xlim(tmax-200, tmax+200)
        
        if i%2==0: ax[i].set_ylabel('Strain (2,%d)'%(i/2))
        else: ax[i].set_ylabel('Psi4 (2,%d)'%(i/2))
    
    
    # plt.title('SXS - %s waveform'%wftype)
    plt.tight_layout()
    plt.show()
    plt.close()


# In[ ]:


def gatech_sxs_psi4_compare(wf_sxs, y, wftype):
   
    psitime_20, re_20, im_20 = np.transpose(wf_sxs['Y_l2_m0.dat'])
    psitime_21, re_21, im_21 = np.transpose(wf_sxs['Y_l2_m1.dat'])
    psitime_22, re_22, im_22 = np.transpose(wf_sxs['Y_l2_m2.dat'])
    
       
    psitmax = psitime_22[amp(re_22,im_22)==np.amax(amp(re_22, im_22))]
    psitime_20 = psitime_20 - psitmax
    psitime_21 = psitime_21 - psitmax
    psitime_22 = psitime_22 - psitmax
    
    ax = []
    time = [psitime_20, psitime_22]
    re = [re_20, re_22]
    im = [im_20, im_22]
    
    plt.figure(figsize=(15,10))
    
    for i in range(2):
        ax.append( plt.subplot2grid((2,2), (i%2,0)))
        ax[i].plot(time[i], amp(re[i], im[i]),'k--', lw=1)
        ax[i].plot(time[i], re[i], 'b', lw=1)
        ax[i].set_xlabel("Time")
        ax[i].set_xlim(-200, 200)
        ax[i].set_title("SXS")
        ax[i].set_ylabel('Psi4 (2,%d)'%(i*4/2))
    
    
    psitime_20, re_20, im_20 = np.transpose(wf_sxs['Y_l2_m0.dat'])
    psitime_21, re_21, im_21 = np.transpose(wf_sxs['Y_l2_m1.dat'])
    psitime_22, re_22, im_22 = np.transpose(wf_sxs['Y_l2_m2.dat'])


    psitmax = psitime_22[amp(re_22,im_22)==np.amax(amp(re_22, im_22))]
    psitime_20 = psitime_20 - psitmax
    psitime_21 = psitime_21 - psitmax
    psitime_22 = psitime_22 - psitmax
    
    ax = []
    time = [psitime_20,  psitime_22]
    re = [re_20, re_22]
    im = [im_20, im_22]
    
    for i in range(2):
        ax.append( plt.subplot2grid((2,2), (i%2,1)))
        ax[i].plot(time[i], amp(re[i], im[i]),'k--', lw=1)
        ax[i].plot(time[i], re[i], 'darkorange', lw=1)
        ax[i].set_xlabel("Time")
        ax[i].set_xlim(-200, 200)
        ax[i].set_title("GT")
        ax[i].set_ylabel('Psi4 (2,%d)'%(4*i/2))
    
    plt.tight_layout()
    plt.show()
    plt.close()


    
def gatech_sxs_strain_compare(str_sxs, y, wftype, w=-1):
   

    lm = []
    for l in range(2,3):
        for m in range(0,l+1):
            lm.append((l,m))
    
    sim = scsearch(keyword=y.simname, verbose=True,institute='gt')[0]
    
    y_20 = gwylm(sim, [2,0], clean=True, w22=w)
    y = gwylm(sim, lm, clean=True)
    
    time_20, hp_20, hx_20 = np.loadtxt(str_sxs['20']).T
    time_21, hp_21, hx_21 = np.loadtxt(str_sxs['21']).T
    time_22, hp_22, hx_22 = np.loadtxt(str_sxs['22']).T
      
    tmax = time_22[amp(hp_22,hx_22)==np.amax(amp(hp_22, hx_22))]
    
    ax = []
    time = [time_20, time_22]
    re = [hp_20, hp_22]
    im = [hx_20, hx_22]
    
    plt.figure(figsize=(15,10))
    
    for i in range(2):
        ax.append( plt.subplot2grid((2,2), (i%2,0)))
        ax[i].plot(time[i], amp(re[i], im[i]),'k--', lw=1)
        ax[i].plot(time[i], re[i], 'b', lw=1)
        ax[i].set_xlabel("Time")
        ax[i].set_xlim(tmax-200, tmax+200)
        ax[i].set_ylabel('Strain (2,%d)'%(4*i/2))
        ax[i].set_title("SXS")
    
    
    time_20, hp_20, hx_20 = y_20.lm[2,0]['strain'].wfarr.T
    time_21, hp_21, hx_21 = y.lm[2,1]['strain'].wfarr.T
    time_22, hp_22, hx_22 = y.lm[2,2]['strain'].wfarr.T
    
    tmax = time_22[amp(hp_22,hx_22)==np.amax(amp(hp_22, hx_22))]
    ax = []
    time = [time_20-tmax, time_22-tmax]
    re = [hp_20, hp_22]
    im = [hx_20, hx_22]
     
    for i in range(2):
        ax.append( plt.subplot2grid((2,2), (i%2,1)))
        ax[i].plot(time[i], amp(re[i], im[i]),'k--', lw=1)
        ax[i].plot(time[i], re[i], 'darkorange', lw=1)
        ax[i].set_xlabel("Time")
        ax[i].set_xlim(-500, 200)
        ax[i].set_title("GT")
        ax[i].set_ylabel('Strain (2,%d)'%(4*i/2))
        
    
    plt.tight_layout()
    plt.show()
    plt.close()

   

def gatech_strain_psi4_freqplots(y, w=-1):

    lm=[]
    for l in range(2,3):
        for m in range(0,l+1):
            lm.append((l,m))
    
    sim = scsearch(keyword=y.simname, verbose=True,institute='gt')[0]
    if w>0:
        y = gwylm(sim, lm, clean=True, w22=w)
        
    freq_21 = y.lm[2,1]['strain'].f
    amp_21 =  y.lm[2,1]['strain'].fd_amp

    
    freq_20 = y.lm[2,0]['strain'].f
    amp_20 =  y.lm[2,0]['strain'].fd_amp

    freq_22 = y.lm[2,2]['strain'].f
    amp_22 =  y.lm[2,2]['strain'].fd_amp

    psifreq_21 = y.lm[2,1]['psi4'].f
    psiamp_21 =  y.lm[2,1]['psi4'].fd_amp

    psifreq_20 = y.lm[2,0]['psi4'].f
    psiamp_20 =  y.lm[2,0]['psi4'].fd_amp

    psifreq_22 = y.lm[2,2]['psi4'].f
    psiamp_22 =  y.lm[2,2]['psi4'].fd_amp

    figure, (ax1, ax2) = plt.subplots(2,1,figsize=(14, 12), sharex=True)
    ax1.loglog(freq_20, 2*np.sqrt(freq_20)*amp_20, color='blue', label="2,0")
    ax1.loglog(freq_21, 2*np.sqrt(freq_21)*amp_21, color='black', ls='--', label="2,1")
    ax1.loglog(freq_22, 2*np.sqrt(freq_22)*amp_22, color='darkorange',  label="2,2")

    ax2.loglog(psifreq_20, 2*np.sqrt(psifreq_20)*psiamp_20, color='blue', label="2,0")
    ax2.loglog(psifreq_21, 2*np.sqrt(psifreq_21)*psiamp_21, color='black', ls='--', label="2,1")
    ax2.loglog(psifreq_22, 2*np.sqrt(psifreq_22)*psiamp_22, color='darkorange',  label="2,2")

    ax2.set_xlabel('Frequency')
    ax1.set_ylabel('Strain Amplitude ')
    ax2.set_ylabel('Psi4 Amplitude ')

    ax1.legend()
    plt.show()
    plt.close()    
    

def cutoff_freqplots(str_sxs, y):
   
    #SXS Waveforms
    time_20, hp_20, hx_20 = np.loadtxt(str_sxs['20']).T
    time_21, hp_21, hx_21 = np.loadtxt(str_sxs['21']).T
    time_22, hp_22, hx_22 = np.loadtxt(str_sxs['22']).T
    
    #Gatech Waveforms
    lm = []
    for l in range(2,3):
        for m in range(0,l+1):
            lm.append((l,m))
    
    sim = scsearch(keyword=y.simname, verbose=True,institute='gt')[0]
    y_win = gwylm(sim, lm, clean=True)

    plt.figure(figsize=(10, 5), dpi= 80, facecolor='w', edgecolor='k')
    
    plt.figure(figsize=(12, 10), dpi= 80, facecolor='w', edgecolor='k')
    plt.subplot(2,3,1)
    plt.plot(time_20, amp(hp_20, hx_20), color='b', label='2,0 SXS')
    #plt.plot(ts_21, amps_21, color='k', label='2,1 SXS')
    #plt.plot(ts_22, amps_22, color='darkorange', label='2,2 SXS')
    plt.xlabel('Time')
    plt.ylabel('Amplitude(Strain)')
    plt.xlim(-200, 200)
    plt.legend(loc ='upper left')
    
    for i in range(2,7):
     
        blockPrint(); 
        y = gwylm(sim, lm=lm, clean=True, w22=(y_win.wstart + 2.*(i+1)/100.));
        enablePrint();
        print y_win.wstart;
        w = y_win.wstart+2.*(i+1)/100.
    
        time_20, hp_20, hx_20 = y.lm[2,0]['strain'].wfarr.T
        time_21, hp_21, hx_21 = y.lm[2,1]['strain'].wfarr.T
        time_22, hp_22, hx_22 = y.lm[2,2]['strain'].wfarr.T
    
        amp_22 = np.sqrt(hp_22**2. + hx_22**2)
        tmax = time_22[amp_22==np.amax(amp_22)]
        
        plt.subplot(2,3,i)
        plt.plot(time_20-tmax, amp(hp_20, hx_20), color='b', label='2,0 GT \n w=%.4f'%w)
        #plt.plot(time_21-tmax, amp_21, color='k', label='2,1 GT')
        #plt.plot(time_22-tmax, amp_22, color='darkorange', label='2,2 GT' )
        plt.xlabel('Time')
        plt.ylabel('Amplitude(Strain)')
        plt.xlim(-200, 200)
        plt.legend(loc='upper left')
        
    plt.tight_layout()    
    plt.show()
    plt.close()