import numpy as np
import os
import matplotlib.pyplot as plt

def orbital_frequency(sim_path):
    #Read in the data
    path1 = sim_path+"/data/ShiftTracker0.asc"
    path2 = sim_path+"/data/ShiftTracker1.asc"
    time1, x1, y1 = np.loadtxt(path1, usecols = (1,2,3), unpack = True)
    time2, x2, y2 = np.loadtxt(path2, usecols = (1,2,3), unpack = True)

    #Crop to the same length
    end_cut = min(len(time1),len(time2))
    time1 = time1[:end_cut]
    x1 = x1[:end_cut]
    y1 = y1[:end_cut]
    time2 = time2[:end_cut]
    x2 = x2[:end_cut]
    y2 = y2[:end_cut]
    
    #Calculate the frequency
    x_diff = x2-x1
    y_diff = y2-y1

    phase = np.unwrap(np.arctan2(y_diff, x_diff))
    frequency = np.gradient(phase)/np.gradient(time1)

    return (time1,frequency)

def gw_frequency(sim_path):
    #Read in the data
    path = sim_path+"/data/Strain/Strain_l2_m2.txt"
    time, phase = np.loadtxt(path, usecols = (0,4), unpack = True)
    
    #Calculate the frequency
    frequency = np.gradient(phase)/np.gradient(time)

    return (time,frequency)

def find_cutoff(sim_path):
    orbital_time, orbital_freq = orbital_frequency(sim_path)
    orbital_time = [x+75 for x in orbital_time]
    orbital_freq = 2*orbital_freq
    gw_time, gw_freq = gw_frequency(sim_path)

    #create new time
    time_min = max(orbital_time[0],gw_time[0])
    time_max = min(orbital_time[-1],gw_time[-1])
    timestep = max(orbital_time[1]-orbital_time[0], gw_time[1]-gw_time[0])
    times = np.arange(time_min, time_max, timestep)    
    
    #interpolate
    interp_orbital = np.interp(times, orbital_time, orbital_freq)
    interp_gw = np.interp(times, gw_time, gw_freq)

    #plt.plot(times, interp_orbital, label = "orbital")
    #plt.plot(times, interp_gw, label = "gw")
    #plt.legend()
    #plt.show()

    difference = abs(interp_gw - interp_orbital) 
    print(times[0:59])
    cutoff = 0

    for i in range(50,len(difference)):
#    while(True):
#        print(i)
#        print(difference[i-50:i+50])
        if(all([x < 0.0015 for x in difference[i-50:i+50]])):
            cutoff = i
            break;
#        i+=1
    return times[cutoff]

start_times = []
waveform_names = []

#path = "/numrel/NumRel/dferguson41/phys44236/localdata2/TEST_OUTPUT/input_dir/AlignedSpin"
path_aligned = "/numrel/NumRel/dferguson41/GENERATING_CATALOG/ascii_data/Remaining/AlignedSpin"
path_nonspinning = "/numrel/NumRel/dferguson41/GENERATING_CATALOG/ascii_data/Remaining/NonSpinning"
path_precessing = "/numrel/NumRel/dferguson41/GENERATING_CATALOG/ascii_data/Remaining/Precessing"
#locations= ["D11_a0.2_q1.00_m103_As"]

output_file = "wf_junkrad_unformatted.txt"
with open(output_file, 'w') as f:
    f.write("Tags \t Time of Junk \n")
#    for location in os.listdir(path_aligned):
#        tag = location
#        time_junk = find_cutoff(path_aligned+"/"+location) 
#        f.write(tag + "\t" + str(time_junk) +"\n")
#    for location in os.listdir(path_nonspinning):
#        tag = location
#        time_junk = find_cutoff(path_nonspinning+"/"+location) 
#        f.write(tag + "\t" + str(time_junk) +"\n")
    for location in os.listdir(path_precessing):
        tag = location
        time_junk = find_cutoff(path_precessing+"/"+location) 
        f.write(tag + "\t" + str(time_junk) +"\n")
        
