simulation_names = {}

with open('wf_junkrad.txt','r') as f:
    for line in f:
        if '#' not in line:
            splitvals = line.split()
            simulation_names[splitvals[1]] = splitvals[0]

ids = []
tags = []
times = []
with open('Deb_wf_junkrad_after_strain.txt','r') as f: 
    for line in f:
        if '#' not in line:
            splitvals = line.split()
            tag = splitvals[0]
            time = splitvals[1]
            id_val = simulation_names[tag]
            ids.append(id_val)
            tags.append(tag)
            times.append(time)

with open('wf_junkrad_new.txt','w') as f:
    f.write("ID \t Tag \t Time  \n")
    for i in range(len(ids)):
        f.write( str(ids[i])+"\t"+str(tags[i])+"\t"+str(times[i])+"\n")
