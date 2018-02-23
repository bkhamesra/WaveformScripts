
import numpy as np
import os
import warnings
import time
	
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


def initial_data(dirpath):

	datadir = os.path.join(dirpath, "data")
	filename = dirpath.split("/")[-1]
    	parfile = os.path.join(datadir, filename + '.par')

	if not os.path.isfile(parfile) :
		raise ValueError(" Parfile missing in the directory. Please check again. \n")
	
	r1 = np.empty(3)
	r2 = np.empty(3)
	p1 = np.empty(3)
	p2 = np.empty(3)
	spin1 = np.empty(3)
	spin2 = np.empty(3)

	mplus = output_data(parfile, 'par_m_plus')
	mminus = output_data(parfile, 'par_m_minus')
	
	par_b = output_data(parfile, 'par_b')
 
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
			
	spin1[0] = output_data(parfile, 'par_s_plus[0]')
	spin1[1] = output_data(parfile, 'par_s_plus[1]')
	spin1[2] = output_data(parfile, 'par_s_plus[2]')
	spin2[0] = output_data(parfile, 'par_s_minus[0]')
	spin2[1] = output_data(parfile, 'par_s_minus[1]')
	spin2[2] = output_data(parfile, 'par_s_minus[2]')
	
	nr_initdata = {}
#	nr_initdata[('BH1')] = {'puncture_mass': mplus, 'position': r1, 'spin': spin1, 'momentum': p1}
#	nr_initdata[('BH2')] = {'puncture_mass': mminus, 'position': r2, 'spin': spin2, 'momentum': p2}

	nr_initdata['puncture_mass1'] = mplus
	nr_initdata['puncture_mass2'] = mminus
	nr_initdata['spin_BH1'] = spin1
	nr_initdata['spin_BH2'] = spin2
	nr_initdata['momentum_BH1'] = p1
	nr_initdata['momentum_BH2'] = p2
	nr_initdata['separation'] = 2.*par_b

	return nr_initdata
