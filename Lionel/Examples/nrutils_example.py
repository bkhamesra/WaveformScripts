# Import useful things
from nrutils import scsearch, gwylm, alert
# Find some simulations
alert('Finding simulations',header=True)
#A = scsearch(verbose=True,institute='gt', setname='q2_LL_D9_a0.6_th1_45_th2_225_m140')
A = scsearch(q=[3,5],nonspinning=True,verbose=True,institute='gt')
#print A
# Choose one, and print its properties
alert('Printing props for a single simulation',header=True)
a = A[0]; 
for k in a.__dict__: print('>> %s = %s'%(k,str(a.__dict__[k])))
# Load select multipoles from this simulation and get rid of junk 
alert('Load simulation data',header=True)
y = gwylm(a,lm=[(2,2),(2,1),(2,0)],verbose=True,clean=True)
# ... Note that the step above automatically calculates stain and news from psi4
# Plot the data in the time domain 
alert('Plotting data ...',header=True)
y.plot()
# Plot the data in the frequency domain
y.plot(domain='freq', show=True)
