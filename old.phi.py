#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import SUGen
import random
import sys



from pylab import rcParams
rcParams['figure.figsize'] = 10, 8

plt.rc('font', size=18)


import PerturbSolve
import ZeroSolve


import matplotlib.pyplot as plt

param=pc.PhysicsConstants(1.0)
#energy=np.logspace(10,12,10)
#energy=np.linspace(10,1000,2)
energy=np.logspace(np.log10(100),np.log10(100000),200)
#theta=np.linspace(3.141592/2,3.141592,100)
probs=np.zeros(len(energy))
ptots=np.zeros(len(energy))




ntype=0

pi=3.141592
#valvec=np.logspace(-20,-14,500)
valvec=np.linspace(0,2*pi,1000)
ampvec=np.zeros(len(valvec))
fracvec=np.zeros(len(valvec))
appampvec=np.zeros(len(valvec))
zeroampvec=np.zeros(len(valvec))

print len(valvec)
for val in range(0,len(valvec)):

	phi=valvec[val]
	#phi=pi/4.0
	theta=0.2318
	#theta=valvec[val]
	#E=valvec[val]*param.GeV
	E=4000*param.GeV

	eig_dcy=np.zeros(2)

	#eig_dcy[1]=0
	#eig_dcy[0]=valvec[val]
	#eig_dcy[0]=valvec[val]
	eig_dcy[0]=1e-16
	eig_dcy[1]=5e-16

	pg=PMNSGen.PMNSGen(param)
	pg.lamb=np.zeros((2,2))
	pg.lamb[0,1]=(-1.0)*theta
	
	
	ug=SUGen.SUGen(param)
	ug.lamb[0,1]=phi

	if ntype==0:
		pg.param.neutype='neutrino'
		ug.param.neutype='neutrino'
		param.neutype='neutrino'
	elif ntype==1:
		pg.param.neutype='antineutrino'
		ug.param.neutype='antineutrino'
		param.neutype='antineutrino'
	else:
		print "BAD NEUTRINO TYPE"
	
	
	
	amp=sp.AtmosphericNeutrinoOscillationProbability(0,0,E,param.PI,param,pg,ug,eig_dcy)
	psolve=PerturbSolve.PerturbSolve(param.dm21sq,E,theta,phi,eig_dcy[0],eig_dcy[1])
	approxamp= psolve.P_ee(2*param.EARTHRADIUS*param.km)
	print "Diameter, APP: ",2*param.EARTHRADIUS*param.km
	zerosolve=ZeroSolve.ZeroSolve(param.dm21sq,E,theta,phi,eig_dcy[0],eig_dcy[1])
	zeroamp= zerosolve.P_ee(2*param.EARTHRADIUS*param.km)
	zeroampvec[val]=zeroamp


	delta=param.dm21sq/(2*E)

	#fracvec[val]=valvec[val]/delta	
	ampvec[val]=amp
	appampvec[val]=approxamp
	print "VAL: ",valvec[val]
	print "FRACVAL: ",fracvec[val]
	print approxamp

fig,ax=plt.subplots()



plt.plot(valvec,zeroampvec,color="red",label="Perturbative: Order 0",lw=2)
plt.plot(valvec,appampvec,color="green",label="Perturbative: Order 1",lw=2)
plt.plot(valvec,ampvec,color="blue",label="Numerical",lw=2)


ax.set_ylim([0,1])
#ax.set_ylim([0.6,1.0])
#ax.set_xlim([0,1.65e4])
ax.set_xlabel(r"$\phi$ (radians)")
ax.set_ylabel("Muon Neutrino Survival Probability")
#figure_title="Decay Eigenvalue 1"
#plt.title(figure_title, y=1.08)

#xlim(x.min()*1.1, x.max()*1.1)
#xticks([(-2 * np.pi), (-3 * np.pi/2), -np.pi, -np.pi/2, 0, np.pi/2, np.pi, (3 * np.pi/2), (2 * np.pi)], [r'$-2\pi$', r'$-3/2\pi$', r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$', r'$3/2\pi$', r'$2\pi$'])
plt.xticks([0, np.pi/2, np.pi, (3 * np.pi/2), (2 * np.pi)], [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

ax.set_xlim([valvec[0],valvec[-1]])


#ax.set_xscale('log')





#legend = ax.legend(bbox_to_anchor=(0.05, 0.25), bbox_transform=plt.gcf().transFigure, shadow=False)




legend = ax.legend(loc='lower right', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(20)

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

#ax.tick_params(axis='x', which='major', labelsize=24)
#ax.tick_params(axis='y', which='major', labelsize=24)


textstr = '$\Delta=%.2e$  (eV)\n$E=%.2e$  (eV)\n'%(delta,E)+r'$\theta=%.2f$  (rad)'%(theta)+'\n$\lambda_1=%.2e$  (eV)\n$\lambda_2=%.2e$  (eV)'%(eig_dcy[0], eig_dcy[1])


# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='cyan', alpha=0.2)


# place a text box in upper left in axes coords
ax.text(0.035, 0.29, textstr, transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props)




plt.show()
#np.savez("4TeVnodcy_matter"+str(cmdargs[0])+"timeplot"+str(cmdargs[1]),d=distvec,p=pvec)
