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

import ZeroSolve
import PerturbSolve
import MatterSolve
import ZeroMatterSolve

import matplotlib.pyplot as plt

param=pc.PhysicsConstants(1.0)




ntype=1

pi=3.141592
valvec=np.logspace(2.5,5,1000)
#valvec=np.linspace(0,2*pi,200)
phivals=[0,1/pi,2/pi,3/pi]
ampvec=np.zeros((len(phivals),len(valvec)))
fracvec=np.zeros(len(valvec))
appampvec=np.zeros((len(phivals),len(valvec)))
appampvec_nomatter=np.zeros(len(valvec))
zeroampvec=np.zeros(len(valvec))
zeroampvec_nomatter=np.zeros(len(valvec))
print len(valvec)

for phi_n in range(0,len(phivals)):
	for val in range(0,len(valvec)):
	
		#phi=valvec[val]
		#phi=pi/4.0
		phi=phivals[phi_n]
		theta=0.2318
		#phi=theta
		#theta=valvec[val]
		E=valvec[val]*param.GeV
		#E=4000*param.GeV
	
		eig_dcy=np.zeros(2)
	
		#eig_dcy[1]=0
		#eig_dcy[0]=valvec[val]
		#eig_dcy[0]=valvec[val]
		eig_dcy[0]=0
		#eig_dcy[0]=1e-16
		eig_dcy[1]=1e-14
		#eig_dcy[1]=5e-16
		#eig_dcy[1]=0
	
		potential=2.65492195354e-13
	
	
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
		
		
		
		amp=sp.AtmosphericNeutrinoOscillationProbability(0,0,E,param.PI,param,pg,ug,eig_dcy,True,True)
		psolve=MatterSolve.MatterSolve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1])
		approxamp= psolve.P_ee(2*param.EARTHRADIUS*param.km)
		psolve_nomatter=PerturbSolve.PerturbSolve(param.dm21sq,E,theta,phi,eig_dcy[0],eig_dcy[1])
		approxamp_nomatter= psolve_nomatter.P_ee(2*param.EARTHRADIUS*param.km)
		zerosolve=ZeroMatterSolve.ZeroMatterSolve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1])
		zeroamp= zerosolve.P_ee(2*param.EARTHRADIUS*param.km)
		zerosolve_nomatter=ZeroSolve.ZeroSolve(param.dm21sq,E,theta,phi,eig_dcy[0],eig_dcy[1])
		zeroamp_nomatter= zerosolve_nomatter.P_ee(2*param.EARTHRADIUS*param.km)
		zeroampvec[val]=zeroamp
		zeroampvec_nomatter[val]=zeroamp_nomatter
		print zeroamp
		print "Diameter, APP: ",2*param.EARTHRADIUS*param.km
	
		delta=param.dm21sq/(2*E)
	
		#fracvec[val]=valvec[val]/delta	
		ampvec[phi_n,val]=amp
		appampvec[phi_n,val]=approxamp
		appampvec_nomatter[val]=approxamp_nomatter
		print "VAL: ",valvec[val]
		print "FRACVAL: ",fracvec[val]
		print approxamp

fig,ax=plt.subplots()


phi_cols=["red","orange","green","blue"]

plt.plot(np.multiply(valvec,param.GeV),zeroampvec,color="black",label="Perturbative: Order 0",lw=2,ls="-.")
plt.plot(np.multiply(valvec,param.GeV),zeroampvec,color="black",label="Perturbative: Order 1",lw=0,ls="-")
#plt.plot(np.multiply(valvec,param.GeV),zeroampvec,color="black",label="Numerical",lw=0,ls="-")
#plt.plot(np.multiply(valvec,param.GeV),zeroampvec_nomatter,color="red",label="Perturbative, No Matter: Order 0",lw=2,ls="--")
for phi_n in range(0,len(phivals)):
	phi_label=" "+r"$\phi=$"+str(phi_n)+r'$\pi/10$'
	plt.plot(np.multiply(valvec,param.GeV),appampvec[phi_n,:],color=phi_cols[phi_n],lw=2,ls="-")
#plt.plot(np.multiply(valvec,param.GeV),appampvec_nomatter,color="green",label="Perturbative, No Matter: Order 1",lw=2,ls="--")
	#plt.plot(np.multiply(valvec,param.GeV),ampvec[phi_n,:],color=phi_cols[phi_n],lw=2,ls="-")


ax.set_ylim([0,1])
#ax.set_ylim([0.6,1.0])
#ax.set_xlim([0,1.65e4])
ax.set_xlabel("Energy (eV)")
ax.set_ylabel("Muon Neutrino Survival Probability")
#figure_title="Decay Eigenvalue 1"
#plt.title(figure_title, y=1.08)


ax.set_xlim([valvec[0]*param.GeV,1e6*param.GeV])

ax.set_xscale('log')





#legend = ax.legend(bbox_to_anchor=(0.05, 0.25), bbox_transform=plt.gcf().transFigure, shadow=False)




legend = ax.legend(loc='lower right', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(18)

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

#ax.tick_params(axis='x', which='major', labelsize=24)
#ax.tick_params(axis='y', which='major', labelsize=24)


textstr = r'$\theta=%.2f$  (rad)'%(theta)+'\n$\lambda_1=%.2e$  (eV)\n$\lambda_2=%.2e$  (eV)'%(eig_dcy[0], eig_dcy[1])


# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='cyan', alpha=0.2)


# place a text box in upper left in axes coords
ax.text(0.825, 0.35, textstr, transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props)




plt.show()
#np.savez("4TeVnodcy_matter"+str(cmdargs[0])+"timeplot"+str(cmdargs[1]),d=distvec,p=pvec)
