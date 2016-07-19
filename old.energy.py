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
energy=np.logspace(np.log10(100),np.log10(100000),10000)
#theta=np.linspace(3.141592/2,3.141592,100)
probs=np.zeros(len(energy))
ptots=np.zeros(len(energy))




ntype=0

pi=3.141592
valvec=np.logspace(3,6,10)
#valvec=np.linspace(0,2*pi,200)
ampvec=np.zeros(len(valvec))
fracvec=np.zeros(len(valvec))
appampvec=np.zeros(len(valvec))
zeroampvec=np.zeros(len(valvec))
print len(valvec)
for val in range(0,len(valvec)):

	#phi=valvec[val]
	phi=pi/4.0
	theta=0.2318
	#phi=theta
	#theta=valvec[val]
	E=valvec[val]*param.GeV
	#E=4000*param.GeV

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
	zerosolve=ZeroSolve.ZeroSolve(param.dm21sq,E,theta,phi,eig_dcy[0],eig_dcy[1])
	zeroamp= zerosolve.P_ee(2*param.EARTHRADIUS*param.km)
	zeroampvec[val]=zeroamp
	print zeroamp
	print "Diameter, APP: ",2*param.EARTHRADIUS*param.km

	delta=param.dm21sq/(2*E)

	#fracvec[val]=valvec[val]/delta	
	ampvec[val]=amp
	appampvec[val]=approxamp
	print "VAL: ",valvec[val]
	print "FRACVAL: ",fracvec[val]
	print approxamp

#fig,ax=plt.subplots()
#fig = plt.figure()
#fig, (ax, a2) = plt.subplots(1,2, gridspec_kw = {'height_ratios':[3, 1]})
axes=[]
axes.append(plt.subplot2grid((4,4),(0,0),rowspan=3,colspan=4))
axes.append(plt.subplot2grid((4,4),(3,0),colspan=4,sharex=axes[0]))
#ax = fig.add_subplot(413)
#ax2 = fig.add_subplot(414, sharex=ax)


axes[0].plot(np.multiply(valvec,param.GeV),zeroampvec,color="red",label="Perturbative: Order 0",lw=2)
axes[0].plot(np.multiply(valvec,param.GeV),appampvec,color="green",label="Perturbative: Order 1",lw=2)
axes[0].plot(np.multiply(valvec,param.GeV),ampvec,color="blue",label="Numerical",lw=2)
axes[1].plot(np.multiply(valvec,param.GeV),np.subtract(ampvec,appampvec),color="black",lw=2)


axes[0].set_ylim([0,1])
#ax2.set_ylim([0,1])
#ax.set_ylim([0.6,1.0])
#ax.set_xlim([0,1.65e4])
axes[1].set_xlabel("Energy (eV)")
axes[0].set_ylabel("Muon Neutrino Survival Probability")
axes[1].set_ylabel("Residual")
#figure_title="Decay Eigenvalue 1"
#plt.title(figure_title, y=1.08)


plt.setp(axes[0].get_xticklabels(), visible=False)

legends=[]
frames=[]

single_axis=[]
single_axis.append(axes[0])
for ax in enumerate(axes):
	ax[1].set_xlim([valvec[0]*param.GeV,valvec[-1]*param.GeV])
	ax[1].set_xscale('log')




for ax in enumerate(single_axis):
#legend = ax.legend(bbox_to_anchor=(0.05, 0.25), bbox_transform=plt.gcf().transFigure, shadow=False)



	legends.append(ax[1].legend(loc='lower right', shadow=False))

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
	frames.append(legends[ax[0]].get_frame())
	frames[ax[0]].set_facecolor('1.0')

	# Set the fontsize
	for label in legends[ax[0]].get_texts():
	    label.set_fontsize(20)
	
	for label in legends[ax[0]].get_lines():
	    label.set_linewidth(1.5)  # the legend line width

#ax.tick_params(axis='x', which='major', labelsize=24)
#ax.tick_params(axis='y', which='major', labelsize=24)


textstr = '$\phi=%.2f$  (rad)\n'%(phi)+r'$\theta=%.2f$  (rad)'%(theta)+'\n$\lambda_1=%.2e$  (eV)\n$\lambda_2=%.2e$  (eV)'%(eig_dcy[0], eig_dcy[1])


# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='cyan', alpha=0.2)


# place a text box in upper left in axes coords
axes[0].text(0.035, 0.325, textstr, transform=axes[0].transAxes, fontsize=20,
        verticalalignment='top', bbox=props)




plt.show()
#np.savez("4TeVnodcy_matter"+str(cmdargs[0])+"timeplot"+str(cmdargs[1]),d=distvec,p=pvec)
