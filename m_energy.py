#! /usr/bin/python
import matplotlib as mpl

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import SUGen
import random
import sys

import PerturbSolve
import MatterSolve
import ZeroMatterSolve
import ZeroSolve

import matplotlib.pyplot as plt

from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
plt.rc('font', size=14)


#-----------------------------------------------------------------#
#Initialize vectors and physics parameters

param=pc.PhysicsConstants(1.0)

pg=PMNSGen.PMNSGen(param)
ug=SUGen.SUGen(param)

ntype=1 #0-neutrino, 1-antineutrino

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

pi=3.141592

valvec=np.logspace(3,6,100)
#valvec=np.linspace(0,pi,100)

ampvec=np.zeros(len(valvec))
fracvec=np.zeros(len(valvec))
appampvec=np.zeros(len(valvec))
zeroampvec=np.zeros(len(valvec))
appampvec_nomatter=np.zeros(len(valvec))
zeroampvec_nomatter=np.zeros(len(valvec))


#-----------------------------------------------------------------#
#Calculate Probabilities
for val in range(0,len(valvec)):

	phi=pi/4.0
	theta=0.2318
	#theta=valvec[val]
	#phi=valvec[val]
	E=valvec[val]*param.GeV
	#E=4000*param.GeV

	eig_dcy=np.zeros(2)

	eig_dcy[0]=1e-16
	#eig_dcy[0]=valvec[val]
	eig_dcy[1]=5e-16
	#eig_dcy[1]=valvec[val]

	potential=0
	#potential=0

	pg.lamb=np.zeros((2,2))
	pg.lamb[0,1]=(-1.0)*theta
	
	ug.lamb=np.zeros((2,2))
    ug.lamb[0,1]=(-1.0)*phi


	amp=sp.AtmosphericNeutrinoOscillationProbability(0,0,E,param.PI,param,pg,ug,eig_dcy)
	psolve=MatterSolve.MatterSolve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1])
	approxamp= psolve.P_ee(2*param.EARTHRADIUS*param.km)
	psolve_nomatter=PerturbSolve.PerturbSolve(param.dm21sq,E,theta,phi,eig_dcy[0],eig_dcy[1])
	approxamp_nomatter= psolve_nomatter.P_ee(2*param.EARTHRADIUS*param.km)
	zerosolve=ZeroMatterSolve.ZeroMatterSolve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1])
	zeroamp= zerosolve.P_ee(2*param.EARTHRADIUS*param.km)
	zerosolve_nomatter=ZeroSolve.ZeroSolve(param.dm21sq,E,theta,phi,eig_dcy[0],eig_dcy[1])
	zeroamp_nomatter= zerosolve_nomatter.P_ee(2*param.EARTHRADIUS*param.km)


	delta=param.dm21sq/(2*E)

	ampvec[val]=amp
	zeroampvec[val]=zeroamp
	zeroampvec_nomatter[val]=zeroamp_nomatter
	appampvec[val]=approxamp
	appampvec_nomatter[val]=approxamp_nomatter


#-----------------------------------------------------------------#
#Plotting data

axes=[]

axes.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=4))
axes[0].plot(np.multiply(valvec,param.GeV),zeroampvec,color="red",label="Perturbative: Order 0",lw=2)
axes[0].plot(np.multiply(valvec,param.GeV),appampvec,color="green",label="Perturbative: Order 1",lw=2)
axes[0].plot(np.multiply(valvec,param.GeV),appampvec_nomatter,color="green",label="Perturbative: Order 1",lw=2,ls="--")
axes[0].plot(np.multiply(valvec,param.GeV),ampvec,color="blue",label="Numerical",lw=2)
axes[0].set_xlim([np.multiply(valvec,param.GeV)[0],np.multiply(valvec,param.GeV)[-1]])
axes[0].set_ylim([0,1])
axes[0].set_xlabel("Energy (eV)")
axes[0].set_ylabel("Muon Neutrino Survival Probability")

axes.append(plt.subplot2grid((4,4),(0,0),colspan=4,sharex=axes[0]))
axes[1].plot(np.multiply(valvec,param.GeV),np.absolute(np.subtract(ampvec,appampvec)),color="black",lw=2)
axes[1].set_yscale('log')
axes[1].set_xlim([np.multiply(valvec,param.GeV)[0],np.multiply(valvec,param.GeV)[-1]])
axes[1].set_ylabel("Abs. Residual (Num-O1)")
plt.setp(axes[1].get_xticklabels(), visible=False)

for ax in axes:
	ax.set_xscale('log')



legends=[]
frames=[]

single_axis=[]
single_axis.append(axes[0])

for ax in enumerate(single_axis):

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


# Get the Bbox
plt.draw()
bb = legends[0].legendPatch.get_bbox().inverse_transformed(axes[0].transAxes)

# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='cyan', alpha=0.2)

# place a text box relative to the legend
textstr = '$\phi=%.2f$  (rad)\n'%(phi)+r'$\theta=%.2f$  (rad)'%(theta)+'\n$\lambda_1=%.2e$  (eV)\n$\lambda_2=%.2e$  (eV)\n$V=%.2e$'%(eig_dcy[0], eig_dcy[1],potential)


axes[0].text(0.665, 0.675, textstr, transform=axes[0].transAxes, fontsize=20,verticalalignment='top', bbox=props)


plt.show()


#np.savez("4TeVnodcy_matter"+str(cmdargs[0])+"timeplot"+str(cmdargs[1]),d=distvec,p=pvec)
