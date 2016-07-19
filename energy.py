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
import O1Solve
import O2Solve

from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
plt.rc('font',family='Times New Roman')
plt.rc('font', size=16)


#-----------------------------------------------------------------#
#Initialize vectors and physics parameters

param=pc.PhysicsConstants(1.0)

pg=PMNSGen.PMNSGen(param)
ug=SUGen.SUGen(param)


numcomp=True
#numcomp=False

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

valvec=np.logspace(3,9,100)
#valvec=np.linspace(0,pi,100)

num_vec=np.zeros(len(valvec))
diag_vec=np.zeros(len(valvec))
fracvec=np.zeros(len(valvec))
O1_vec=np.zeros(len(valvec))
O2_vec=np.zeros(len(valvec))
pert_vec0=np.zeros(len(valvec))
pert_vec1=np.zeros(len(valvec))
pert_vec2=np.zeros(len(valvec))

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

	dcy_switch=True
	mtr_switch=False
	potential=0

	pg.lamb=np.zeros((2,2))
	pg.lamb[0,1]=(-1.0)*theta
	
	ug.lamb=np.zeros((2,2))
	ug.lamb[0,1]=(-1.0)*phi




	#diag_amp = sp.DiagProb(0,0,E,param.PI,param,pg,ug,eig_dcy,dcy_switch,mtr_switch)
	diag_amp = sp.DiagProb(0,0,E,param.PI,param,pg,ug,eig_dcy,dcy_switch,mtr_switch,2*param.EARTHRADIUS*param.km)
	#diag_amp = sp.DiagProb(0,0,E,param.PI,param,pg,ug,eig_dcy,dcy_switch,mtr_switch)
	print "DIAGAMP: ",diag_amp
	diag_vec[val]=diag_amp


	print "PERTURBLEN: ",2*param.EARTHRADIUS*param.km


	O1_solve = O1Solve.O1Solve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1],mtr_switch)
	#O1_amp = O1_solve.P_ee(2*param.EARTHRADIUS*param.km)
	#O1_vec[val] = O1_amp

	O2_solve = O2Solve.O2Solve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1],mtr_switch)
	#O2_amp = O2_solve.P_ee(2*param.EARTHRADIUS*param.km)
	#O2_vec[val] = O2_amp
	
	num_amp = sp.AtmosphericNeutrinoOscillationProbability(0,0,E,param.PI,param,pg,ug,eig_dcy,dcy_switch,mtr_switch)
	num_vec[val] = num_amp 

	pert_solve0 = PerturbSolve.PerturbSolve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1],0,mtr_switch)
	pert_amp0 = pert_solve0.P_ee(2*param.EARTHRADIUS*param.km)
	pert_vec0[val] = pert_amp0

	pert_solve1 =PerturbSolve.PerturbSolve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1],1,mtr_switch)
	pert_amp1 = pert_solve1.P_ee(2*param.EARTHRADIUS*param.km)
	pert_vec1[val] = pert_amp1

	pert_solve2 = PerturbSolve.PerturbSolve(param.dm21sq,E,potential,theta,phi,eig_dcy[0],eig_dcy[1],2,mtr_switch)
	pert_amp2 = pert_solve2.P_ee(2*param.EARTHRADIUS*param.km)
	pert_vec2[val] = pert_amp2

	delta=param.dm21sq/(2*E)


#-----------------------------------------------------------------#
#Plotting data

axes=[]

axes.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=4))
if (numcomp==False):
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),pert_vec0,color="red",label="Standard Vacuum Solution",lw=2,ls='-')
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),pert_vec1,color="green",label="",lw=2,ls='--')
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),pert_vec2,color="green",label="",lw=2,ls='-.')
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),O1_vec,color="green",label="Perturbation O1",lw=2,ls='-')
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),O2_vec,color="blue",label="Perturbation O2",lw=2,ls='-')
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),diag_vec,color="black",label="Truncation O1",lw=0,ls='--')
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),diag_vec,color="black",label="Truncation O2",lw=0,ls='-.')
axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),num_vec,color="purple",label="Numerical Solution",lw=2,ls='-')

if (numcomp==True):
	axes[0].plot(np.multiply(valvec,param.GeV/param.TeV),diag_vec,color="orange",label="Diagonal Solution",lw=2,ls='--')

axes[0].set_xlim([np.multiply(valvec,param.GeV/param.TeV)[0],np.multiply(valvec,param.GeV/param.TeV)[-1]])
axes[0].set_ylim([0.5,1.1])
axes[0].set_xlabel("Neutrino Energy [TeV]",labelpad=0)
axes[0].set_ylabel("Muon Neutrino Survival Fraction")
axes[0].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

axes.append(plt.subplot2grid((4,4),(0,0),colspan=4,sharex=axes[0]))

if (numcomp==False):
	axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(pert_vec0,num_vec),num_vec),color="red",lw=2,ls='-')
	axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(pert_vec1,num_vec),num_vec),color="green",lw=2,ls='--')
	axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(pert_vec2,num_vec),num_vec),color="green",lw=2,ls='-.')
	axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(O1_vec,num_vec),num_vec),color="green",lw=2,ls='-')
#axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(diag_vec,num_vec),num_vec),color="pink",lw=2,ls='-')
	axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(O2_vec,num_vec),num_vec),color="blue",lw=2,ls='-')
axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(num_vec,num_vec),num_vec),color="purple",lw=2,ls='-')

if (numcomp==True):
	axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(num_vec,num_vec),num_vec),color="orange",lw=2,ls='--')

#axes[1].set_yscale('log')
axes[1].set_xlim([np.multiply(valvec,param.GeV/param.TeV)[0],np.multiply(valvec,param.GeV/param.TeV)[-1]])
axes[1].set_ylabel("Relative Residual")
axes[1].set_ylim([-0.2,0.2])
axes[1].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
plt.setp(axes[1].get_xticklabels(), visible=False)

for ax in axes:
	ax.set_xscale('log')



legends=[]
frames=[]

single_axis=[]
single_axis.append(axes[0])

for ax in enumerate(single_axis):

	legends.append(ax[1].legend(loc='lower right', shadow=False,fancybox=True))

	# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
	frames.append(legends[ax[0]].get_frame())
	frames[ax[0]].set_facecolor('1.0')

	# Set the fontsize
	for label in legends[ax[0]].get_texts():
	    label.set_fontsize(16)
	
	for label in legends[ax[0]].get_lines():
	    label.set_linewidth(2.0)  # the legend line width

for ax in axes:
	ax.tick_params(axis='x', which='major', width=1,length=5)
	ax.tick_params(axis='x', which='minor', width=1,length=3)
	ax.tick_params(axis='y', which='major', width=1,length=5)
	ax.tick_params(axis='y', which='minor', width=1,length=3)


# Get the Bbox
plt.draw()
bb = legends[0].legendPatch.get_bbox().inverse_transformed(axes[0].transAxes)

# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='cyan', alpha=0.2)


def extract_man_exp(val):
	#Extract mantissae and exponents for pretty printing!
	mystr="%.2e"%(val)
	mantissa=[]
	exponent=[]
	print mystr
	hit_e=False
	for c in mystr:
		if (c=='e'):
			hit_e=True
		else:
			if (hit_e==False):
				mantissa.append(c)
			else:
				exponent.append(c)
	man_str=''.join(mantissa)
	exp_str=''.join(exponent)
	
	ret=[float(man_str),int(exp_str)]
	return ret
		

def pretty_exp_note(parstring,unitstring,val,f_force=False,last=False):
	if (('e' not in '%g'%(val)) or (f_force==True)):
		ret = parstring+r'$=%.2f$'%(val)+'  ['+unitstring+']'
	else:
		splitval=extract_man_exp(val)
		ret = parstring+r'$=%.2f \times 10^{%d}$'%(splitval[0],splitval[1])+' ['+unitstring+']'
	if (last==False):
		ret+='\n'
	return ret



# place a text box relative to the legend

textlist=[
r'$\phi=\pi/4$  [rad]'+'\n', #Special value for phi
pretty_exp_note(r'$\theta$','rad',theta,f_force=True),
pretty_exp_note(r'$V$','TeV',potential/param.TeV),
pretty_exp_note(r'$\lambda_1$','TeV',eig_dcy[0]/param.TeV),
pretty_exp_note(r'$\lambda_2$','TeV',eig_dcy[1]/param.TeV,last=True)]
textstr=''.join(textlist)


axes[0].text(0.25, 0.37, textstr, transform=axes[0].transAxes, fontsize=16,verticalalignment='top', bbox=props)


plt.show()


#np.savez("4TeVnodcy_matter"+str(cmdargs[0])+"timeplot"+str(cmdargs[1]),d=distvec,p=pvec)
