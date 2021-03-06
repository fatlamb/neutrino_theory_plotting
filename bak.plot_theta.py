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

param=pc.PhysicsConstants(1.0)


if (len(sys.argv)==1):
	numcomp=False
else:
	if (sys.argv[1]=="numcomp"):
		numcomp=True
	else:
		print "BAD INPUT"

data=np.load('theta_nomatter.npz')


phi=data['_phi']
theta=data['_theta']
l1=data['_l1']
l2=data['_l2']
E=data['_E']
potential=data['_potential']
valvec=data['_valvec']
num_vec=data['_num_vec']
diag_vec=data['_diag_vec']
frac_vec=data['_frac_vec']
O1_vec=data['_O1_vec']
O2_vec=data['_O2_vec']
pert_vec0=data['_pert_vec0']
pert_vec1=data['_pert_vec1']
pert_vec2=data['_pert_vec2']
delta=data['_delta']
ntype=data['_ntype']


pi=3.141592
#-----------------------------------------------------------------#
#Plotting data


axes=[]

xax_vector=valvec

axes.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=4))
if (numcomp==False):
	std_vac=	axes[0].plot(xax_vector,pert_vec0,color="red",label="Standard Vacuum Solution",lw=2,ls='-')
	axes[0].plot(xax_vector,pert_vec1,color="green",label="",lw=2,ls='--')
	axes[0].plot(xax_vector,pert_vec2,color="green",label="",lw=2,ls='-.')
	pert1=	axes[0].plot(xax_vector,O1_vec,color="green",label="Perturbation O1",lw=2,ls='-')
	pert2=	axes[0].plot(xax_vector,O2_vec,color="blue",label="Perturbation O2",lw=2,ls='-')
	trunc1=	axes[0].plot(xax_vector,diag_vec,color="black",label="Truncation O1",lw=0,ls='--')
	trunc2=	axes[0].plot(xax_vector,diag_vec,color="black",label="Truncation O2",lw=0,ls='-.')
numsol=axes[0].plot(xax_vector,num_vec,color="purple",label="Numerical Solution",lw=2,ls='-')

if (numcomp==True):
	diagsol=axes[0].plot(xax_vector,diag_vec,color="orange",label="Diagonal Solution",lw=2,ls='--')

axes[0].set_xlim([xax_vector[0],xax_vector[-1]])
axes[0].set_ylim([0.5,1.0])
axes[0].set_xlabel(r"$\theta$ (radians)")
axes[0].set_ylabel("Muon Neutrino Survival Fraction")
axes[0].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

axes.append(plt.subplot2grid((4,4),(0,0),colspan=4,sharex=axes[0]))

if (numcomp==False):
	axes[1].plot(xax_vector,np.divide(np.subtract(pert_vec0,num_vec),num_vec),color="red",lw=2,ls='-')
	axes[1].plot(xax_vector,np.divide(np.subtract(pert_vec1,num_vec),num_vec),color="green",lw=2,ls='--')
	axes[1].plot(xax_vector,np.divide(np.subtract(pert_vec2,num_vec),num_vec),color="green",lw=2,ls='-.')
	axes[1].plot(xax_vector,np.divide(np.subtract(O1_vec,num_vec),num_vec),color="green",lw=2,ls='-')
#axes[1].plot(xax_vector,np.divide(np.subtract(diag_vec,num_vec),num_vec),color="pink",lw=2,ls='-')
	axes[1].plot(xax_vector,np.divide(np.subtract(O2_vec,num_vec),num_vec),color="blue",lw=2,ls='-')
axes[1].plot(xax_vector,np.divide(np.subtract(num_vec,num_vec),num_vec),color="purple",lw=2,ls='-')

if (numcomp==True):
	axes[1].plot(xax_vector,np.divide(np.subtract(num_vec,num_vec),num_vec),color="orange",lw=2,ls='--')

#axes[1].set_yscale('log')
axes[1].set_xlim([xax_vector[0],xax_vector[-1]])
axes[1].set_ylabel("Relative Residual")
axes[1].set_ylim([-0.2,0.2])
axes[1].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
plt.setp(axes[1].get_xticklabels(), visible=False)

plt.xticks([0, np.pi/2, np.pi, (3 * np.pi/2), (2 * np.pi)], [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])





legends=[]
frames=[]

if(numcomp==True):
	legends.append(axes[0].legend(handles=[numsol[0],diagsol[0]],loc='lower right', shadow=False,fancybox=True))
if(numcomp==False):
	legends.append(axes[0].legend(handles=[std_vac[0],pert1[0],pert2[0],numsol[0]],loc='lower right', shadow=False,fancybox=True))
	#legends.append(axes[0].legend(handles=[trunc1[0],trunc2[0]],loc='lower left', shadow=False,fancybox=True))

for leg in range(0,len(legends)):
	if (leg<len(legends)-1):
		axes[0].add_artist(legends[leg])
	# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
	frames.append(legends[leg].get_frame())
	frames[leg].set_facecolor('1.0')
	
	# Set the fontsize
	for label in legends[leg].get_texts():
	    label.set_fontsize(16)
	
	for label in legends[leg].get_lines():
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
#pretty_exp_note(r'$\theta$','rad',theta,f_force=True),
pretty_exp_note(r'$V$','TeV',potential/param.TeV),
pretty_exp_note(r'$E$','TeV',E/param.TeV),
pretty_exp_note(r'$\lambda_1$','TeV',l1/param.TeV),
pretty_exp_note(r'$\lambda_2$','TeV',l2/param.TeV,last=True)]
textstr=''.join(textlist)


axes[0].text(0.03, 0.37, textstr, transform=axes[0].transAxes, fontsize=16,verticalalignment='top', bbox=props)


plt.show()


#np.savez("4TeVnodcy_matter"+str(cmdargs[0])+"timeplot"+str(cmdargs[1]),d=distvec,p=pvec)
