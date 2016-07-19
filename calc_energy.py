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

valvec=np.logspace(0,6,100)
#valvec=np.linspace(0,pi,100)

frac_vec=np.zeros(len(valvec))
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
	E=valvec[val]*param.TeV
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



np.savez("energy",_phi=phi,_theta=theta,_l1=eig_dcy[0],_l2=eig_dcy[1],_E=E,_potential=potential,_valvec=valvec,_num_vec=num_vec,_diag_vec=diag_vec,_frac_vec=frac_vec,_O1_vec=O1_vec,_O2_vec=O2_vec,_pert_vec0=pert_vec0,_pert_vec1=pert_vec1,_pert_vec2=pert_vec2,_delta=delta,_ntype=ntype)

