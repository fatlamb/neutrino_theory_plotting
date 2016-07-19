import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import random
import cmath
import math
import mathematica_caller as matcall


class PerturbSolve(object):

	def __init__(self,dm2,E,V,theta,phi,l1,l2,order,matter):

		delta=dm2/(2*E)

		if (matter==False):
			eta=theta
			xi=delta			
		elif (matter==True):
			eta=0.5*math.atan((delta*math.sin(2*theta))/(delta*math.cos(2*theta)-V))	
			xi=math.sqrt((delta*math.cos(2*theta)-V)**2+(delta**2)*(math.sin(2*theta))**2)

		if (order==0):
			zeta=0
			ep=0		
			b1=0
			b2=0
		elif (order>0):
			zeta=2*(eta-phi)
			ep=(l1-l2)/xi
			b1=l1*math.cos(zeta/2.0)**2+l2*math.sin(zeta/2.0)**2
			b2=l2*math.cos(zeta/2.0)**2+l1*math.sin(zeta/2.0)**2
			

		if (order==2):
			O2=1.0
		elif (order!=2):
			O2=0.0
	


		
		self.par=dict(zip(['theta','eta','xi','b1','b2','ep','zeta','O2'],[theta,eta,xi,b1,b2,ep,zeta,O2]))

		print "PARAMS"
		print self.par

	def P_ee(self,x):
	
		D1 = math.exp(-2*self.par['b1']*x)*math.cos(self.par['eta'])**2*(math.cos(self.par['eta'])**2 + self.par['O2']*(math.sin(self.par['eta'])**2)*math.sin(self.par['zeta'])**2*self.par['ep']**2) 
		D2 = math.exp(-2*self.par['b2']*x)*math.sin(self.par['eta'])**2*(math.sin(self.par['eta'])**2 + self.par['O2']*(math.cos(self.par['eta'])**2)*math.sin(self.par['zeta'])**2*self.par['ep']**2)

		OSC = 0.5*math.exp(-1*(self.par['b1']+self.par['b2'])*x)*math.sin(2*self.par['eta'])*(math.sin(2*self.par['eta'])*math.cos(self.par['xi']*x)*(1-0.5*self.par['O2']*self.par['ep']**2*math.sin(self.par['zeta'])**2) - 2*self.par['ep']*math.sin(self.par['zeta'])*math.sin(self.par['xi']*x))

		return D1+D2+OSC

