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


class MatterSolve(object):

	def __init__(self,dm2,E,V,theta,phi,l1,l2):

		delta=dm2/(2*E)
		eta=0.5*math.atan((delta*math.sin(2*theta))/(delta*math.cos(2*theta)-V))	
		xi=math.sqrt((delta*math.cos(2*theta)-V)**2+(delta**2)*(math.sin(2*theta))**2)
		ep=(l1-l2)/xi
		zeta=2*(eta-phi)
		b1=l1*math.cos(zeta/2.0)**2+l2*math.sin(zeta/2.0)**2
		b2=l2*math.cos(zeta/2.0)**2+l1*math.sin(zeta/2.0)**2

		
		self.par=dict(zip(['theta','eta','xi','b1','b2','ep','zeta'],[theta,eta,xi,b1,b2,ep,zeta]))


		print "PARAMS"
		print self.par
		print phi
	def P_ee(self,x):

		return (1/(xi**2))*math.exp(-2*x*(self.par['


	
		return math.exp(-2*self.par['b1']*x)*math.cos(self.par['eta'])**4 + math.exp(-2*self.par['b2']*x)*math.sin(self.par['eta'])**4 + 0.5*math.exp(-1*(self.par['b1']+self.par['b2'])*x)*math.sin(2*self.par['eta'])*(math.sin(2*self.par['eta'])*math.cos(self.par['xi']*x) - 2*self.par['ep']*math.sin(self.par['zeta'])*math.sin(self.par['xi']*x))
		#return 1.0-(math.sin(2*self.par['theta'])**2)*(math.sin(self.par['delta']*x/2)**2)
