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


class ZeroSolve(object):

	def __init__(self,dm2,E,theta,phi,l1,l2):

		delta=dm2/(2*E)
		ep=(l1-l2)/delta
		zeta=2*(theta-phi)
		b1=l1*math.cos(zeta/2.0)**2+l2*math.sin(zeta/2.0)**2
		b2=l2*math.cos(zeta/2.0)**2+l1*math.sin(zeta/2.0)**2

		
		self.par=dict(zip(['theta','delta','b1','b2','ep','zeta'],[theta,delta,b1,b2,ep,zeta]))


		print "PARAMS"
		print self.par
		print phi
		print "THETA: ",theta
	def P_ee(self,x):
	
		#return math.exp(-2*self.par['b1']*x)*math.cos(self.par['theta'])**4 + math.exp(-2*self.par['b2']*x)*math.sin(self.par['theta'])**4 + 0.5*math.exp(-1*(self.par['b1']+self.par['b2'])*x)*math.sin(2*self.par['theta'])*(math.sin(2*self.par['theta'])*math.cos(self.par['delta']*x) - 2*self.par['ep']*math.sin(self.par['zeta'])*math.sin(self.par['delta']*x))
		return 1.0-(math.sin(2*self.par['theta'])**2)*(math.sin(self.par['delta']*x/2)**2)
