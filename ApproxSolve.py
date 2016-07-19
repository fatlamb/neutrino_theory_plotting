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


class ApproxSolve(object):

	def __init__(self,theta,phi,l1,l2):

		ep=(l1-l2)/delta
		zeta=2*(theta-phi)
		b1=l1*math.cos(zeta)
		b2=l2*math.cos(zeta)

		
		self.par=dict(zip(['theta','delta','b1','b2','ep','zeta'],[theta,delta,b1,b2,ep,zeta]))


	def P_ee(self,x):
	
		return math.exp(-2*self.par['b1']*x)*math.cos(self.par['theta'])**4 + math.exp(-2*self.par['b2']*x)*math.sin(self.par['theta'])**4 + 0.5*math.exp(-1*(self.par['b1']+self.par['b2'])*x)*math.sin(2*self.par['theta'])*(math.sin(2*self.par['theta'])*math.cos(self.par['delta']*x) - 2*self.par['ep']*math.sin(self.par['zeta'])*math.sin(self.par['delta']*x))

