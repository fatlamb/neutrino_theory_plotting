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


class O1Solve(object):

	def __init__(self,dm2,E,V,theta,phi,l1,l2,matter):

		delta=dm2/(2*E)

		if (matter==False):
			eta=theta
			xi=delta			
		elif (matter==True):
			eta=0.5*math.atan((delta*math.sin(2*theta))/(delta*math.cos(2*theta)-V))
			xi=math.sqrt((delta*math.cos(2*theta)-V)**2+(delta**2)*(math.sin(2*theta))**2)
		
		self.par=dict(zip(['eta','phi','xi','l1','l2'],[eta,phi,xi,l1,l2]))

		print "PARAMS"
		print self.par

	def P_ee(self,x):

		return matcall.math_eval(self.par['eta'],self.par['phi'],self.par['xi'],self.par['l1'],self.par['l2'],x,"o1_eval.m")

