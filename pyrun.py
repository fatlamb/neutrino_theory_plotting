#! /usr/bin/python

import subprocess


def math_eval(theta,phi,xi,l1,l2,t,scriptname):

	arglist=[theta,phi,xi,l1,l2,t]
	arg_str_list=[]
	
	for arg in arglist:
		arg_str_list.append(' '+str(arg))
	
	arg_str=''.join(arg_str_list)
	print arg_str
	
	input_string = "./"+scriptname+arg_str
	print input_string
	
	out = subprocess.check_output([input_string], shell=True)
	return float(out)
