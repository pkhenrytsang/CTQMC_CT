# Adapted by Pak Ki Henry Tsang Jun 2020 from Author: Kristjan Haule, March 2007-2017

from scipy import *
from fileio import *
import os,sys,subprocess
from scipy.optimize import curve_fit

def CreateInputFile(params,prefix):
	import os
	" Creates input file (PARAMS) for CT-QMC solver"
	CreateFolder(prefix)
	path = os.path.join(prefix,'PARAMS')  
	f = open(path, 'w')
	print >> f, '# Input file for continuous time quantum Monte Carlo'
	for p in params:
		print >> f, p, params[p][0], '\t', params[p][1]
	f.close()
    

def Diff(fg1, fg2):
    data1 = loadtxt(fg1).transpose()
    data2 = loadtxt(fg2).transpose()
    diff = sum(abs(data1-data2))/(shape(data1)[0]*shape(data1)[1])
    return diff

def Diff0(fg1, fg2):
    data1 = loadtxt(fg1).transpose()
    data2 = loadtxt(fg2).transpose()
    diff0 = data2[2][0]-data1[2][0]
    return diff0

def RaileighRitz(fg1,fg2,fg3):
    data1 = loadtxt(fg1)
    data2 = loadtxt(fg2)
    data3 = loadtxt(fg3)
    import numpy as np
    return 1-sum(np.sqrt(np.sum((data3-data2)**2,axis=1))*np.sqrt(np.sum((data2-data1)**2,axis=1)))/sum(np.sqrt(np.sum((data3-data2)**2,axis=1)))
    #return 1-sum(abs(data3-data2)*abs(data2-data1))/sum(abs(data3-data2))**2
    

def CalcSlope(fSig):

    Sigma = loadtxt(fSig).T

    def func(t,a,b,c):
        return a*(t-b)**2+c        
    try: #fit for quadratic curve using 3 points
        popt, pcov = curve_fit(func, Sigma[0][0:3], Sigma[2][0:3])
        slope = -2*popt[0]*popt[1]
        print "slope (quadratic) = ",slope
    except: #find linear slope (exact)
        m = (Sigma[2][1]-Sigma[2][0])/(Sigma[0][1]-Sigma[0][0])
        c = Sigma[2][1] - m*Sigma[0][1]
        slope = m
        print "slope (linear) = ",slope

    return slope
