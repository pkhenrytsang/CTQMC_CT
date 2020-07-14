#!/usr/bin/env python
#Python Wrapper for Kristjan Haule's CTQMC code (March 2007-2017). To solve the CT Model by Pak Ki Henry Tsang, Jun 2020

import argparse

parser = argparse.ArgumentParser()
#Physical Parameters
parser.add_argument("--T", help="temperature", required=True)
parser.add_argument("--U", help="Hubbard U", required=True)
parser.add_argument("--mu", help="chemical potential", required=True)
parser.add_argument("--ef", help="f-site energy", required=True)
parser.add_argument("--t", help="c-site hopping", required=True)
parser.add_argument("--V", help="c- f- orbital hybridization", required=True)

parser.add_argument("--gridsize", help="omega-grid size", required=True)

parser.add_argument("--cmix", help="mixing parameter (0.0-1.0)")
parser.add_argument("--Msteps", help="number of Monte Carlo steps")
parser.add_argument("--tsample", help="tsample")

parser.add_argument("--mode", help="Select initial conditions", choices=['fromNI','fromIC','refine'], required=True)

parser.add_argument("--Niter", help="number of iterations")
parser.add_argument("--NC", help="number of cores", nargs='?', type=int, const=1, default=1)
parser.add_argument("--nom", help="number of Matsubara points to sample")
parser.add_argument("--adiabatic", help="adiabatic switching", choices=['true', 'false', 'True','False','1','0','T','F'])
parser.add_argument("--spath", help="path to store files (under results/)", required=True)
args = parser.parse_args()


#Note that the following parameters are "global"
mu = float(args.mu)
print "Working on mu = ", mu

T = float(args.T)
print "Working on T = ",T

U = float(args.U)
print "Working on U = ",U

V = float(args.V)
print "Working on V = ", V

ef = float(args.ef)
print "Working on ef = ", ef

t = float(args.t)
print "Working on t = ", t

gridsize = int(args.gridsize)

if args.cmix:
    try:
        cmix = float(args.cmix)
        if not(cmix<1.0 and cmix>0.0):
            cmix = 1.0
            print "Invalid cmix input"
    except:
        cmix = 1.0
        print "Invalid cmix input"
    print "Mixing parameter set to : ",cmix
else:
    cmix = 1.0

Niterdefault = 10
# Number of DMFT iterations        
if args.Niter:
    try:
        Niter = int(args.Niter)
        if not(Niter>0):
            Niter = Niterdefault
            print "Invalid Niter input"
    except:
        Niter = Niterdefault
        print "Invalid Niter input"
    print "Niter set to : ",Niter
else:
    Niter = Niterdefault
    
tsampledefault = 200
# Number of DMFT iterations        
if args.tsample:
    try:
        tsample = int(args.tsample)
        if not(tsample>0):
            tsample = tsampledefault
            print "Invalid tsample input"
    except:
        tsample = tsampledefault
        print "Invalid tsample input"
    print "tsample set to : ",tsample
else:
    tsample = tsampledefault
    
nomdefault = 80
# Number of DMFT iterations        
if args.nom:
    try:
        nom = int(args.nom)
        if not(nom>0):
            nom = nomdefault
            print "Invalid nom input"
    except:
        nom = nomdefault
        print "Invalid nom input"
    print "nom set to : ",nom
else:
    nom = nomdefault
    
Mstepsdefault = 5e6
if args.Msteps:
    try:
        Msteps = float(args.Msteps)
        if not(Msteps>1):
            Msteps = Mstepsdefault
            print "Invalid Msteps input"
    except:
        Msteps = Mstepsdefault
        print "Invalid Msteps input"
    print "Msteps set to : ",Msteps
else:
    Msteps = Mstepsdefault
    
refine = False
fromNI = False
fromIC = False
if (args.mode=='refine'):
    refine = True
elif (args.mode=='fromIC'):
    fromIC = True
else:
    fromNI = True
    
if (args.adiabatic=='True' or args.adiabatic=='T' or args.adiabatic=='1'):
    adiabatic = True
else:
    adiabatic = False
    
    

##################        
#   Begin code   #
##################

from numpy import *
import matplotlib.pyplot as plt
import os,sys,subprocess
from fileio import *
from ctqmc_routines import *
from data import *
import glob,re
import scipy.signal as signal

ncores = int(args.NC)


icix="""# Cix file for cluster DMFT with CTQMC
# cluster_size, number of states, number of baths, maximum matrix size
1 4 2 1
# baths, dimension, symmetry, global flip
0       1 0 0
1       1 0 0
# cluster energies for unique baths, eps[k]
0 0
#   N   K   Sz size F^{+,dn}, F^{+,up}, Ea  S
1   0   0    0   1   2         3        0   0
2   1   0 -0.5   1   0         4        0   0.5
3   1   0  0.5   1   4         0        0   0.5
4   2   0    0   1   0         0        0   0
# matrix elements
1  2  1  1    1    # start-state,end-state, dim1, dim2, <2|F^{+,dn}|1>
1  3  1  1    1    # start-state,end-state, dim1, dim2, <3|F^{+,up}|1>
2  0  0  0
2  4  1  1   -1    # start-state,end-state, dim1, dim2, <4|F^{+,up}|2>
3  4  1  1    1
3  0  0  0
4  0  0  0
4  0  0  0
HB2                # Hubbard-I is used to determine high-frequency
# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]
0 0 0 0 0.0
# number of operators needed
0
"""

#############################
#   DMFT Self consistency   #
#############################

def DMFT_SCC(params,t,fDelta,cmix = 0.5):
    """This subroutine creates Delta.inp from Gc.out for DMFT on bethe lattice: Delta=t^2*Gc
    If Gc.out does not exist, it creates Gc.out which corresponds to the non-interacting model
    In the latter case also creates the inpurity cix file, which contains information about
    the atomic states.
    """
    fileGc = 'Gc.out'
    if (os.path.exists(fileGc)): # If output file exists, start from previous iteration
        # In the new Python, io.readarray is dropped and we should use loadtxt instead!
        Gc = loadtxt(fileGc)
    else: # otherwise start from non-interacting limit
        print 'Starting from non-interacting model (metal)'
        # Gc is metallic
        Gc=[]
        for n in range(gridsize):
            iom = (2*n+1)*pi/params['beta'][0]
            gc = 0.5*t**-2*1j*(iom-sqrt(iom**2+4*t**2))
            Gc.append([iom, gc.real, gc.imag])
        Gc = array(Gc)
        
        
    # creating impurity cix file
    filecix = params['cix'][0]
    if not(os.path.exists(filecix)):
        f = open(filecix, 'w')
        print >> f, icix
        f.close()
        
    if (os.path.exists(fDelta)): # If output file exists, start from previous iteration
        print "Mixing solutions"
        Mixing = True
        old_Delta = loadtxt(fDelta)
    else:
        Mixing = False
        old_Delta = array([])
        
    # Compute Deltaf
    iomega = (2*arange(len(Gc))+1)*pi*T
    Deltac = t**2*(Gc[:,1]+1.j*Gc[:,2])
    Deltaf = V**2/(1.j*iomega+mu-Deltac)
    
    #Save Deltac into file
    f = open('Deltac.inp', 'w')
    for i in xrange(len(Gc)):
        print >> f, Gc[i,0] , t**2*Gc[i,1], t**2*Gc[i,2]
    f.close()
        
    # Preparing input file Delta.inp
    f = open(fDelta, 'w')
    if ((len(Gc.T)==len(old_Delta.T)) and Mixing == True):
        for i in xrange(len(Gc)):
            print >> f, iomega[i], (1-cmix)*Deltaf.real[i]+cmix*Deltaf.real[i], (1-cmix)*Deltaf.imag[i]+cmix*Deltaf.imag[i]
    else:
        for i in xrange(len(Gc)):
            print >> f, iomega[i] , Deltaf.real[i], Deltaf.imag[i]
    f.close()
    
def Calc_n(omega,reG,T,points=20):

    def func(t,a,b):
            return a/(t**2+b**2)
    x=omega[-points:]
    y=reG[-points:]

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func, x, y)
    
    A=popt[0]
    C=popt[1]
    
    n=0
    
    n+=sum(reG)*2*T
    n+=A/C*(cosh(C/T)/sinh(C/T)-0.5*cosh(C/T/2)/sinh(C/T/2))-T*sum(2*A/(omega**2+C**2))
    n+=0.5
    
    return n
    
#Doped Hubbard Model Scan mu from insulator
beta = 1/T
print "MC steps = %d"%int(ncores*Msteps)
prefix = os.path.join("results",args.spath)
print "prefix = ",prefix
CreateFolder(prefix)

params = {"exe":   ["mpirun -np %d ctqmc"%ncores,          "# Path to executable"],
        "U":     [U,                 "# Coulomb repulsion (F0)"],
        "mu":    [mu-ef,              "# Chemical potential"],
        "beta":  [beta,                "# Inverse temperature"],
        "M" :    [Msteps,                "# Number of Monte Carlo steps"],
        "mode":  ["SM",               "# S stands for self-energy sampling, M stands for high frequency moment tail"],
        "cix":   ["one_band.imp",     "# Input file with atomic state"],
        "Delta": ["Delta.inp",        "# Input bath function hybridization"],
        "tsample":[tsample,               "# how often to record the measurements" ],
        "nom":   [nom,                 "# number of Matsubara frequency points to sample"],
      "svd_lmax":[30,                 "# number of SVD functions to project the solution"],
        "aom":   [1,                  "# number of frequency points to determin high frequency tail"],
        "GlobalFlip":[1000000,         "# how often to perform global flip"],
        }

if (fromIC==True): #Copy initial condition
    if (os.path.exists("initialcondition/Gc.out")):
        cmd = 'cp initialcondition/Gc.out run/Gc.out'
        subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying insulator Gc as initial condition
        print "mode: fromIC - Copying from initial condition"
        
# refine dmft run (need refine=True)
prev_it = 0
existpath = os.path.join(prefix,"Gc.out")

iters_list = []
G0_list = []
slope_list = []

if (refine and os.path.exists(existpath)):
    print "mode: refine - Continue DMFT from previously obtained data"
    cpfrompath = os.path.join(prefix,"*")
    cmd = 'cp %s run/'%cpfrompath
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) 

    fGcs = glob.glob(os.path.join(prefix,"Gc.out.*"))

    if (len(fGcs)>0):
        sfGcs= sorted(fGcs,key=lambda f: int(filter(str.isdigit, f)))[-1]
        prev_it = int(re.findall(r"\d+", sfGcs)[-1])+1
        
        it_ex,G0_list_ex = loadtxt(os.path.join(prefix,'G0.iter'), dtype=complex).T
        it_ex,slope_list_ex = loadtxt(os.path.join(prefix,'slope.iter')).T
        
        iters_list = it_ex.tolist()
        G0_list = G0_list_ex.tolist() 
        slope_list = slope_list_ex.tolist()
    
        
#clear existing data if not refining
if (not(refine) and os.path.exists(existpath)):
    cpfrompath = os.path.join(prefix,"*")
    cmd = 'rm %s'%cpfrompath
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) 
    
# Creating parameters file PARAMS for qmc execution
CreateInputFile(params,"run")

cwd=os.getcwd()
os.chdir("run")

for it in range(Niter):
    
    import time
    t0=time.time()

    # Constructing bath Delta.inp from Green's function
    DMFT_SCC(params,t,params['Delta'][0],cmix=cmix)
    
    # Running ctqmc
    print 'Running ---- qmc itt.: ', it+prev_it, '-----'

    subprocess.call(params['exe'][0], shell=True,stdout=sys.stdout,stderr=sys.stderr)
    
    # START PROCESS TO COMPUTE Gc
    
    #Get Gf
    om,reGf,imGf = loadtxt("Gf.out").T

    #Get Sigma
    om,reSigma,imSigma = loadtxt("Gf.out").T
    Sigma = reSigma+1.j*imSigma
    
    #Get Deltac
    om,reDeltac,imDeltac = loadtxt("Deltac.inp").T
    Deltac = reDeltac+1.j*imDeltac
    
    #Compute Phi
    Phi  = V**2/(1.j*om+mu-ef-Sigma)
    
    #Compute Gc
    Gc = 1/(1.j*om+mu-Phi-Deltac)
    reGc = Gc.real
    
    #Save Gc into file
    f = open('Gc.out', 'w')
    for i in xrange(len(Gc)):
        print >> f, om[i] , Gc.real[i], Gc.imag[i]
    f.close()
    
    G0 = reGf[0]+1.j*imGf[0]
    G0_list.append(G0)
    
    nf = Calc_n(om,reGf,T,points=20)*2
    nc = Calc_n(om,reGc,T,points=20)*2
    ntot = nf+nc
    print "mu = %f , ntot = %f \n nf = %f ,nc = %f"%(mu,ntot,nf,nc)
    
    # Some copying to store data obtained so far (at each iteration)
    cmd = 'cp Gf.out Gf.out.'+str(it+prev_it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying Gf
    
    cmd = 'cp Gc.out Gc.out.'+str(it+prev_it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying Gf

    cmd = 'cp Sig.out Sig.out.'+str(it+prev_it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying Sig

    cmd = 'cp ctqmc.log ctqmc.log.'+str(it+prev_it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying log file
    
    # Filter the results
    N  = 1    # Filter order
    Wn = 0.05 # Cutoff frequency
    B, A = signal.butter(N, Wn, output='ba')
    
    #if it>6:
    G0_diff_it,G0_diff = pdiff(abs(array(G0_list)))
    savetxt('diff.iter',vstack((G0_diff_it,G0_diff)).T )


    slope = CalcSlope("Sig.out")
    slope_list.append(slope)
    iters_list.append(it+prev_it)
    
    try:
        G0_diff_filtered = signal.filtfilt(B,A, G0_diff)
        savetxt('diff.filt.iter',vstack((G0_diff_it,G0_diff_filtered)).T )
        slope_list_filtered = signal.filtfilt(B,A, array(slope_list))
        savetxt('slope.filt.iter',vstack((array(iters_list),slope_list_filtered)).T)
        print "it %d : Diff = %le"%(it+prev_it , G0_diff_filtered[-1])
    except:    
        print "not enough data point to smoothen data"
        
    savetxt('slope.iter',vstack((array(iters_list),array(slope_list))).T)
    savetxt('G0.iter',vstack((array(iters_list),array(G0_list))).T)
    if slope > 0:
        print "Insulator found"
    else:
        Z = 1/(1-slope)
        print "Metal found, Z = %f"%Z
        
    print "Time elapsed : %fs"%(time.time()-t0)


os.chdir(cwd)

if (adiabatic):
    print "adiabatic - saving results to initialcondition/Gc.out"
    cmd = 'cp run/Gc.out initialcondition/Gc.out'
    os.system(cmd)

cmd = 'cp run/* %s'%prefix
os.system(cmd)

cmd = 'rm run/*'
os.system(cmd)

