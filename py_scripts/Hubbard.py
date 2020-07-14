#!/usr/bin/env python
#Python Interface to Kristjan Haule's CTQMC code (March 2007-2017). To solve the Doped Hubbard Model by Pak Ki Henry Tsang, Jun 2020

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--T", help="temperature", required=True)
parser.add_argument("--U", help="Hubbard U", required=True)
parser.add_argument("--mu", help="chemical potential", required=True)
parser.add_argument("--cmix", help="mixing parameter (0.0-1.0)")
parser.add_argument("--Msteps", help="number of Monte Carlo steps")

parser.add_argument("--mode", help="Select initial conditions", choices=['fromNI','fromIC','refine'], required=True)

parser.add_argument("--Niter", help="number of iterations")
parser.add_argument("--NC", help="number of cores", nargs='?', type=int, const=1, default=1)
parser.add_argument("--nom", help="number of Matsubara points to sample")
parser.add_argument("--adiabatic", help="adiabatic switching", choices=['true', 'false', 'True','False','1','0','T','F'])
parser.add_argument("--spath", help="path to store files (under results/)", required=True)
args = parser.parse_args()


mu = float(args.mu)
print "Working on mu = ", mu

T = float(args.T)
print "Working on T = ",T

U = float(args.U)
print "Working on U = ",U

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
    
    

    
        
#Begin code

from numpy import *
import matplotlib.pyplot as plt
import os,sys,subprocess
from fileio import *
from ctqmc_routines import *
from data import *
import glob,re

ncores = int(args.NC)
#storeas = "DopedHubbard"


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

def DMFT_SCC(params,t,fDelta,cmix = 0.5):
    """This subroutine creates Delta.inp from Gf.out for DMFT on bethe lattice: Delta=t^2*G
    If Gf.out does not exist, it creates Gf.out which corresponds to the non-interacting model
    In the latter case also creates the inpurity cix file, which contains information about
    the atomic states.
    """
    fileGf = 'Gf.out'
    if (os.path.exists(fileGf)): # If output file exists, start from previous iteration
        # Gf = io.read_array(fileGf, columns=(0,-1), lines=(1,-1))
        # In the new Python, io.readarray is dropped and we should use loadtxt instead!
        Gf = loadtxt(fileGf)
    else: # otherwise start from non-interacting limit
        print 'Starting from non-interacting model'
        Gf=[]
        for n in range(5000):
            iom = (2*n+1)*pi/params['beta'][0]
            gf = 0.5*t**-2*1j*(iom-sqrt(iom**2+4*t**2))
            Gf.append([iom, gf.real, gf.imag])
        Gf = array(Gf)
        
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
        
    # Preparing input file Delta.inp
    f = open(fDelta, 'w')
    if ((len(Gf.T)==len(old_Delta.T)) and Mixing == True):
        for i in xrange(len(Gf)):
            print >> f, Gf[i,0], (1-cmix)*t**2*Gf[i,1]+cmix*old_Delta[i,1], (1-cmix)*t**2*Gf[i,2]+cmix*old_Delta[i,2] # This is DMFT SCC: Delta = t**2*G (with t=1/2)
    else:
        for i in xrange(len(Gf)):
            print >> f, Gf[i,0], t**2*Gf[i,1], t**2*Gf[i,2] # This is DMFT SCC: Delta = t**2*G (with t=1/2)
    f.close()
    

#Doped Hubbard Model Scan mu from insulator
beta = 1/T
print "MC steps = %d"%int(ncores*Msteps)
#prefix = os.path.join("results",spath,"T.%.04f"%T,"mu.%.04f"%mu)
prefix = os.path.join("results",args.spath)
print "prefix = ",prefix
CreateFolder(prefix)
t=0.5

params = {"exe":   ["mpirun -np %d ctqmc"%ncores,          "# Path to executable"],
        "U":     [U,                 "# Coulomb repulsion (F0)"],
        "mu":    [mu,              "# Chemical potential"],
        "beta":  [beta,                "# Inverse temperature"],
        "M" :    [Msteps,                "# Number of Monte Carlo steps"],
        "mode":  ["SM",               "# S stands for self-energy sampling, M stands for high frequency moment tail"],
        "cix":   ["one_band.imp",     "# Input file with atomic state"],
        "Delta": ["Delta.inp",        "# Input bath function hybridization"],
        "tsample":[200,               "# how often to record the measurements" ],
        "nom":   [nom,                 "# number of Matsubara frequency points to sample"],
      "svd_lmax":[30,                 "# number of SVD functions to project the solution"],
        "aom":   [1,                  "# number of frequency points to determin high frequency tail"],
        "GlobalFlip":[1000000,         "# how often to perform global flip"],
        }

if (fromIC==True): #Copy insulator initial condition
    if (os.path.exists("initialcondition/Gf.out")):
        cmd = 'cp initialcondition/Gf.out run/Gf.out'
        subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying insulator Gf as initial condition
        print "mode: fromIC - Copying from initial condition"

# refine dmft run (need refine=True)
prev_it = 0
existpath = os.path.join(prefix,"Gf.out")
if (refine and os.path.exists(existpath)):
    print "mode: refine - Continue DMFT from previously obtained data"
    cpfrompath = os.path.join(prefix,"*")
    cmd = 'cp %s run/'%cpfrompath
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) 

    fGfs = glob.glob(os.path.join(prefix,"Gf.out.*"))

    if (len(fGfs)>0):
        sfGfs= sorted(fGfs,key=lambda f: int(filter(str.isdigit, f)))[-1]
        prev_it = int(re.findall(r"\d+", sfGfs)[-1])+1
#clear existing data if not refining
if (not(refine) and os.path.exists(existpath)):
    cpfrompath = os.path.join(prefix,"*")
    cmd = 'rm %s'%cpfrompath
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) 

# Creating parameters file PARAMS for qmc execution
CreateInputFile(params,"run")

cwd=os.getcwd()
os.chdir("run")

G0_list = []
slope_list = []

for it in range(Niter):
    
    import time
    t0=time.time()

    # Constructing bath Delta.inp from Green's function
    DMFT_SCC(params,t,params['Delta'][0],cmix=cmix)

    # Running ctqmc
    print 'Running ---- qmc itt.: ', it+prev_it, '-----'

    subprocess.call(params['exe'][0], shell=True,stdout=sys.stdout,stderr=sys.stderr)

    # Some copying to store data obtained so far (at each iteration)
    cmd = 'cp Gf.out Gf.out.'+str(it+prev_it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying Gf

    cmd = 'cp Sig.out Sig.out.'+str(it+prev_it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying Sig

    cmd = 'cp ctqmc.log ctqmc.log.'+str(it+prev_it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying log file

    #if it+prev_it>1:
        #diff = Diff('Gf.out', 'Gf.out.'+str(it+prev_it-1))
        #print 'Diff=', diff
        #diff0 = Diff0('Gf.out', 'Gf.out.'+str(it+prev_it-1))
        #print 'log(Diff0)=', log(abs(diff0))
    #if it+prev_it>2:
        #rrlambda = RaileighRitz('Gf.out.'+str(it+prev_it-2),'Gf.out.'+str(it+prev_it-1),'Gf.out')
        #savetxt("lambda.out",[rrlambda])
        #cmd = 'cp lambda.out lambda.out.'+str(it+prev_it)
        #subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying Gf
        #print 'lambda(RaileighRitz)=',rrlambda

    #Get Gf
    om,reGf,imGf = loadtxt("Gf.out").T
    G0 = reGf[0]+imGf[0]
    G0_list.append(G0)

    if it>6:
        iterrs, avg_G0abs, std_G0abs = paverage(abs(array(G0_list)),npoints=5,order=1)
        G0_diff_it,G0_diff = pdiff(avg_G0abs)
        G0_diff_avg_it,G0_diff_avg,G0_diff_std = paverage(abs(G0_diff),npoints=5,order=2)
        try:
            print "it %d : Diff = %le"%(it , G0_diff_avg[-1]/abs(G0_list[-1]))
        except:
            print ""
    slope = CalcSlope("Sig.out")
    slope_list.append(slope)
    if slope > 0:
        print "Insulator found"
    else:
        Z = 1/(1-slope)
        print "Metal found, Z = %f"%Z
        
    print "Time elapsed : %fs"%(time.time()-t0)


os.chdir(cwd)

if (adiabatic):
    print "adiabatic - saving results to initialcondition/Gf.out"
    cmd = 'cp run/Gf.out initialcondition/Gf.out'
    os.system(cmd)

cmd = 'cp run/* %s'%prefix
os.system(cmd)

cmd = 'rm run/*'
os.system(cmd)
