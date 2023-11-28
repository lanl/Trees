#============================================================================
#   This script does post-processing treatments for trees program output.
#   Input:  .dat output from trees, number of fuel types
#   Output: .dat trees, adjusted/treated values of density, moisture, and/or depth 
#   CURRENTLY THIS SCRIPT WILL REMOVE ALL CANOPY FUELS FROM THE LEFT HAND SIDE OF THE DOMAIN
#   Original file: Julia Oliveto, joliveto@lanl.gov, November 28th, 2023
#=======#============================================================================
# © 2022. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#============================================================================
import numpy as np
from scipy.io import FortranFile  
import os


#======================= USER INPUTS =========================            
#patfiles------------------
#pf = pathfile where .dat file lives
pf = './Inputs/'
#of = where the treated .dat files will be saved
of = pf

#trees/viewing parameters------------------
datfile1 = 'treesrhof.dat' #which .dat to make png
datfile2 = 'treesmoist.dat'
datfile3 = 'treesfueldepth.dat'
datfile4 = 'treesss.dat'
nfuel    = 4 #see bottom of trees output, number of output fuels

#grid parameters------------------
Nx    = 100 
Ny    = 100 
Nz    = 41 
dx    = 2.0
dy    = 2.0
dz    = 15.0
aa1   = 0.1

#======================= DEFINE FUNCTIONS =========================                                                                                            
def readfield(fuelfile, Nx, Ny, Nz):
    np.frombuffer(fuelfile.read(4),'f')
    return np.frombuffer(fuelfile.read(Nx*Ny*Nz*4), 'f').reshape((Nx,Ny,Nz),order='F')

def readFile(pf, datfile, nfuel, Nx, Ny, Nz):
    rhoffile = open(pf+datfile,'rb')
    print('Reading fuel: ',datfile)
    rhof = np.zeros(nfuel*Nx*Ny*Nz).reshape(nfuel,Nx,Ny,Nz)
    for ift in range(nfuel):
        rhof[ift,:,:,:] = readfield(rhoffile,Nx,Ny,Nz)
    rhoffile.close()
    print(rhof.shape)
    return rhof

def writeFile(rhof, of, datfile, nfuel):
    if os.path.isfile(of+datfile):
        os.remove(of+datfile)
    f = FortranFile(of+datfile, 'w') #open fortran file
    #f.write_record(rhof.T )
    for i in np.arange(0,nfuel):
        trhof = rhof[i,:,:,:].astype('float32')
        f.write_record(trhof.T )
    f.close()
    return 0

     
#STEP 1======================
#open fuel files and store in numpy array from processing
#These are 4D arrays: 1st position is species ID, 2nd is x, 3rd is y, 4th is z
rhof  = readFile(pf, datfile1, nfuel, Nx, Ny, Nz) #FUEL DENSITY ARRAY
moist = readFile(pf, datfile2, nfuel, Nx, Ny, Nz) #FUEL MOISTURE ARRAY
fd    = readFile(pf, datfile3, nfuel, Nx, Ny, Nz) #FUEL DEPTH ARRAY
ss    = readFile(pf, datfile4, nfuel, Nx, Ny, Nz) #FUEL SIZE SCALE ARRAY

#STEP 2=======================
#Here is where you can adjust fuels and create/apply treatments
#The structure of the fuel arrays is [nfuel, Nx, Ny, Nz]
#Note when manipulating arrays, put locations in terms of CELLS, NOT METERS
#The example below assumes we have 1 grass fuel type and nfuel-1 canopy fuel types
#We will remove canopy fuels from the left hand side of the domain:
if nfuel>1:
    for i in range(1,nfuel):
        rhof[i,:25,:,:] = 0.0
        moist[i,:25,:,:] = 0.0
        fd[i,:25,:,:] = 0.0

#STEP 3=======================
#Write the output files to your specificed location
writeFile(rhof , of, datfile1, nfuel)
writeFile(moist, of, datfile2, nfuel)
writeFile(fd   , of, datfile3, nfuel)
writeFile(ss   , of, datfile4, nfuel)