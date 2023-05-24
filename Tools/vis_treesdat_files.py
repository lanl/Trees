#============================================================================
#   This script generates two figures of .dat files obtained from trees program.
#   Figure 1 is an X/Y plot of fuels sliced on z=plane
#   Figure 2 is an X/Y plot of fuels where species are summed along z
#   Input:  .dat output from trees, number of fuel types
#   Output: X/Y figures of .dat output from trees    
#   Original file:       Julia Oliveto, joliveto@lanl.gov, April 5th, 2022
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
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os.path
import sys
import struct

#====================INPUTS====================
#pathfiles------------------
#pf = pathfile where .dat file lives
pf = './'
#of = out pathfile where you save image (as .png)
of = './'
#b = the name of the fuels case you are visualizing
b = 'test_case'
#topofile = pathfile and name of topography (if used), if no topography input: ''
topofile = ''#'/Users/joliveto/Desktop/trees_fix_topo/Inputs/rampHill.dat'

#trees/viewing parameters------------------
datfile  = 'treesrhof.dat' #which .dat to make png
nfuel    = 5 #see bottom of trees output, number of output fuels
sum_flag = 1 #if = 0, the topdown figures will be z-slices specified by "plane" ; if = 1 plots topdown figures as a sum in z
plane    = 0 #z-index slice of plotting (0=ground/bottom layer), will *not* be used if sum_flag = 1

#grid parameters------------------
Nx      = 350 
Ny      = 250 
Nz      = 41 
dx      = 2.0
dy      = 2.0
dz      = 15.0
aa1     = 0.1
f1      = 0.0
stretch = 2 

#======================= DEFINE FUNCTIONS =========================                                                                                            

def readfield(fuelfile, Nx, Ny, Nz):
    np.frombuffer(fuelfile.read(4),'f')
    return np.frombuffer(fuelfile.read(Nx*Ny*Nz*4), 'f').reshape((Nx,Ny,Nz),order='F')

def zheight(ZI):
    # generates array of cell heights from z-index array                                                                                
    Z = np.copy(ZI, order='K')
    ZItemp = Z[0,0,:]
    ZItemp[0] = ZItemp[0] * 2
    for ii in range(1,len(ZItemp)):
        ZItemp[ii] = (ZItemp[ii] - sum(ZItemp[:ii]))*2
    for ii in range(len(ZItemp)):
        Z[:,:,ii] = ZItemp[ii]
    return Z

def metrics(topofile, Nx, Ny, Nz, dx, dy, dz, a1, f0, Stretch):
    # --- read topo file if present ---
    if os.path.isfile(topofile):
        #topo = numpy.zeros((Nx,Ny))
        f = open(topofile, 'rb')
        pre1 = struct.unpack("i",f.read(4))[0]
        pre2 = struct.unpack("i",f.read(4))[0]
        topo=np.fromstring(f.read(Nx*Ny*4),'f').reshape((Nx,Ny), order = 'F')
        f.close()

    # --- build base grid ---
    x = np.zeros((Nx))
    y = np.zeros((Ny))
    z = np.zeros((Nz))
    zedge = np.zeros((Nz+1))
    XI = np.zeros((Nx,Ny,Nz))
    YI = np.zeros((Nx,Ny,Nz))
    ZI = np.zeros((Nx,Ny,Nz))
    for i in range(Nx):
        x[i] = i*dx - 0.5*Nx*dx
    for j in range(Ny):
        y[j] = j*dy - 0.5*Ny*dy
    for k in range(Nz):
        z[k] = k*dz + 0.5*dz
        zedge[k] = k*dz
    zedge[Nz] = Nz*dz
    # --- using no stretching ---
    if Stretch == 0:
        print('not using stretching')
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    XI[i,j,k] = x[i]
                    YI[i,j,k] = y[j]
                    ZI[i,j,k] = z[k]    
    # --- using hyperbolic tangent stretching --- 
    if Stretch == 1:
        print('using hyperbolic tangent stretching')
        print('this part does not work yet! exiting!')
        sys.exit()
    # --- using cubic polynomial stretching ---
    if Stretch == 2:
        print('using cubic polynomial stretching')
    # --- set cubic polynomial 2nd and 3rd term coefficients ---
        a2 = f0*(1.0-a1)/zedge[Nz]
        a3 = (1.0-a2*zedge[Nz]-a1)/(zedge[Nz]**2.0)
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    XI[i,j,k] = x[i]
                    YI[i,j,k] = y[j]
                    ZI[i,j,k] = (a3*(z[k]**3.0)+a2*(z[k]**2.0)+a1*z[k])*(zedge[Nz]-zedge[0])/zedge[Nz]+zedge[0]
        if os.path.isfile(topofile):
            print("Modifying coordinate to be terrain following!")
            for i in range(Nx):
                for j in range(Ny):
                    for k in range(Nz):
                        ZI[i,j,k] = ZI[i,j,k]*(zedge[Nz]-topo[i,j])/zedge[Nz] + topo[i,j]
    Z = zheight(ZI)
    volume = np.multiply(dx,dy,Z)
    return XI, YI, ZI, volume

def plotTopdown(fig,axs,arr,title,X,Y,sum_flag,plane):
    if sum_flag == 1:
        arr = np.sum(arr,axis=2)
        sp1 = axs.pcolormesh(X[:,:,0],Y[:,:,0],arr,cmap='Greens',shading='auto', )
        #sp1 = axs.imshow(arr,cmap='Greens')
    else:
        arr = arr[:,:,plane]
        sp1 = axs.pcolormesh(X[:,:,0],Y[:,:,0],arr,cmap='Greens',shading='auto')
        #sp1 = axs.imshow(arr,cmap='Greens')
    cbar = fig.colorbar(sp1, ax=axs)
    #cbar.ax.set_ylabel(ylabel, rotation=270)
    axs.set_title(title)
    axs.set_xlabel('X [m]')
    axs.set_ylabel('Y [m]')
    return 0

def plotVertical(fig,axs,arr,title,X,Z):
    sp1 = axs.pcolormesh(X[:,int(Ny/2),:],Z[:,int(Ny/2),:],arr[:,int(Ny/2),:],cmap='Greens',shading='auto')
    cbar = fig.colorbar(sp1, ax=axs)
    #cbar.ax.set_ylabel(ylabel, rotation=270)
    axs.set_title(title)
    axs.set_xlabel('X [m]')
    axs.set_ylabel('Z [m]')
    return 0


    
#STEP 1======================
#open rhof
#rhof is a 4D array: 1st position is species ID, 2nd is x, 3rd is y, 4th is z
rhof= np.zeros(nfuel*Nx*Ny*Nz).reshape(nfuel,Nx,Ny,Nz)
rhoffile = open(pf+datfile,'rb')
for ift in range(nfuel):
    print('Reading fuel type:',ift)
    rhof[ift,:,:,:] = readfield(rhoffile,Nx,Ny,Nz)
    trhof = rhof[ift,:,:,:]
    print( 'SPECIES ',ift+1,' MIN = ',np.min(trhof) ,' ; MAX = ',np.max(trhof))
rhoffile.close()
print(rhof.shape)

#STEP 2======================
#Create X,Y grid based on Nx,Ny,dx,dy
XI, YI, ZI, vol = metrics(topofile, Nx, Ny, Nz, dx, dy, dz, aa1, f1, stretch)
X = XI[:,0,0]
Y = YI[0,:,0]
Z = ZI[0,0,:]

#STEP 3======================
#Visualilize the X/Y plane of the .dat file
name = datfile[:-4] #remove ".dat" from string
#find correct number of subplots for nfuel
rows = int(np.floor(np.sqrt(nfuel+1)))
cols = int(rows + np.ceil(np.sqrt(nfuel+1) - rows))
if (rows<(nfuel+1)/cols):
    rows += 1
print(rows, cols)

fig,axs = plt.subplots(rows,cols, figsize=(15,10))
for n in range(nfuel+1):
    if (n==0):
        arr = np.sum(rhof, axis=0)
        t = 'Sum Species'
    else:
        arr = rhof[n-1,:,:,:]
        t = 'Species '+str(n)
    if (rows!=1):
        inpaxs = axs[int(np.floor((n)/cols)),((n)%cols)]
    else:
        inpaxs = axs[n]
    plotTopdown(fig,inpaxs,arr,t,XI+Nx,YI+Ny,sum_flag,plane)  
#Hide unused subplots
for nn in range(nfuel+1,rows*cols):
    axs[int(np.floor((nn)/cols)),((nn)%cols)].axis('off')
#figure title    
fig.suptitle(name)
plt.tight_layout()
plt.savefig(of+name+'_'+b+'_topdownview.png')
plt.close()    
#=================================

fig,axs = plt.subplots(rows,cols, figsize=(15,10))

for n in range(nfuel+1):
    if (n==0):
        arr = np.sum(rhof, axis=0)
        t = 'Sum Species'
        max_z = int(np.max(np.nonzero(np.sum(arr,axis=(0,1)))))+2
    else:
        arr = rhof[n-1,:,:,:]
        t = 'Species '+str(n)
    if (rows!=1):
        inpaxs = axs[int(np.floor((n)/cols)),((n)%cols)]
    else:
        inpaxs = axs[n]
    plotVertical(fig,inpaxs,arr[:,:,:max_z],t,XI[:,:,:max_z]+Nx,ZI[:,:,:max_z])  
#Hide unused subplots
for nn in range(nfuel+1,rows*cols):
    axs[int(np.floor((nn)/cols)),((nn)%cols)].axis('off')
#figure title    
fig.suptitle(name)
plt.tight_layout()
plt.savefig(of+name+'_'+b+'_verticalview.png')
plt.close()    
#=================================
