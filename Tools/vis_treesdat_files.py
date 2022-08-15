#============================================================================
#   This script generates two figures of .dat files obtained from trees program.
#   Figure 1 is an X/Y plot of fuels sliced on z=plane
#   Figure 2 is an X/Y plot of fuels where species are summed along z
#   Input:  .dat output from trees, number of fuel types
#   Output: X/Y figures of .dat output from trees    
#   Original file:       Julia Oliveto, joliveto@lanl.gov, April 5th, 2022
#============================================================================

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os.path

#patfiles------------------
#pf = pathfile where .dat file lives
pf = '/Users/joliveto/Desktop/Projects/Rod_Projects/JG_fuel_hetero/Flagstaff_newtrees_fuels/Inputs_Flagstaff/'
#of = outfile where you save image (as .png)
of = '/Users/joliveto/Desktop/Projects/Rod_Projects/JG_fuel_hetero/fuels_pngs/'
#b = the name of the fuels case you are visualizing
b = 'test_fuels'

#trees/viewing parameters------------------
datfile = 'treesrhof.dat' #which .dat to make png
nfuel = 3 #see bottom of trees output, number of output fuels
plane = 0 #z-index slice of plotting (0=ground/bottom layer)

#grid parameters------------------
Nx    = 300 
Ny    = 300 
Nz    = 41 
dx    = 2.0
dy    = 2.0


#======================= DEFINE FUNCTIONS =========================                                                                                            

def readfield(fuelfile, Nx, Ny, Nz):
    np.frombuffer(fuelfile.read(4),'f')
    return np.frombuffer(fuelfile.read(Nx*Ny*Nz*4), 'f').reshape((Nx,Ny,Nz),order='F')

def createGrid(Nx, Ny, dx, dy):
    X = np.zeros(Nx+1)
    for i in range(Nx+1):
        X[i] = dx*i
    Y = np.zeros(Ny+1)
    for i in range(Ny+1):
        Y[i] = dy*i
    return X,Y

def plotTopdown(fig,axs,arr,title,X,Y,plane):
    if plane!=9999:
        sp1 = axs.pcolormesh(X,Y,np.transpose(arr[:,:,plane]),cmap='viridis',shading='auto', vmin=0)
    else:  
        sp1 = axs.pcolormesh(X,Y,np.transpose(arr[:,:]),cmap='viridis',shading='auto', vmin=0)
    cbar = fig.colorbar(sp1, ax=axs)
    #cbar.ax.set_ylabel(ylabel, rotation=270)
    axs.set_title(title)

def plotTopdownSum(fig,axs,arr,title,X,Y):
    sp1 = axs.pcolormesh(X,Y,np.transpose(np.sum(arr, axis=2)),cmap='viridis',shading='auto', vmin=0)
    cbar = fig.colorbar(sp1, ax=axs)
    #cbar.ax.set_ylabel(ylabel, rotation=270)
    axs.set_title(title)


    
#STEP 1======================
#open rhof
#rhof is a 4D array: 1st position is species ID, 2nd is x, 3rd is y, 4th is z
rhof= np.zeros(nfuel*Nx*Ny*Nz).reshape(nfuel,Nx,Ny,Nz)
rhoffile = open(pf+datfile,'rb')
for ift in range(nfuel):
    print('Reading fuel type:',ift)
    rhof[ift,:,:,:] = readfield(rhoffile,Nx,Ny,Nz)
rhoffile.close()
print(rhof.shape)

#STEP 2======================
#Create X,Y grid based on Nx,Ny,dx,dy
X,Y = createGrid(Nx, Ny, dx, dy)

#STEP 3======================
#Visualilize the X/Y plane of the .dat file
name = datfile[:-4] #remove ".dat" from string

#find correct number of subplots for nfuel
rows = int(np.floor(np.sqrt(nfuel+1)))
cols = int(rows + np.ceil(np.sqrt(nfuel+1) - rows))
if (rows<(nfuel+1)/cols):
    rows += 1
print(rows, cols)

#fig 1: single zplane .dat file=========================== 
fig,axs = plt.subplots(rows,cols, figsize=(15,10))
print('FIG1')
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
    plotTopdown(fig,inpaxs,arr,t,X,Y,plane)  
#Hide unused subplots
for nn in range(nfuel+1,rows*cols):
    axs[int(np.floor((nn)/cols)),((nn)%cols)].axis('off')
#figure title    
fig.suptitle(name)
plt.tight_layout()
plt.savefig(of+name+'_'+b+'_topdown_view_zplane_'+str(plane)+'.png')
plt.close()   
#=================================== 

#fig 2: sum .dat file===========================
fig,axs = plt.subplots(rows,cols, figsize=(15,10))
print('FIG2')

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
    plotTopdownSum(fig,inpaxs,arr,t,X,Y)  
#Hide unused subplots
for nn in range(nfuel+1,rows*cols):
    axs[int(np.floor((nn)/cols)),((nn)%cols)].axis('off')
#figure title    
fig.suptitle(name)
plt.tight_layout()
plt.savefig(of+name+'_'+b+'_topdown_view_sumZ.png')
plt.close()    
#=================================
