import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os.path

#patfiles------------------
#pf = pathfile where .dat file lives
pf = '/Users/joliveto/Desktop/Projects/CreepyFire/Blodgett/blodgett1/' 
#of = outfile where you save image (as .png)
of = '/Users/joliveto/Desktop/Projects/CreepyFire/Blodgett/blodgett1/'

#trees/viewing parameters------------------
datfile = 'treesrhof.dat' #which .dat to make png
nfuel = 6 #see bottom of trees output, number of output fuels
plane = 0 #z-index slice of plotting (0=ground/bottom layer)

#grid parameters------------------
Nx    = 250 
Ny    = 250 
Nz    = 41  
dx    = 2.0
dy    = 2.0
dz    = 15.0
aa1   = 0.1


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
    sp1 = axs.pcolormesh(X,Y,np.transpose(arr[:,:,plane]),cmap='Greens',shading='auto')
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

#STEP 2======================
#Create X,Y grid based on Nx,Ny,dx,dy
X,Y = createGrid(Nx, Ny, dx, dy)

#STEP 3======================
#Visualilize the X/Y plane of the .dat file
name = datfile[:-4] #remove ".dat" from string

for n in range(nfuel+1):
    fig,axs = plt.subplots(figsize=(10,8))
    if (n==0):
        arr = np.sum(rhof, axis=0)
        t = 'Sum.Species'
    else:
        arr = rhof[n-1,:,:,:]
        t = 'Species.'+str(n)
    plotTopdown(fig,axs,arr,t,X,Y,plane)
    fig.suptitle(name)
    plt.tight_layout()
    plt.savefig(of+name+'_'+t+'_topdown_view.png')
    plt.close()    
         
             
 
    

