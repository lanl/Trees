import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os

def readfield(fname, Nx, Ny):
    infile = open(fname, 'rb')
    return np.fromstring(infile.read(Nx*Ny*4), 'f').reshape((Nx,Ny),order='F')

def plotTopo(z, z1, xlabel, N, d, title, fname):
    fig,axs = plt.subplots(1,2)
    img = axs[0].imshow(z, origin='lower', extent=[0,Nx*dx, 0, Ny*dy])
    plt.colorbar(img, ax=axs[0])
    axs[0].set_xlabel('x'+' [m]')
    axs[0].set_ylabel('y'+' [m]')

    axs[1].plot(np.arange(0,N*d,d),z1)
    axs[1].set_xlabel(xlabel+' [m]')
    axs[1].set_ylabel('z'+' [m]')

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(fname)
    return 0

def createFile(fname, z):
    if os.path.isfile(fname):
        os.remove(fname)
    file = FortranFile(fname, 'w')
    file.write_record(z)
    file.close()
    return 0

#params
Nx = 150
Ny = 100
dx = 2
dy = 2
of = './'
topofile = ''
if os.path.isfile(topofile):
    topo = readfield(topofile, Nx, Ny)
else:
    topo = np.zeros((Nx,Ny))

for i in np.arange(Nx):
    topo[i,:] = 0.1*i#25*( (1/(1+(np.exp( (-i+50)/10 )) )) )-0

z = np.transpose(topo)
z = z.astype('float32')

plotTopo(topo, topo[:,int(Ny/2)], 'x', Nx, dy, 'Topo', of+'Topo_plot.png')  

createFile(of+'topotest.dat',z)
#'''
