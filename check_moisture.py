import numpy as np


def readfield(infile, Nx, Ny, Nz):
  return np.frombuffer(infile.read(Nx*Ny*Nz*4), 'f').reshape((Nx,Ny,Nz),order='F')


def read_rhof(datfile : str, Nx : int, Ny : int, Nz : int, nfuel : int):
    rhof= np.zeros(nfuel*Nx*Ny*Nz).reshape(nfuel,Nx,Ny,Nz)
    rhoffile = open(datfile,'rb')
    for ift in range(nfuel):
        print('Reading fuel type:',ift)
        rhof[ift,:,:,:] = readfield(rhoffile,Nx,Ny,Nz)
        trhof = rhof[ift,:,:,:]
        print( 'SPECIES ',ift+1,' MIN = ',np.min(trhof) ,' ; MAX = ',np.max(trhof))
    rhoffile.close()
    return rhof

if __name__ == "__main__":
    moist_file  = "control/treesmoist.dat"
    moist = read_rhof(moist_file, 300,300,41,1)
    print(moist.shape)
    moist = moist[0,:,:,:]
    print(np.mean(moist[:, :, 0]))

