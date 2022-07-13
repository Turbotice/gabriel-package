import numpy as np

def NextPowerOfTwo(number):
    # Returns next power of two following 'number'
    return np.ceil(np.log2(number))

def Zeropadding1D(Z, A = 4, B = 2):
    [nt, ny]= Z.shape
    Nzero = int(A*pow(B,NextPowerOfTwo(ny)))
    Zeropadded1D = np.zeros((nt,ny+Nzero))
    Zeropadded1D[:,:ny] = Z
    return Zeropadded1D


def Zeropadding2D(Z,A=2,B=2):
    [nx, ny]= Z.shape
    Nzx = int(A*pow(2,NextPowerOfTwo(nx)))
    Nzy = int(A*pow(2,NextPowerOfTwo(ny)))

    #rajouter une condition sur np.mod
    nx0 = int((Nzx-nx)/2)
    ny0 = int((Nzy-ny)/2)
    
    Zeropadded2D = np.zeros((Nzx,Nzy))+0*1j
    Zeropadded2D[nx0:nx0+nx,ny0:ny0+ny] = Z
    return Zeropadded2D