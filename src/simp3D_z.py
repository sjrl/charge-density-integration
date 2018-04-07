# Created by Sebastian Lee

import numpy as np

# U --> 3D tensor(grid) containing the value of the density at each x,y,z point
# hx --> step size in x direction
# hy --> step size in y direction
# hz --> step size in z direction
# NX --> number of x points (need to have an odd amount)
# NY --> number of y points (need to have an odd amount)
# NZ --> number of z points (need to have an odd amount)

# Note: Number of intervals NX,NY,NZ must be even
def simp3D_z(U,hx,hy,hz,NX,NY,NZ):

    # Initialize the z vector
    zVec = np.zeros([NZ])

    xSimpWeight = np.zeros([NX])
    xSimpWeight[1::2] = 4
    xSimpWeight[::2] = 2
    xSimpWeight[0] = 1
    xSimpWeight[-1] = 1

    ySimpWeight = np.zeros([NY])
    ySimpWeight[1::2] = 4
    ySimpWeight[::2] = 2
    ySimpWeight[0] = 1
    ySimpWeight[-1] = 1

    zSimpWeight = np.zeros([NZ])
    zSimpWeight[1::2] = 4
    zSimpWeight[::2] = 2
    zSimpWeight[0] = 1
    zSimpWeight[-1] = 1

    for i in range(NZ):
        xVec = np.dot(U[:,:,i],ySimpWeight)
        zVec[i] = np.dot(xSimpWeight,xVec)

    # Finally add all the contributions and multiply by the step sizes hx, hy and hz
    # and a factor 1/27 (1/3 in each direction).

    zVec = (zVec*hx*hy)/9.0

    out = (np.dot(zSimpWeight,zVec)*hz)/3.0

    return(zVec,out)
