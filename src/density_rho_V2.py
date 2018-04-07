#!/usr/bin/env python

import numpy as np
import os

import sys
sys.path.insert(0,os.getcwd())
from simp3D_x import *
from simp3D_y import *
from simp3D_z import *
#from simp3D import *

class CHD():
    def __init__(self):
        #This simply allocates the different data structures we need:
        self.natoms = 0 # to record number of atoms
        self.grid = np.zeros([0,0,0]) # 3D array for storing grid
        self.v = np.zeros([3,3]) # to record volume info
        self.N = np.zeros([3]) # to record amount of points in each direction
        self.position = np.zeros([0,0]) # to record atom type and positions
        self.origin = np.zeros([3]) # to record the origin of the density

    def integrate(self):
        #This allows us to integrate the stored charge density
        hx = self.v[0,0]
        hy = self.v[1,1]
        hz = self.v[2,2]
        NX = int(self.N[0])
        NY = int(self.N[1])
        NZ = int(self.N[2])
        U = self.grid
        #print("hx: " + str(hx))
        #print("hy: " + str(hy))
        #print("hz: " + str(hz))
        #print("NX: " + str(NX))
        #print("NY: " + str(NY))
        #print("NZ: " + str(NZ))
        #total = simp3D(U,hx,hy,hz,NX,NY,NZ)
        (xVec,out) = simp3D_x(U,hx,hy,hz,NX,NY,NZ)
        (yVec,out2) = simp3D_y(U,hx,hy,hz,NX,NY,NZ)
        (zVec,out3) = simp3D_z(U,hx,hy,hz,NX,NY,NZ)
        #return(xVec,out,total)
        return(xVec,yVec,zVec,out)


#The following function reads the charge density from a cube file
#TODO Make gaussian cube file auto detection
def read(cubefile,saveDest,gaussianLogical):
    density = CHD()
    f = open(cubefile,'r')

    #skip two header lines
    next(f)
    next(f)

    line = next(f).split()
    #Get the number of atoms if we want to store it
    density.natoms = int(line[0])
    #Get the origin of the density
    for i in range(0,3):
        density.origin[i] = float(line[i+1])

    #This gets the nx,ny,nz info of the charge density
    #As well as the differential volume
    for i in range(0,3):
        line = next(f).split()
        density.N[i] = int(line[0]) #FIXME Despite calling int() number is being stored as float
        for j in range(1,4):
            density.v[i][j-1] = float(line[j])

    density.position = np.zeros([density.natoms,4])

    #Get the type and positions of the atoms
    for i in range(0,density.natoms):
        line = next(f).split()
        for j in range(1,5):
            density.position[i][j-1] = float(line[j])

    # Initializes the size of the density.grid vector
    density.grid = np.zeros([int(density.N[0]),int(density.N[1]),int(density.N[2])])


    zLines = int(np.ceil(density.N[2]/6.0))
    for i in range(int(density.N[0])):
        for j in range(int(density.N[1])):
            zcount=0
            for k in range(zLines):
                line = next(f).split()
                for l in line:
                    density.grid[i,j,zcount] = float(l)
                    zcount += 1
    f.close()

    return density


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments, run as ./density.py CUBEFILELOCATION SAVEDESTINATION gaussianLogical")
        print("gaussianLogical: Specify if cube file is gaussian cube file (TRUE or FALSE)")
        sys.exit(6)

    density = read(sys.argv[1],sys.argv[2],sys.argv[3])

    xOrig = density.origin[0]
    xStepSize = density.v[0,0]

    yOrig = density.origin[1]
    yStepSize = density.v[1,1]

    zOrig = density.origin[2]
    zStepSize = density.v[2,2]

    #print("Density using Riemann Sum: " + str(np.sum(density.grid)*xStepSize*yStepSize*zStepSize))
    (xVec,yVec,zVec,out) = density.integrate()
    print("Total density using Simpson's Rule New Code: "+str(out))

    cubeInput = sys.argv[1].split('/')
    cubeFile = cubeInput[-1]
    filename = cubeFile.split('.cube')[0] + '_rho_x.txt'
    filename2 = cubeFile.split('.cube')[0] + '_rho_y.txt'
    filename3 = cubeFile.split('.cube')[0] + '_rho_z.txt'

    saveDest = sys.argv[2]

    with open(os.path.join(saveDest,filename),'w') as file1:
        file1.write(cubeFile+'\n')
        file1.write("Total density using Simpson's Rule: " + str(out) + '\n')
        file1.write('{0:6s}  {1:>10s}  {2:>12s}\n'.format('X Step','X Coord','Rho'))
        for i in range(len(xVec)):
            xCoord = float(xOrig) + float(xStepSize)*float(i)
            file1.write('{0:6d}  {1:10.4f}  {2:12.6f}\n'.format(i,xCoord,xVec[i]))

    with open(os.path.join(saveDest,filename2),'w') as file1:
        file1.write(cubeFile+'\n')
        file1.write("Total density using Simpson's Rule: " + str(out) + '\n')
        file1.write('{0:6s}  {1:>10s}  {2:>12s}\n'.format('Y Step','Y Coord','Rho'))
        for i in range(len(yVec)):
            yCoord = float(yOrig) + float(yStepSize)*float(i)
            file1.write('{0:6d}  {1:10.4f}  {2:12.6f}\n'.format(i,yCoord,yVec[i]))

    with open(os.path.join(saveDest,filename3),'w') as file1:
        file1.write(cubeFile+'\n')
        file1.write("Total density using Simpson's Rule: " + str(out) + '\n')
        file1.write('{0:6s}  {1:>10s}  {2:>12s}\n'.format('Z Step','Z Coord','Rho'))
        for i in range(len(zVec)):
            zCoord = float(zOrig) + float(zStepSize)*float(i)
            file1.write('{0:6d}  {1:10.4f}  {2:12.6f}\n'.format(i,zCoord,zVec[i]))


