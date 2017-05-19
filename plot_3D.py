import numpy as np,sys
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

if len(sys.argv)<3:
    sys.exit("Syntax: python plot_3D.py input.dat output.png")

coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
coord2D = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))

triangles = np.loadtxt("triangles.dat")

try:
	phi = np.transpose(np.loadtxt(sys.argv[1]))
except:
    sys.exit("Error: input file does not exist")

rgmap=[]

for i in range(256):
	rgmap.append([255,i,0,255])

from mayavi import mlab
s = mlab.triangular_mesh(coord3D[0],coord3D[1],coord3D[2],triangles,scalars = (phi[0]+1)/2,vmax=1,vmin=0)
s.module_manager.scalar_lut_manager.lut.table = rgmap
mlab.view(90, 0)
mlab.savefig(sys.argv[2],size=(1920, 1080))

mlab.show()

