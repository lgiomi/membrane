import numpy as np,sys
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

if len(sys.argv)<3:
	sys.exit("Syntax: python plot_3D.py input.dat output.png")

try: 
	coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
	coord2D = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))
	triangles = np.loadtxt("triangles.dat")
except:
	sys.exit("Error: missing geometry files")

if len(sys.argv)==4:
	try:
		ofs=float(sys.argv[3])
	except:
		sys.exit("Error: not a valid offset")
	if ofs <= 0:
		sys.exit("Error: offset has to be a positive number")
	for p in range(len(np.transpose(coord3D))):
		if coord3D[0,p]<0:
			coord3D[0,p]+=ofs
		else:
			coord3D[0,p]-=ofs
try:
	phi = np.transpose(np.loadtxt(sys.argv[1]))
except:
	sys.exit("Error: input file does not exist")

rgmap=[]

for i in range(256):
#	rgmap.append([i,255-i,44*(1-i/255)+i*128/255,255])
	rgmap.append([255-i,i,0,255])

from mayavi import mlab
mlab.figure(bgcolor=(1,1,1))
s = mlab.triangular_mesh(coord3D[0],coord3D[1],coord3D[2],triangles,scalars = (phi[0]+1)/2,vmax=1,vmin=0,colormap='PiYG', mode='rgba')
#s.module_manager.scalar_lut_manager.lut.table = rgmap
#mlab.view(90, 0)
mlab.view(45, 45, 1.5)
mlab.savefig(sys.argv[2],size=(1920, 1080))

#mlab.show()

