import numpy as np,sys
import matplotlib.pyplot as plt
import glob
from mayavi import mlab

listconf = glob.glob("gc*.dat")

coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
coord2D = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))

triangles = np.loadtxt("triangles.dat")

rgmap=[]

for i in range(256):
	rgmap.append([255-i,i,0,255])

counter=0

for f in listconf:

	n=f.split("_")[1]
	n=n.split(".")[0]

	print "\033[F"+"                                                                                    "
	print "\033[F"+"Processing file "+f+" ("+str(counter+1)+" of "+str(len(listconf))+")"
    	sys.stdout.flush()

	phi = np.transpose(np.loadtxt(f))

	s = mlab.triangular_mesh(coord3D[0],coord3D[1],coord3D[2],triangles,scalars = (phi[0]+1)/2,vmax=1,vmin=0)
	s.module_manager.scalar_lut_manager.lut.table = rgmap
	mlab.view(90, 0)
	mlab.savefig("t_"+n+".png",size=(1920, 1080))

        mlab.close()

	counter=counter+1


