import numpy as np,sys
import matplotlib.pyplot as plt
import glob
from mayavi import mlab

mlab.options.offscreen = True

listconf = sorted(glob.glob("gc*.dat"))

coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
coord2D = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))

triangles = np.loadtxt("triangles.dat")

rgmap=[]

for i in range(256):
#	rgmap.append([255-i,i,0,255])
#	rgmap.append([0,255-i,i,255])
#	rgmap.append([255,i,0,255])
	rgmap.append([i,255-i,44*(1-i/255)+i*128/255,255])

f=listconf[0]
phi = np.transpose(np.loadtxt(f))
mlab.figure(bgcolor=(1,1,1))
s = mlab.triangular_mesh(coord3D[0],coord3D[1],coord3D[2],triangles,scalars = (phi[0]+1)/2,vmax=1,vmin=0,colormap='PiYG')
#s.module_manager.scalar_lut_manager.lut.table = rgmap
#mlab.view(45, 45, 1.5)
mlab.view(0, 180,4)
#mlab.view(45, -45,2.1)
n=f.split("_")[1]
n=n.split(".")[0]
mlab.savefig("t_"+n+".png",size=(1920, 1080))

counter=1

s.scene.anti_aliasing_frames = 0
s.scene.disable_render = True

print ""

for f in listconf:

	n=f.split("_")[1]
	n=n.split(".")[0]

	print "\033[F"+"                                                                                 "
	print "\033[F"+"Processing file "+f+" ("+str(counter+1)+" of "+str(len(listconf))+")"

	phi = np.transpose(np.loadtxt(f))
	s.mlab_source.scalars = (phi[0]+1)/2
	#mlab.view(counter, 45,1.3)
	#mlab.view(counter, 135,2.5)
	mlab.savefig("t_"+n+".png",size=(1920, 1080))

	counter=counter+1

mlab.close()
