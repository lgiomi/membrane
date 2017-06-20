import numpy as np,sys

if len(sys.argv)<2:
    sys.exit("Syntax: python array_plot_3D.py folder_path")

import glob,os

if os.path.exists(sys.argv[1]):
	print "Opening folder "+sys.argv[1]
else:
	sys.exit("Error: input directory does not exist")

if not os.path.exists("img_temp"):
    os.makedirs("img_temp")

import matplotlib.pyplot as plt
from mayavi import mlab

listconf = sorted(glob.glob(sys.argv[1]+"/*last.dat"))

listk = list(set([i.split("_")[-5] for i in listconf]))
listk = sorted([float(i) for i in listk])
listkb = list(set([i.split("_")[-4] for i in listconf]))
listkb = sorted([float(i) for i in listkb])
listc = list(set([i.split("_")[-3] for i in listconf]))
listc = sorted([float(i) for i in listc])
listn = list(set([i.split("_")[-2] for i in listconf]))
listn = sorted([int(i) for i in listn])

coord3D = np.transpose(np.loadtxt("../geometry.dat",usecols=(1,2,3),delimiter='\t'))
coord2D = np.transpose(np.loadtxt("../geometry.dat",usecols=(7,8),delimiter='\t'))

triangles = np.loadtxt("../triangles.dat")

rgmap=[]

for i in range(256):
#	rgmap.append([255-i,i,0,255])
	rgmap.append([i,255-i,44*(1-i/255)+i*128/255,255])

f=listconf[0]
phi = np.transpose(np.loadtxt(f))
mlab.figure(bgcolor=(0,0,0))
s = mlab.triangular_mesh(coord3D[0],coord3D[1],coord3D[2],triangles,scalars = (phi[0]+1)/2,vmax=1,vmin=0,colormap='PiYG')
#s.module_manager.scalar_lut_manager.lut.table = rgmap
mlab.view(180, 0, 1.5)
n=f.split("_")
mlab.savefig("img_temp/t_"+str(listk.index(float(n[-5])))+"_"+str(listkb.index(float(n[-4])))+"_"+str(listc.index(float(n[-3])))+"_"+str(listn.index(int(n[-2])))+".png",size=(1920, 1080))

counter=1

s.scene.anti_aliasing_frames = 0
s.scene.disable_render = True

print ""

for f in listconf:

	n=f.split("_")
	print "\033[F"+"                                                                                   				 "
	print "\033[F"+"Processing file t_"+n[-5]+"_"+n[-4]+"_"+n[-3]+"_"+n[-2]+" ("+str(counter+1)+" of "+str(len(listconf))+")"
	counter=counter+1
	if os.path.isfile("img_temp/"+"t_"+str(listk.index(float(n[-5])))+"_"+str(listkb.index(float(n[-4])))+"_"+str(listc.index(float(n[-3])))+"_"+str(listn.index(int(n[-2])))+".png"):
		continue

	phi = np.transpose(np.loadtxt(f))
	s.mlab_source.scalars = (phi[0]+1)/2
	mlab.savefig("img_temp/"+"t_"+str(listk.index(float(n[-5])))+"_"+str(listkb.index(float(n[-4])))+"_"+str(listc.index(float(n[-3])))+"_"+str(listn.index(int(n[-2])))+".png",size=(1920, 1080))

print "\nThere are (",len(listk),",",len(listkb),") couplings,",len(listc),"concentrations and",len(listn),"series"
