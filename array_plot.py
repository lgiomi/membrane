import numpy as np,sys,os
import matplotlib.pyplot as plt
import glob

if len(sys.argv)<1:
    sys.exit("Syntax: python array_plot.py folder")

if os.path.exists(sys.argv[1]):
	print "Opening folder "+sys.argv[1]+""
else:
	sys.exit("Error: input directory does not exist")

if not os.path.exists("img_temp"):
    os.makedirs("img_temp")

from mayavi import mlab
mlab.options.offscreen = True

listconf = sorted(glob.glob(sys.argv[1]+"/conf*.dat"))

triangles = np.loadtxt("triangles.dat")

coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
curvatures = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))

mlab.figure(bgcolor=(1,1,1))
s = mlab.triangular_mesh(coord3D[0],coord3D[1],coord3D[2],triangles,scalars = curvatures[0],vmax=1,vmin=0,colormap='PiYG')

mlab.view(0, 180,5)

if os.path.isfile("img_temp/H2.png")==False:
	mlab.savefig("img_temp/H2.png",size=(1920, 1080))

if os.path.isfile("img_temp/KG.png")==False:
	s.mlab_source.scalars = curvatures[1]
	mlab.savefig("img_temp/KG.png",size=(1920, 1080))

counter=0

for f in listconf:

	n=f.split("_")
	namefile="img_temp/conf_"+n[1]+"_"+n[2]+"_"+n[3]+"_"+n[4]+"_"+n[5]+"_"+n[6]+"_"+n[7]+"_"+n[8]+".png"

	counter=counter+1

	if os.path.isfile(namefile):
		continue

	phi = np.transpose(np.loadtxt(f))

	print "\033[F"+"Processing file "+str(counter)+" of "+str(len(listconf))

	s.mlab_source.scalars = (phi[0]+1)/2
	mlab.savefig(namefile,size=(1920, 1080))

mlab.close()
