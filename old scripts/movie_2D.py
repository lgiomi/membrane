import numpy as np,sys
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import glob

listconf = sorted(glob.glob("gc*.dat"))

coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
coord2D = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))

triangles = np.loadtxt("triangles.dat")

cdict = {'red': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),'green':((0.0, 0.0, 0.0),(1.0, 1.0, 1.0)),'blue':((0.0, 0.0, 0.0),(1.0, 0.0, 0.0))}
red_to_green = LinearSegmentedColormap('BlueRed1', cdict)

f=listconf[0]
phi = np.transpose(np.loadtxt(f))
fig, ax = plt.subplots()
ax.scatter(coord2D[0], coord2D[1], c=(phi[0]+1)/2, s=30, cmap=red_to_green, vmin=0,vmax=1,edgecolor='')

plt.xlim([-3.141592,3.141592])
plt.ylim([-3.141592/2,3.141592/2])

n=f.split("_")[1]
n=n.split(".")[0]

plt.savefig(n+".png",bbox_inches='tight',dpi=200)

counter=1

print ""

for f in listconf:

	n=f.split("_")[1]
	n=n.split(".")[0]

	print "\033[F"+"                                                                                 "
	print "\033[F"+"Processing file "+f+" ("+str(counter+1)+" of "+str(len(listconf))+")"

	phi = np.transpose(np.loadtxt(f))

	ax.scatter(coord2D[0], coord2D[1], c=(phi[0]+1)/2, s=30, cmap=red_to_green, vmin=0,vmax=1,edgecolor='')
	plt.savefig(n+".png",bbox_inches='tight',dpi=120)

	plt.cla()

	counter=counter+1

