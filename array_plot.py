import numpy as np,sys
import matplotlib.pyplot as plt
import glob
from matplotlib.colors import LinearSegmentedColormap

listconf = glob.glob("testconf/*last.dat")

listg = list(set([i.split("_")[1] for i in listconf]))
listg = sorted([float(i) for i in listg])
listc = list(set([i.split("_")[2] for i in listconf]))
listc = sorted([float(i) for i in listc])
listn = list(set([i.split("_")[3] for i in listconf]))
listn = sorted([int(i) for i in listn])

coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
coord2D = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))

cdict = {'red': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),'green':((0.0, 0.0, 0.0),(1.0, 1.0, 1.0)),'blue':((0.0, 0.0, 0.0),(1.0, 0.0, 0.0))}
red_to_green = LinearSegmentedColormap('BlueRed1', cdict)

fig, ax = plt.subplots(frameon=False)

print "\n"
counter=0

for f in listconf:

	n=f.split("_")
	print "\033[F"+"                                                                                    "
	print "\033[F"+"Processing file f_"+n[1]+"_"+n[2]+"_"+n[3]+" ("+str(counter+1)+" of "+str(len(listconf))+")"
    	sys.stdout.flush()

	phi = np.transpose(np.loadtxt(f))

	fig.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0, wspace=0)

	plt.xlim([-3.141592,3.141592])
	plt.ylim([-3.141592/2,3.141592/2])

	ax.axis('off')

	ax.scatter(coord2D[0], coord2D[1], c=(phi[0]+1)/2, s=30, cmap=red_to_green, vmin=0,vmax=1,edgecolor='')

	fig.savefig("testconf/"+"f_"+str(listg.index(float(n[1])))+"_"+str(listc.index(float(n[2])))+"_"+str(listn.index(int(n[3])))+".png",pad_inches=0,dpi=200)

        plt.cla()

	counter=counter+1


print "\nThere are %g coupling points and %g concentration points and %g series" % len(listg),len(listc),len(listn)
