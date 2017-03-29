import numpy as np,sys
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

if len(sys.argv)<3:
    sys.exit("Syntax: python plot.py input.dat output.png")


coord3D = np.transpose(np.loadtxt("geometry.dat",usecols=(1,2,3),delimiter='\t'))
coord2D = np.transpose(np.loadtxt("geometry.dat",usecols=(7,8),delimiter='\t'))

triangles = np.loadtxt("triangles.dat")

try:
	phi = np.transpose(np.loadtxt(sys.argv[1]))
except:
    sys.exit("Error: input file does not exist")


cdict = {'red': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),'green':((0.0, 0.0, 0.0),(1.0, 1.0, 1.0)),'blue':((0.0, 0.0, 0.0),(1.0, 0.0, 0.0))}

red_to_green = LinearSegmentedColormap('BlueRed1', cdict)

fig, ax = plt.subplots()
ax.scatter(coord2D[0], coord2D[1], c=(phi[0]+1)/2, s=30, cmap=red_to_green, vmin=0,vmax=1,edgecolor='')

plt.xlim([-3.141592,3.141592])
plt.ylim([-3.141592/2,3.141592/2])

plt.savefig(sys.argv[2],bbox_inches='tight',dpi=200)
