import matplotlib
import numpy as np
import scipy.interpolate
import sys,glob,os
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def tuples(A):
    try: return tuple(tuples(a) for a in A)
    except TypeError: return A

try:
	data = np.loadtxt("fig6a.dat")
except:
	sys.exit("Error: input file does not exist")

X = np.linspace(min(data[:,0]), max(data[:,0]), 100)
Y = np.linspace(min(data[:,1]), max(data[:,1]), 100)
X, Y = np.meshgrid(X, Y)

#print X,Y

kinetic_energy_rbf 	= scipy.interpolate.Rbf(data[:,0], data[:,1], data[:,2], function='quintic')
kinetic_energy 		= kinetic_energy_rbf(X,Y)

full_plot = plt.figure(facecolor='white')
plt.rc('text', usetex=True)
font_title=16
font_contours=6
font_axeslabel=12
lw=.8
n_contours=6
c_map='jet'

full_plot.add_subplot(111)

#plt.title(r"$\mathrm{Interface\;length}$", fontsize=font_title, color='black')
plt.xlabel(r"$\varphi_-$",fontsize=font_axeslabel)
plt.ylabel(r"$\eta_k$",fontsize=font_axeslabel)
plt.imshow(kinetic_energy, cmap=c_map, origin='lower', extent=[data[:,0].min(), data[:,0].max(), data[:,1].min(), data[:,1].max()],vmin=0, vmax=2)
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
ax.set_aspect('auto')
plt.colorbar(cax=cax)

plt.savefig("Length_mathematica.pdf",bbox_inches='tight',dpi=80)



