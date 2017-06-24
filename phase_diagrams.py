import matplotlib
import numpy as np
import scipy.interpolate
import sys,glob,os
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

if len(sys.argv)<2:
    sys.exit("Syntax: python phase_diagrams.py data_file.dat")

try:
	data = np.loadtxt(sys.argv[1])
except:
	sys.exit("Error: input file does not exist")

def tuples(A):
    try: return tuple(tuples(a) for a in A)
    except TypeError: return A

# The structure of final.dat is of 13 columns, with the values:
# 1	|2	|3	|4		|5			|6		|7		|8		|9		|10		|11		|12		|13
# eta_k	|eta_kb	|c0	|kinetic_energy	|phase_field_potential	|full_potential	|<phi^2>_c	|<phi H^2>_c	|<phi K_G>_c	|Lagrange	|N_domains	|total_time	|C0
# (note that c0 is an input from run_jobs and ranges [0,1], while C0 is computed from an integral of phi and ranges [-1,1])


#if not os.path.exists("pd_images"):
#    os.makedirs("pd_images")

listkkb = sorted(set(tuples(data[:,0:2])))
listkkb = [list(t) for t in listkkb]

sdata=np.zeros((len(listkkb),len(data[0])))
n=0

for i in listkkb:
	nj=0
	for j in data:
		if j[0] == i[0] and j[1] == i[1]:
			sdata[n] = j+sdata[n]
			nj+=1
	sdata[n]/=nj
	n+=1

X = np.linspace(min(sdata[:,1]), max(sdata[:,1]), 100)
Y = np.linspace(min(sdata[:,0]), max(sdata[:,0]), 100)

# In the older  version, no need of this line
X, Y = np.meshgrid(X, Y)

# In an older version of the code, we used to griddata to mesh the scatter plot:
#kinetic_energy 	= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,3], X, Y, interp='linear')
#phase_field_potential 	= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,4], X, Y, interp='linear')
#full_potential 	= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,5], X, Y, interp='linear')
#phi2c			= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,6], X, Y, interp='linear')
#phiH2			= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,7], X, Y, interp='linear')
#phiKG			= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,8], X, Y, interp='linear')
#lagrange		= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,9], X, Y, interp='linear')
#N_domains		= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,10], X, Y, interp='linear')
#total_time		= mlab.griddata(sdata[:,1], sdata[:,0], sdata[:,11], X, Y, interp='linear')

# Which then could have been used to plot contours as:
# plt.pcolormesh(X, Y, phiH2, cmap = plt.get_cmap('viridis'))
# contour_phiH2 = plt.contour(X, Y, phiH2, n_contours, inline=1, fontsize=font_contours, linewidths=lw, colors='k')
# contour_N_domains = plt.contour(X, Y, N_domains, [1,2,3,4,5,6], inline=1, fontsize=font_contours, linewidths=lw, colors='k')
# plt.clabel(contour_phiH2, inline=1, fmt='%2.2f',fontsize=font_contours)
# plt.title(r"$\left<\phi H^2\right>_c$", fontsize=font_title, color='black')
# plt.xlabel(r"$\eta_{\bar{k}}$",fontsize=font_axeslabel)
# plt.ylabel(r"$\eta_k}$",fontsize=font_axeslabel)

N_domains_rbf 		= scipy.interpolate.Rbf(sdata[:,1], sdata[:,0], sdata[:,10], function='linear')
N_domains 		= N_domains_rbf(X,Y)
# Alternatively, griddata could have been used:
#N_domains 		= scipy.interpolate.griddata((sdata[:,1], sdata[:,0]), sdata[:,10], (X, Y), method='cubic')

kinetic_energy_rbf 	= scipy.interpolate.Rbf(sdata[:,1], sdata[:,0], sdata[:,3], function='thin_plate')
kinetic_energy 		= kinetic_energy_rbf(X,Y)

phiH2_rbf 		= scipy.interpolate.Rbf(sdata[:,1], sdata[:,0], sdata[:,7], function='linear')
phiH2 			= phiH2_rbf(X,Y)

phiKG_rbf 		= scipy.interpolate.Rbf(sdata[:,1], sdata[:,0], sdata[:,8], function='linear')
phiKG	 		= phiKG_rbf(X,Y)

full_plot = plt.figure(facecolor='white')
plt.rc('text', usetex=True)
font_title=16
font_contours=6
font_axeslabel=12
lw=.8
n_contours=6

full_plot.add_subplot(221)

plt.title(r"$\mathrm{Number\;of\;domains}$", fontsize=font_title, color='black')
plt.xlabel(r"$\eta_{\bar{k}}$",fontsize=font_axeslabel)
plt.ylabel(r"$\eta_k}$",fontsize=font_axeslabel)
plt.imshow(N_domains, cmap='viridis', origin='lower', extent=[sdata[:,1].min(), sdata[:,1].max(), sdata[:,0].min(), sdata[:,0].max()])
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax, ticks=[1,2,3,4])

full_plot.add_subplot(222)

plt.title(r"$\mathrm{Interface\;length}$", fontsize=font_title, color='black')
plt.xlabel(r"$\eta_{\bar{k}}$",fontsize=font_axeslabel)
plt.ylabel(r"$\eta_k}$",fontsize=font_axeslabel)
plt.imshow(kinetic_energy, cmap='viridis', origin='lower', extent=[sdata[:,1].min(), sdata[:,1].max(), sdata[:,0].min(), sdata[:,0].max()])
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)

full_plot.add_subplot(223)

plt.title(r"$\left<\phi H^2\right>_c$", fontsize=font_title, color='black')
plt.xlabel(r"$\eta_{\bar{k}}$",fontsize=font_axeslabel)
plt.ylabel(r"$\eta_k}$",fontsize=font_axeslabel)
plt.imshow(phiH2, cmap='viridis', origin='lower', extent=[sdata[:,1].min(), sdata[:,1].max(), sdata[:,0].min(), sdata[:,0].max()])
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)

full_plot.add_subplot(224)

plt.title(r"$\left<\phi K_G\right>_c$", fontsize=font_title, color='black')
plt.xlabel(r"$\eta_{\bar{k}}$",fontsize=font_axeslabel)
plt.ylabel(r"$\eta_k}$",fontsize=font_axeslabel)
plt.imshow(phiKG, cmap='viridis', origin='lower', extent=[sdata[:,1].min(), sdata[:,1].max(), sdata[:,0].min(), sdata[:,0].max()])
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)

#plt.show()
plt.tight_layout()
plt.savefig("report.pdf",bbox_inches='tight',dpi=80)
