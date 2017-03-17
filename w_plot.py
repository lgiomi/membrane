#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main():
	coord3D = np.transpose(np.loadtxt("w_3D.dat",delimiter='\t'))
	coord2D = np.transpose(np.loadtxt("w_2D.dat",delimiter='\t'))
	#fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	#plt.scatter(coord3D[0], coord3D[1], coord3D[2], edgecolor='')
	#data = np.memmap(filename, datatype, 'r') 
	#plt.plot(data['floati'],data['floatq'],'r,')
	#plt.show()

if __name__ == "__main__":
	main() 

