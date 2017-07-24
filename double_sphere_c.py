import matplotlib
import numpy as np
import sys,glob,os
from math import sqrt,pi

if len(sys.argv)<3:
    sys.exit("Syntax: python double_sphere_c.py geometry.dat folder/")

if os.path.exists(sys.argv[2]):
	print "Opening folder "+sys.argv[2]+""
else:
	sys.exit("Error: input directory does not exist")

listconf= sorted(glob.glob(sys.argv[2]+"/*last.dat"))

listk = list(set([i.split("_")[-5] for i in listconf]))
listk = sorted([float(i) for i in listk])
listkb = list(set([i.split("_")[-4] for i in listconf]))
listkb = sorted([float(i) for i in listkb])
listc = list(set([i.split("_")[-3] for i in listconf]))
listc = sorted([float(i) for i in listc])
listn = list(set([i.split("_")[-2] for i in listconf]))
listn = sorted([int(i) for i in listn])

try:
	xAvalues= np.loadtxt(sys.argv[1],usecols=(1,4),delimiter='\t')
except:
	sys.exit("Error: geometry file does not exist")

A1=0
A2=0

for i in xAvalues:

	if i[0]>0:
		A1+=i[1]
	else:
		A2+=i[1]

R1=sqrt(A1/4/pi)
R2=sqrt(A2/4/pi)

print "The mesh has",len(xAvalues),"mesh points. The two spheres have areas (",A1,",",A2, ") and radii (",R1,",",R2,")"
print ""

counter=0

dsfile = open('double_sphere_p.dat', 'w')

for f in listconf:
	c1=0
	c2=0
	n=f.split("_")
	print "\033[F"+"Processing file t_"+n[-5]+"_"+n[-4]+"_"+n[-3]+"_"+n[-2]+" ("+str(counter+1)+" of "+str(len(listconf))+")"
	counter=counter+1

	phi = np.loadtxt(f)

	for i in range(len(xAvalues)):
		if xAvalues[i,0]>0:
			c1+=xAvalues[i,1]*phi[i,0]
		else:
			c2+=xAvalues[i,1]*phi[i,0]
	dsfile.write(n[-5]+"\t"+n[-4]+"\t"+n[-3]+"\t"+str(c1)+"\t"+str(c2)+os.linesep)

dsfile.close()




