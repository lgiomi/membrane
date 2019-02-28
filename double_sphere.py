import matplotlib
import numpy as np
import math
import sys,glob,os
from math import sqrt,pi

if len(sys.argv)<5 or len(sys.argv)>6:
    sys.exit("Syntax: python double_sphere.py geometry.dat X1 [X2] XYZ folder/")

if len(sys.argv)==5 and os.path.exists(sys.argv[4]):
	xA=float(sys.argv[2])
	xAf=xA
	dimA=int(sys.argv[3])
	print "Opening folder "+sys.argv[4]+""
	listconf= sorted(glob.glob(sys.argv[4]+"/*last.dat"))
elif len(sys.argv)==6 and os.path.exists(sys.argv[5]):
	xA=float(sys.argv[2])
	xAf=float(sys.argv[3])
	dimA=int(sys.argv[4])
	print "Opening folder "+sys.argv[5]+""
	listconf= sorted(glob.glob(sys.argv[5]+"/*last.dat"))
else:
	sys.exit("Error: input directory does not exist")

try:
	xAvalues= np.loadtxt(sys.argv[1],usecols=(1,2,3,4),delimiter='\t')
except:
	sys.exit("Error: geometry file does not exist")

A1=0
A2=0
A3=0

for i in xAvalues:

	if i[dimA]>xA:
		A1+=i[3]
	elif i[dimA]<xAf:
		A2+=i[3]
	else:
		A3+=i[3]

R1=sqrt(A1/4/pi)
R2=sqrt(A2/4/pi)
x1=A1/(A1+A2)
x2=A2/(A1+A2)

print "The mesh has",len(xAvalues)," mesh points."
print "Areas: (",A1,",",A2,",",A3,")"
print "Radii: (",R1,",",R2,")"
print "Area fractions: (",x1,",",x2,")"
print ""

counter=0

dsfile = open('extracted_phi_values.dat', 'w')

for f in listconf:
	phi1=0
	phi2=0
	phi3=0
	phiTot=0
	n=f.split("_")

	counter=counter+1

	phi = np.loadtxt(f)

	if math.isnan(phi[1,0]):
		continue

	for i in range(len(xAvalues)):
		phiTot+=xAvalues[i,3]*phi[i,0]
		if xAvalues[i,dimA]<xA:
			phi1+=xAvalues[i,3]*phi[i,0]
		elif xAvalues[i,dimA]>xAf:
			phi2+=xAvalues[i,3]*phi[i,0]
		else:
			phi3+=xAvalues[i,3]*phi[i,0]

	phi1=phi1/A1
	phi2=phi2/A2
	phi3=phi2/A3
	phiTot=phiTot/(A1+A2)

	phi1=(phi1+1)/2
	phi2=(phi2+1)/2
	phiTot=(phiTot+1)/2

	if math.isnan(phi1) or math.isnan(phi2):
		continue
	
	print "\033[F"+"Processing file "+str(counter)+" of "+str(len(listconf))
	dsfile.write(n[1]+"\t"+n[2]+"\t"+n[3]+"\t"+n[4]+"\t"+n[5]+"\t"+n[6]+"\t"+n[7]+"\t"+n[8]+"\t"+str(phi1)+"\t"+str(phi2)+"\t"+str(phiTot)+"\t"+str(x1)+"\t"+str(x2)+os.linesep)

dsfile.close()




