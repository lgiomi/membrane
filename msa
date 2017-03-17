#!/bin/bash

rm -f membrane
rm -f *.dat
rm -f *.png
rm -f *.m
rm -f *.mkv
rm -fR finalconfigurations
make

# Save current dir in the PATH environment
DIR=`pwd`
export PATH=${DIR}:$PATH

# Total number of steps in the concentration, going from 0 to 1 (both values included, so the total number is acutally c0steps+1). 
n_c0=3

# Coupling counter (from -g_i to g_i with n_g steps)
n_g=2
g_i=.1

# Number of runs per phase diagram point
rpp=2

# Initialising counters
c0_c=0
g_c=0
rpp_c=0

while [ $g_c -le $n_g ];
	do
	g=$(echo "scale=3; -$g_i+$g_c/$n_g*2*$g_i" |bc)
	while [  $c0_c -le $n_c0 ]; 
	do
		while [  $rpp_c -lt $rpp ]; 
		do
		c0=$(echo "scale=3; $c0_c/$n_c0" |bc)

		mkdir -p temp_${g}_${c0}_${rpp_c}
		cd temp_${g}_${c0}_${rpp_c}

		sleep .5
		membrane -e .08 -m ../mesh/cube.msh -t 30 -I 5 -r $RANDOM $c0 -C 0 $g 0 -L 0 -P 0 0 0 -T 1E-6
		echo "Parameters at" $c0 $g $rpp_c 

		cd ..
		let rpp_c=rpp_c+1
		done 
	rpp_c=0
	let c0_c=c0_c+1 
	done
	c0_c=0
let g_c=g_c+1
done

echo "Finished running. Now collecting data together."

mkdir -p finalconfigurations

# Re-initialising counters
c0_c=0
g_c=0
rpp_c=0

while [ $g_c -le $n_g ];
	do
	g=$(echo "scale=3; -$g_i+$g_c/$n_g*2*$g_i" |bc)
	while [  $c0_c -le $n_c0 ]; 
	do
		while [  $rpp_c -lt $rpp ]; 
		do
		c0=$(echo "scale=3; $c0_c/$n_c0" |bc)

		mv temp_${g}_${c0}_${rpp_c}/last.dat finalconfigurations/${g}_${c0}_${rpp_c}_last.dat
		cat temp_${g}_${c0}_${rpp_c}/final.dat | awk -v var1="$g" -v var2="$c0"  -F"\t" '{print var1,"\t", var2,"\t", $2,"\t", $3,"\t", $4,"\t", $5,"\t", $6,"\t", $7,"\t", $8,"\t",$9}'  >> finalconfigurations/final.dat

		let rpp_c=rpp_c+1
		done 
	rpp_c=0
	let c0_c=c0_c+1 
	done
	c0_c=0
let g_c=g_c+1
done

rm -R temp_*
