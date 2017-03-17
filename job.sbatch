#!/bin/bash
#SBATCH -p weak-computation # queue (aka partition) name
#SBATCH -c 1 # for multithreaded applications
#SBATCH --mem-per-cpu 1800 # memory
#SBATCH -o slurm_%A_%a.out # stdout
#SBATCH -e slurm_%A_%a.out # stderr
#SBATCH --mail-type=FAIL,START,END # notifications for job start,end and failure
#SBATCH --mail-user=fonda@lorentz.leidenuniv.nl # send-to address

index=${SLURM_ARRAY_TASK_ID}
echo "Launching Job ${index}"

srun membrane -e .08 -m ../mesh/cube.msh -t 30 -I 5 -r $RANDOM $1 -C 0 $2 0 -L 0 -P 0 0 0 -T 1E-6 > /dev/null

