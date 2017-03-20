#!/bin/bash
#SBATCH -p weak-computation # queue (aka partition) name
#SBATCH -c 1 # for multithreaded applications
#SBATCH --mem-per-cpu 1800 # memory
#SBATCH -o slurm_%A_%a.out # stdout
#SBATCH -e slurm_%A_%a.out # stderr
#SBATCH --mail-type=FAIL # notifications for job start,end and failure
#SBATCH --mail-user=fonda@lorentz.leidenuniv.nl # send-to address

aindex=${SLURM_ARRAY_TASK_ID}
jindex=${SLURM_ARRAY_JOB_ID}
echo "Launching task ${aindex} of job ${jindex} with arguments $1 and $2 (which is the $3 of the series)"

srun membrane -e .08 -m ../mesh/cube.msh -t 30 -I 5 -r $RANDOM $1 -C 0 $2 0 -L 0 -P 0 0 0 -T 1E-6 > /dev/null

