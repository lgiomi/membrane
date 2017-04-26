#!/bin/bash
#SBATCH -p lowmem # queue (aka partition) name
#SBATCH -c 1 # for multithreaded applications
#SBATCH --mem-per-cpu 1800 # memory
#SBATCH -o slurm_%A_%a.out # stdout
#SBATCH -e slurm_%A_%a.out # stderr
#SBATCH --mail-type=FAIL # notifications for job start,end and failure
#SBATCH --mail-user=fonda@lorentz.leidenuniv.nl # send-to address
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

aindex=${SLURM_ARRAY_TASK_ID}
jindex=${SLURM_ARRAY_JOB_ID}

echo "Launching task ${aindex} of job ${jindex} with arguments $1 and $2 (which is the $3 of the series)" >> ../run_stat.dat

srun membrane -m ../mesh/cube_2.msh -t 120 -I 5 -r $RANDOM $1 -C 0 $2 0 -L 0 -a 1 -k 2 > /dev/null

