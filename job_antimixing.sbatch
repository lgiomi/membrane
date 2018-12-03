#!/bin/bash
#SBATCH -p computation # queue (aka partition) name
#SBATCH -c 1 # for multithreaded applications
#SBATCH --mem-per-cpu 1800 # memory
#SBATCH -o slurm_%A_%a.out # stdout
#SBATCH -e slurm_%A_%a.out # stderr
#SBATCH --mail-type=FAIL # notifications for job start,end and failure
#SBATCH --mail-user=fonda@lorentz.leidenuniv.nl # send-to address
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Launching task with arguments ($1 ; $2, $3, $4, $5, $6, $7)" >> ../run_stat.dat

srun membrane -m ../mesh/dumbbell_13.msh -t 200 -I 5 -r $RANDOM $1 -L 0 -M -u $2 $3 $4 $5 $6 $7 -w 100
