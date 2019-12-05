#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=24GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=rk2805@nyu.edu
#SBATCH --output=rajeev_%j.out
  
module purge
module load matlab/2019a

RUNDIR=$SCRATCH/FIBR/HO_Simulations/
cd $RUNDIR

matlab -nodisplay -nodesktop -r "run SimulationLOS.m"
