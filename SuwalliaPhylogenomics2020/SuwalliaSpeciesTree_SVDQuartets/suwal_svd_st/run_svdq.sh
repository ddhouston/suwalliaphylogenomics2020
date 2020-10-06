#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH -t 124:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 32   # total number of processor cores in this job
#SBATCH -J "svdq-50p"   # job name
#SBATCH --mail-user=jsatler@iastate.edu   # email address
#SBATCH --mail-type=END


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#code
paup mafft-nexus-edge_internal-trimmed-gblocks-clean-50p.nexus

