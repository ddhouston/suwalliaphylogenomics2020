#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH -t 96:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 6   # total number of processor cores in this job
#SBATCH -J "iq-suw-50p"   # job name
#SBATCH --no-requeue
#SBATCH --mail-user=jsatler@iastate.edu   # email address
#SBATCH --mail-type=END


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#code
cd mafft-nexus-edge_internal-trimmed-gblocks-clean-50p/

for i in *.nexus; do iqtree -s $i -bb 1000 -nt AUTO -ntmax 6; done

mkdir ../trees_50p

cp *.treefile ../trees_50p


