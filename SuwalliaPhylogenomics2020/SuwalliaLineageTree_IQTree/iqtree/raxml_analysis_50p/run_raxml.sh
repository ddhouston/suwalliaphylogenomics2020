#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH -t 96:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 32   # total number of processor cores in this job
#SBATCH --mem=64G
#SBATCH -J "raxml-unpart-all-50p"   # job name
#SBATCH --no-requeue
#SBATCH --mail-user=jsatler@iastate.edu   # email address
#SBATCH --mail-type=END


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load raxml/8.2.11-openmpi3-ffoixg2
#code

raxmlHPC-PTHREADS-SSE3 \
              -f a -m GTRGAMMA -N 100 -T 32 \
              -x 12345 -p 12345 \
              -s mafft-nexus-edge_internal-trimmed-gblocks-clean-50p.phylip \
              -n results_unpart


