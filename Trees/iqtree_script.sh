#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module load Bioinformatics/Software/vital-it
module add Phylogeny/iqtree/1.7.beta17;

iqtree -s $1 -m LG -bb 1000 -nt 8
