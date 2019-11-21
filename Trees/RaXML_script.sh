#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --time 72:00:00
#SBATCH --partition=long
module load Bioinformatics/Software/vital-it
module add Phylogeny/raxml/8.2.12;

raxmlHPC-PTHREADS -T 15 -f a -m PROTGAMMALG -p 12345 -x 12345 -# 100 -s $1 -n $2
