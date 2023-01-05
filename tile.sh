#!/bin/bash
#SBATCH --job-name=tile_degronID 
#SBATCH --time=3:0:0
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL


module load gcc/9.2.0 R/4.2.1
Rscript makeTiles.R $fasta $len $out

# sbatch --export=fasta="example.fasta",len=12,out="example_tiles_12.fasta" tile.sh
