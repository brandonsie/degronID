#!/bin/bash
#SBATCH --job-name=compute_scores 
#SBATCH --time=0:45:0 
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL

echo $tile
echo $mat
echo $pos_weight
echo $aa_weight
module load gcc/9.2.0 R/4.2.1 #updated from R/4.1.2
Rscript compute_scores.R $tile $mat $pos_weight $aa_weight

# --export=tile="example_tiles_12.fasta",mat="example_pssm.rds",pos_weight="even",aa_weight="even" compute_scores.sh
