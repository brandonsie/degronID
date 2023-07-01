# DegronID 



[![DOI](https://zenodo.org/badge/585416651.svg)](https://zenodo.org/badge/latestdoi/585416651)


This repository contains original code used in (paper title) Zhang et al. 2023 to score similarity between finely mapped degron sequences. 

More specifically, `DegronID` takes as input (1) a position specific scoring matrix (PSSM) describing saturation mutagenesis of a degron peptide motif and (2) a database of amino acid sequences of interest. For each amino acid sequence in the query database, DegronID generates a numeric score representing how similar or different the sequences is from the degron motif described by the PSSM. The details of this scoring are described in Zhang et al.

In our experience, degrons that are similar to each other will tend to overlap with each other in which query sequences from the database are scored favorably. Thus, the output from DegronID can be harnessed to cluster degron motifs into related groups. Illustration of the clustering result from Zhang et al. and related visualizations are accessible at https://brandonsie.shinyapps.io/DegronID/.

## Basic Usage
A walkthrough with example data is provided in `example.R`. 

## SLURM Usage
This algorithm was written in R and run on a computing cluster that uses Slurm Workload Manager. Necessary R files and example shell scripts are included.  

- Creating a query database of amino acid sequences.   Save sequences as a .fasta file. `example.fasta`
- Generate all tiles of each sequence of length equal to the length of the PSSM (in this example, length = 12aa).  

```{bash}
sbatch --export=fasta="example.fasta",len=12,out="example_tiles_12.fasta" tile.sh
```

- Store position specific scoring matrix as an .rds file. Save a position specific scoring matrix describing a degron footprint as an .rds file. Each column should represent a position in the amino acid sequence, and each row should represent a possible amino acid substitution. Rows should be named with the single letter code for the corresponding amino acid.  

- Score tiles from the query database against the PSSM of interest.  
```{bash}
sbatch --export=tile="example_tiles_12.fasta",mat="example_pssm.rds",pos_weight="even",aa_weight="even" compute_scores.sh
``` 


## Other Applications 

Although we used this method to characterize saturation mutagenesis of degrons, this algorithm could conceivably be applied to various other applications involving saturation mutagenesis, such as for investigation of TCR cross reactivity.
