# example usage ofthe DegronID scoring algorithm in R


# read fasta of amino acid sequences of interest
input_seqs <- Biostrings::readAAStringSet("data/example.fasta")

# convert aa sequences into database of peptide tiles of length L (for this example, L=12)
source("imports.R")
source("tile.R")

tile_length <- 12
tile_database <- tile(input_seqs, tile_length, parallel = FALSE)
Biostrings::writeXStringSet(tile_database, "data/example_tiles_12.fasta")

# load rds pssm
this_pssm <- readRDS("data/example_pssm.rds")

# compute DegronID Scores
source("pep_scoring.R")

output_scores <- pep_scoring(tiles = tile_database, scoring_matrix = this_pssm, method = "mean", pb = TRUE)

# Examine best scoring output (lowest tile_score
`%>%` <- magrittr::`%>%`
output_scores %>% dplyr::arrange(tile_score)
