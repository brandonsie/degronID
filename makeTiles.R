source("imports.R")
source("tile.R")

Sys.time()

# load arguments
args <- commandArgs(TRUE)
fasta_path <- args[[1]] %>% as.character()
tile_length <- args[[2]] %>% as.numeric()
output_name <- args[[3]] %>% as.character()
input_seqs <- Biostrings::readAAStringSet(fasta_path)

print(paste0("fasta path: ", fasta_path))
print(paste0("tile_length: ", tile_length))
print(paste0("output name: ", output_name))


# call tile() function
output_tiles <- tile(input_seqs, tile_length, parallel = FALSE)
Biostrings::writeXStringSet(output_tiles, output_name)
Sys.time()

