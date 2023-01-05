source("pep_scoring.R")
`%>%` <- magrittr::`%>%`
args <- commandArgs(TRUE)
tile_path <- args[[1]] %>% as.character()
matrix_path <- args[[2]] %>% as.character()
pos_weight_path <- args[[3]] %>% as.character()
aa_weight_path <- args[[4]] %>% as.character()

# load tiles, ignore tiles with stop codon * (no longer ignore, updated pep_scoring 2022-11-15)
tile_set <- Biostrings::readAAStringSet(tile_path)
# tile_set <- tile_set[!grepl("\\*", tile_set)]


scoring_mat <- readr::read_rds(matrix_path)

if(pos_weight_path == "even"){
  pos_weight_vals <- 1
} else{
  pos_weight_vals <- readr::read_rds(pos_weight_path)
}
if(aa_weight_path == "even"){
  aa_weight_vals <- 1
} else{
  aa_weight_vals <- readr::read_rds(aa_weight_path)
}

output_name <- paste0(
  "scores__",
  tile_path %>% basename() %>% gsub("\\.[A-z]+$", "", .),
  "__",
  matrix_path %>% basename() %>% gsub("\\.[A-z]+$", "", .),
  "__",
  pos_weight_path %>% basename() %>% gsub("\\.[A-z]+$", "", .),
  "__",
  aa_weight_path %>% basename() %>% gsub("\\.[A-z]+$", "", .),
  ".csv"
)

print(tile_path)
print(length(tile_set))
print(head(tile_set))

print(matrix_path)
print(dim(scoring_mat))

print(pos_weight_path)
print(pos_weight_vals)
print(aa_weight_path)
print(aa_weight_vals)

print(output_name)

Sys.time()
output_scores <- pep_scoring(tiles = tile_set, scoring_matrix = scoring_mat, pos_weight = pos_weight_vals, aa_weight = aa_weight_vals, method = "mean", pb = TRUE)



# data.table::fwrite(output_scores, output_name)
Sys.time()


print("writing data subsets")



this_scores <- output_scores %>% tibble::as_tibble()

# add protein information
# if difffernent split if 4 or 2 field delimited by doubleunderscore __
num_fields <- this_scores$name[1] %>% strsplit("__") %>% unlist %>% length

print(paste(num_fields))
print(this_scores$name[1] %>% strsplit("__") %>% unlist)
print(this_scores$name[2] %>% strsplit("__") %>% unlist)

if(num_fields == 4){
    this_scores_abovethreshold <- this_scores %>%
    dplyr::mutate(
      entry = name %>% gsub("^(.*)__(.*)__(.*)__([0-9]+)$", "\\1", .),
      # sublibrary = name %>% gsub("^(.*)__(.*)__(.*)__([0-9]+)$", "\\2", .),
      gene = name %>% gsub("^(.*)__(.*)__(.*)__([0-9]+)$", "\\3", .),
      start = name %>% gsub("^(.*)__(.*)__(.*)__([0-9]+)$", "\\4", .),
    ) %>%
    dplyr::relocate(
      c(entry, #sublibrary, 
        gene, start), .after = name
    ) %>%
    dplyr::arrange(tile_score) #_norm_min
} else if(num_fields == 3){ #koren
    this_scores_abovethreshold <- this_scores %>%
    dplyr::mutate(
      gene = name %>% gsub("^(.*)__(.*)__([0-9]+)$", "\\1", .),
      entry = name %>% gsub("^(.*)__(.*)__([0-9]+)$", "\\2", .),
      start = name %>% gsub("^(.*)__(.*)__([0-9]+)$", "\\3", .),
    ) %>%
    dplyr::relocate(
      c(entry, #sublibrary, 
        gene, start), .after = name
    ) %>%
    dplyr::arrange(tile_score) #_norm_min
} else if(num_fields == 2){
    this_scores_abovethreshold <- this_scores %>%
    dplyr::mutate(
      entry = name %>% gsub("^(.*)__(.*)$", "\\1", .),
      start = name %>% gsub("^(.*)__(.*)$", "\\2", .),
    ) %>%
    dplyr::relocate(
      c(entry, start), .after = name
    ) %>%
    dplyr::arrange(tile_score) #_norm_min

} else{
    this_scores_abovethreshold <- this_scores
    
}



# write top peptides
data.table::fwrite(this_scores_abovethreshold, paste0("all_tiles_", output_name,".gz"))

# write top for each protein summary
if(num_fields == 4){
    this_protscores_ab_th <- this_scores_abovethreshold %>%
    dplyr::group_by(entry) %>%
    dplyr::summarise(
      best_score = min(tile_score), #_norm_min
      best_fragment_name = name[tile_score == best_score], #_norm_min
      # sublibrary = sublibrary[tile_score == best_score], #_norm_min
      gene = gene[tile_score == best_score], #_norm_min
      best_start = start[tile_score == best_score], #_norm_min
      best_tile = tile[tile_score == best_score] #_norm_min
    ) %>%
    dplyr::arrange(best_score)
} else if(num_fields == 3){
    this_protscores_ab_th <- this_scores_abovethreshold %>%
    dplyr::group_by(entry) %>%
    dplyr::summarise(
      best_score = min(tile_score), #_norm_min
      best_fragment_name = name[tile_score == best_score], #_norm_min
      # sublibrary = sublibrary[tile_score == best_score], #_norm_min
      best_start = start[tile_score == best_score], #_norm_min
      best_tile = tile[tile_score == best_score] #_norm_min
    ) %>%
    dplyr::arrange(best_score)
    
} else if(num_fields == 2){
    this_protscores_ab_th <- this_scores_abovethreshold %>%
    dplyr::group_by(entry) %>%
    dplyr::summarise(
      best_score = min(tile_score), #_norm_min
      best_fragment_name = name[tile_score == best_score], #_norm_min
      # sublibrary = sublibrary[tile_score == best_score], #_norm_min
      best_start = start[tile_score == best_score], #_norm_min
      best_tile = tile[tile_score == best_score] #_norm_min
    ) %>%
    dplyr::arrange(best_score)
    
}

data.table::fwrite(this_protscores_ab_th, paste0("proteins_", output_name,".gz"))


# write unique tiles (with an example of what protein)
top_unique_tiles <- this_scores_abovethreshold %>%
dplyr::group_by(tile) %>%
dplyr::summarise(
  # tile_score_norm_min = tile_score_norm_min[1],
  tile_score = tile_score[1],
  fragment_ids = paste0(name, collapse = "___")
) %>%
dplyr::arrange(tile_score) #_norm_min
data.table::fwrite(top_unique_tiles, paste0("unique_tiles_", output_name,".gz"))




Sys.time()