pep_scoring <- function(
  tiles, scoring_matrix, pos_weight = 1, aa_weight = 1, method = c("mean", "sum"), pb = FALSE,
  stop_value = 1 
){

  # take AAStringSet of tiles
  # convert to data frame
  # vector: name
  # vector: full sequence
  # matrix: column for each aa position (1-9)
  # matrix: score for each aa position (1-9)
  # vector: score for each row (peptide tile)
  # vector: score normalized to max score
  # join columns into data frame output

  # pos_weight either 1 or a numeric vector equal to motif length (number of columns of scoring_matrix)
  # aa_weight either 1 or numeric vector equal to amino acid alphabet length (20) indicating relative weight. element names = amino acid single letter codes, justlike matrix
  # for each substitution, take the max of pos_weight or aa_weight
  # method. mean mean scores of individual amino acid substitutions. sum sum scores of individual amino acid substitutions. 
  # pb show progress bar

  # if input AAStringSet, coerce to character vector.
  # then split each string into individual amino acids
  # matrix_aa is a character matrix. one row per tile. one column per amino acid postion
  
  # stop_value value to assign to stop codon
  
  matrix_aa <- tiles %>% as.character %>% strsplit(., "") %>%
    setNames(names(tiles)) %>%
    dplyr::bind_rows() %>% t %>% as.matrix
  colnames(matrix_aa) <- paste0("aa_", c(1:ncol(scoring_matrix)))

  # initialize numeric matrix matrix_score to store the amino acid substitution scores for each position in each tile
  matrix_score <- matrix(0, nrow = nrow(matrix_aa), ncol = ncol(matrix_aa))
  colnames(matrix_score) <- paste0("score_", c(1:ncol(scoring_matrix)))

  # prepare position weight variable. if single value provided, convert to numeric vector of length euqal to nchar(tiles[[1]])
  if(length(pos_weight) == 1){
      pos_weight <- rep(pos_weight, length(tiles[[1]]))
  }
  if(length(aa_weight) == 1){
            aa_weight <- rep(aa_weight, nrow(scoring_matrix)) %>% setNames(rownames(scoring_matrix))
  }

  # for each aa position (1-9), for each possible aa (1-20), input score into matrix_score
  if(pb) {prog <- txtProgressBar(0, ncol(scoring_matrix), style = 3)}
  for(pos in 1:ncol(scoring_matrix)){
    if(pb) {(setTxtProgressBar(prog, pos))}
    for(aa in rownames(scoring_matrix)){
      print(paste0(pos, ", ", aa))
      # define final weight as max of pos_weight and aa_weight
      final_weight <- max(aa_weight[aa], pos_weight[pos])
      
      # apply score
      matrix_score[,pos][matrix_aa[,pos] == aa] <- scoring_matrix[aa,pos] * final_weight
    }
  }
  
  # add default value for stop codon
  print("*")
  matrix_score[matrix_aa == "*"] <- stop_value


  if(method[1] == "mean"){
    tile_score <- matrix_score %>% apply(1,mean) #%>% apply(1, function(x) x %>% "*"(weight) %>% mean)
  }  else {
    tile_score <- matrix_score %>% apply(1,sum) #%>% apply(1, function(x) x %>% "*"(weight) %>% sum)
  } 
  
  



  output_df <- dplyr::bind_cols(
    name = names(tiles),
    tile = tiles %>% as.character,
    tile_score = tile_score,
    matrix_aa,
    matrix_score
  )

  output_df

}
