#' adjust_scoring_matrix
#' set max value for scoring matrix. if only one value provided, apply to all coloumns. if multiple values provided, use one for each column
#' @export

adjust_scoring_matrix <- function(scoring_matrix, col_max = 1){
  # update scoring matrix based on max value specification
  if(!is.null(col_max)){
    if(length(col_max) == 1){
      # if one value provided, apply to all columns
      scoring_matrix[scoring_matrix > col_max] <- col_max
    } else if(length(col_max) != ncol(scoring_matrix)){
      # if multiple values provided, first check that length is equal to number of columns in scoring matrix.
      error("Error: (pep_scoring): col_max must be a numeric of either length (1) or length equal to number of columns in scoring_matrix")
    } else{
      #then apply to each column in corresponding order
      for(i in 1:length(col_max)){
        scoring_matrix[,i][scoring_matrix[,i] > col_max[i]] <- col_max[i]
      }
    }
  }
  scoring_matrix


}
