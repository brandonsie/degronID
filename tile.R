#' tile.R
#' generate all possible peptide tiles from longer sequences.
#' @param input_seqs named character vector or AAStringSet of ORF amino acid sequences
#' @param tile_length number of amino acids to include per tile
#' @param parallel logical. indicates whether to use parallel processing.
#' @param collapse_list. logical. if true (default), return single aastringset with tiles named orfname1_1, orfname1_2, ..., orfname2_1 .... if false, return list of character vectors. each list element is orfname1, orfname2, ..., and each character vector element is named with amino acid start position.
#'
#' @export

tile <-
  function(input_seqs,
           tile_length,
           parallel = TRUE,
           collapse_list = TRUE) {

    # coerce to character
    input_class <- class(input_seqs)[1]
    if (input_class == "AAStringSet") {
      input_names <- names(input_seqs)
      input_seqs <- as.character(input_seqs)
      names(input_seqs) <- input_names
      input_class <- "character"
    }

    # remove orfs shorter than tile length
    input_seqs <- input_seqs[nchar(input_seqs) >= tile_length]



    if (parallel) {
      no_cores <- future::availableCores() - 1
      future::plan(workers = no_cores)

      tile_sample <- furrr::future_map(input_seqs,  function(x) {
        nchar_x <- x %>% as.character %>% nchar
        output_startpos <- 1:(nchar_x - tile_length + 1)
        output_subseq <- output_startpos %>% sapply(function(y) {
          x %>% Biostrings::subseq(start = y, width = tile_length)
        }) %>% Biostrings::AAStringSet()
        if (input_class == "character") {
          output_subseq <- output_subseq %>% as.character
        }
        names(output_subseq) <- c(1:length(output_subseq))
        output_subseq
      })


    } else{
      tile_sample <- purrr::map(input_seqs,  function(x) {
        nchar_x <- x %>% as.character %>% nchar
        output_startpos <- 1:(nchar_x - tile_length + 1)
        output_subseq <- output_startpos %>% sapply(function(y) {
          x %>% Biostrings::subseq(start = y, width = tile_length)
        }) %>% Biostrings::AAStringSet()
        if (input_class == "character") {
          output_subseq <- output_subseq %>% as.character
        }
        names(output_subseq) <- c(1:length(output_subseq))
        output_subseq
      })

    }


    # unlist output
    if (collapse_list) {
      tile_sample <- tile_sample %>% unlist()
      names(tile_sample) <-
        names(tile_sample) %>% gsub("\\.([0-9]+)$", "__\\1", .)
      tile_sample <- tile_sample %>% Biostrings::AAStringSet()
    } else{
      tile_sample
    }

  }
