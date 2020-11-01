#' Extract the TF inferred activities from the output of BITFAM main function
#'
#' @param BITFAM_list A list generate by BITFAM main function
#' @return TF inferred activities of each cell

BITFAM_activities <- function(BITFAM_list){
  result_matrix <- apply(extract(BITFAM_list@fit.vb,result,"Z")[[1]], c(2,3), mean)
  rownames(result_matrix) <- BITFAM_list@Genes
  colnames(result_matrix) <- BITFAM_list@Cell_names
  return(result_matrix)
}
