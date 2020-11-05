#' Extract the TF-gene weights from the output of BITFAM main function
#'
#' @param BITFAM_list A list generate by BITFAM main function
#' @return TF-gene weights for each pair
#' @export
#' @import rstan

BITFAM_weights <- function(BITFAM_list){
  result_matrix <- apply(extract(BITFAM_list$Model,result,"W")[[1]], c(2,3), mean)
  rownames(result_matrix) <- BITFAM_list$Genes
  colnames(result_matrix) <- BITFAM_list$TF_used
  return(result_matrix)
}
