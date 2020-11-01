BITFAM_extract <- function(BITFAM_list, result = "Z"){
  result_matrix <- apply(extract(stan.fit.vb.real.beta.prior,result)[[1]], c(2,3), mean)
  return(result_matrix)
}