#' Extract the genes that have scATAC-seq peaks on their promoter regions
#'
#' @param scATAC_obj A preprocessed Seurat object of scATAC-seq data.
#' @return the genes that have scATAC-seq peaks on their promoter regions
#' @export
#' @import rstan
#' @import Seurat

BITFAM_scATAC <- function(scATAC_obj){
	return(rownames(scATAC_obj))
}

