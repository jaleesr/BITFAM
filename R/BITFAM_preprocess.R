#' Generate the normalized scRNA-seq data by Seurat 3.0
#'
#' @param raw_data A matrix or dataframe of scRNA-seq raw counts table
#' @return Normalized scRNA-seq data
#' @export
#' @import rstan
#' @import Seurat

BITFAM_preprocess <- function(raw_data){
  process_data <- CreateSeuratObject(counts = raw_data, min.cells = 3, min.features = 200)
  process_data <- NormalizeData(object = process_data)
  process_data <- FindVariableFeatures(object = process_data, nfeatures = 5000)

  data_normalized <- as.matrix(GetAssayData(object = process_data)[VariableFeatures(process_data), ])
  rownames(data_normalized) <- VariableFeatures(process_data)
  colnames(data_normalized) <- colnames(GetAssayData(object = process_data))
  return(data_normalized)
}
