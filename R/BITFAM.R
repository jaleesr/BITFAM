#' BITFAM main function. BITFAM will infer the transcription factor activities from single cell RNA-seq data based on the ChIP-seq data
#'
#' @param data A matrix or dataframe, normalized single cell RNA-seq data
#' @param species mouse or human
#' @param interseted_TF Transcription factors of interests
#' @param scATAC_obj A preprocessed Seurat object of scATAC-seq data
#' @param number of CPU cores
#' @param number of max iteration
#' @param convergence tolerance on the relative norm of the objective
#' @return sampling results of TF inferred activities and TF-gene weights
#' @export
#' @import rstan
#' @import Seurat

BITFAM <- function(data, species, interseted_TF = NA, scATAC_obj = NA,ncores, iter = 8000, tol_rel_obj=0.005){
  if(species == "mouse"){
    TF_targets_dir <- "mouse/"
  }else if(species == "human"){
    TF_targets_dir <- "human/"
  }else{
    stop("The species must be either mouse or human.")
  }

  if(dim(data)[1] > 5000){
    variable_genes <- Seurat::FindVariableFeatures(data)
    variable_genes <- variable_genes[which(x = variable_genes[, 1, drop = TRUE] != 0), ]
    variable_genes <- variable_genes[order(variable_genes$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
    variable_genes <- head(x = rownames(x = variable_genes), n = 5000)
    data <- data[variable_genes, ]
  }
  
  All_TFs <-system.file("extdata", paste0(TF_targets_dir, "all_TFs.txt"), package = "BITFAM")
  All_TFs <- read.table(All_TFs, stringsAsFactors = F)$V1
  TF_used <- rownames(data)[rownames(data) %in% All_TFs]
  rownames(data) <- toupper(rownames(data))
  if(is.na(interseted_TF)){
  }else{
    TF_used <- unique(c(TF_used, interseted_TF))
  }
  gene_list <- list()
  for(i in TF_used){
    TF_targets_path <-system.file("extdata", paste0(TF_targets_dir, i), package = "BITFAM")
    tmp_gene <- read.table(TF_targets_path, stringsAsFactors = F)
    tmp_gene <- toupper(tmp_gene$V1)
    gene_list[[which(TF_used == i)]] <- rownames(data)[rownames(data) %in% tmp_gene]
  }

  TF_used <- TF_used[ unlist(lapply(gene_list, length)) > 10]

  gene_list <- list()
  for(i in TF_used){
    TF_targets_path <-system.file("extdata", paste0(TF_targets_dir, i), package = "BITFAM")
    tmp_gene <- read.table(TF_targets_path, stringsAsFactors = F)
    tmp_gene <- toupper(tmp_gene$V1)
    gene_list[[which(TF_used == i)]] <- rownames(data)[rownames(data) %in% tmp_gene]
  }
  
  if(is.na(scATAC_obj)){
  }else{
    for(i in TF_used){
      gene_list[[which(TF_used == i)]] <- gene_list[[which(TF_used == i)]][gene_list[[which(TF_used == i)]] %in% BITFAM_scATAC(scATAC_obj)]
    }
  }
  
  X <- t(as.matrix(data))
  chipseq_weight <- matrix(1, nrow = length(colnames(X)), ncol = length(TF_used))
  for(i in 1:length(TF_used)){
    chipseq_weight[, i] <- ifelse(colnames(X) %in% gene_list[[i]], 1, 0)
  }


  Mask_matrix <- chipseq_weight
  X <- t(as.matrix(data))
  N <- dim(X)[1]
  D <- dim(X)[2]
  K <- length(TF_used)
  data_to_model <- list(N = N, D = D, K = K, X = X, Mask = Mask_matrix)


  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = ncores)

  set.seed(100)
  pca_beta_piror <- "
data {
int<lower=0> N; // Number of samples
int<lower=0> D; // The original dimension
int<lower=0> K; // The latent dimension
matrix[N, D] X; // The data matrix
matrix[D, K] Mask; // The binary mask of prior knowledge indicate the target of TFs
}

parameters {
matrix<lower=0, upper=1>[N, K] Z; // The latent matrix
matrix[D, K] W; // The weight matrix
real<lower=0> tau; // Noise term
vector<lower=0>[K] alpha; // ARD prior
}

transformed parameters{
matrix<lower=0>[D, K] t_alpha;
real<lower=0> t_tau;
for(wmd in 1:D){
for(wmk in 1:K){
t_alpha[wmd, wmk] = Mask[wmd, wmk] == 1 ? inv(sqrt(alpha[wmk])) : 0.01;
}
}
t_tau = inv(sqrt(tau));
}
model {
tau ~ gamma(1,1);
to_vector(Z) ~ beta(0.5, 0.5);
alpha ~ gamma(1e-3,1e-3);
for(d in 1:D){
for(k in 1:K){
W[d,k] ~ normal(0, t_alpha[d, k]);
}
}
to_vector(X) ~ normal(to_vector(Z*W'), t_tau);
} "

  m_beta_prior <- stan_model(model_code = pca_beta_piror)
  suppressWarnings(fit.vb <- vb(m_beta_prior, data = data_to_model, algorithm = "meanfield",
                                  iter = iter, output_samples = 300, tol_rel_obj = tol_rel_obj))
  BITFAM_list <- list(Model = fit.vb,
                      TF_used = TF_used,
                      Genes = rownames(data),
                      Cell_names = colnames(data))
  return(BITFAM_list)
}





