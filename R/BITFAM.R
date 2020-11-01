BITFAM <- function(data, species, interseted_TF = NA, ncores){
  if(species == "mouse"){
    TF_targets_dir <- "TF/mouse/"
  }else if(species == "human"){
    TF_targets_dir <- "TF/human/"
  }else{
    stop("The species must be either mouse or human.")
  }

  gene_list <- list()
  for(i in TF_used){
    tmp_gene <- read.table(paste0(TF_targets_dir, i), stringsAsFactors = F)
    gene_list[[which(TF_used == i)]] <- VariableFeatures(process_data)[VariableFeatures(process_data) %in% tmp_gene$V1]
  }

  TF_used <- TF_used[ unlist(lapply(gene_list, length)) > 10]
  if(is.na(interseted_TF)){
  }else{
    TF_used <- unique(c(TF_used, interseted_TF))
  }

  gene_list <- list()
  for(i in TF_used){
    tmp_gene <- read.table(paste0(TF_targets_dir, i), stringsAsFactors = F)
    gene_list[[which(TF_used == i)]] <- VariableFeatures(process_data)[VariableFeatures(process_data) %in% tmp_gene$V1]
  }

  data_matrix_normalized <- t(as.matrix(GetAssayData(object = process_data)[VariableFeatures(process_data), ]))
  data_matrix_normalized <- data_matrix_normalized[, -grep(pattern = "gRNA", x = VariableFeatures(process_data))]

  chipseq_weight <- matrix(1, nrow = length(colnames(data_matrix_normalized)), ncol = length(TF_used))
  for(i in 1:length(TF_used)){
    chipseq_weight[, i] <- ifelse(colnames(data_matrix_normalized) %in% gene_list[[i]], 1, 0)
  }


  Mask_matrix <- chipseq_weight
  X <- data_matrix_normalized
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
  stan.fit.vb.real.beta.prior <- vb(m_beta_prior, data = data_to_model, algorithm = "meanfield",
                                  iter = 8000, output_samples = 300)
  BITFAM_list <- list(Model = stan.fit.vb.real.beta.prior,
                      TF_used = TF_used,
                      Genes = VariableFeatures(process_data))
  return(BITFAM_list)
}




