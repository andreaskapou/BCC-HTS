#' Generate synthetic data for the BCC model
#' 
#' BCC model is a generative model, so we can generate data from the model and 
#' then fit the model to the data back again. This function generates data from 
#' the following observation models: Normal (Gaussian), Binomial, Poisson and 
#' Binomial Probit Regresssion model (BPR).
#' 
#' @param N Total number of objects.
#' @param M Total number of sources.
#' @param K Total number of clusters.
#' @param C_true Ground truth overall clustering assignments.
#' @param alpha Ground truth adherence parameter.
#' @param distr Vector of characters of size M denoting the distribution that 
#'   each data source follows. Currently, 'G' for Normal, 'B' for Binomial, 'P' 
#'   for Poisson, and 'BPR' for the BPR model.
#' @param mu_init Mean parameters when generating data from Gaussian source.
#' @param std_init Standard deviation parameters when generating data from
#'   Gaussian source.
#' @param l_init Mean/Variance parameters when generating data from Poisson
#'   source.
#' @param p_init Probability of success parameters when generating data from
#'   Binomial source.
#' @param w_init Coefficient parameters when generating data from BPR source.
#' @param basis Basis object generated from BPRMeth package 
#'   (https://github.com/andreaskapou/BPRMeth) and is only required if you have 
#'   data that follow the BPR model, currently we support RBF (recommended) and 
#'   polynomial basis.
#' @param r The total number of trials, used only when we have source following
#'   a Binomial distribution.
#'   
#'   
gen_bcc_data <- function(N=200, M=2, K=2, C_true, alpha, distr = c("G","G"), 
                         mu_init = NULL, std_init = NULL, l_init = NULL, 
                         p_init = NULL, w_init = NULL, basis, r){
  B = P = G = BPR <- 1  # Initialize source counts
  
  #---------------------------------
  # Source specific clustering probabilities
  #---------------------------------
  nu_true  <- matrix(0, nrow=N, ncol=K)
  for(m in 1:M){
    for(k in 1:K){
      # Set to a when C=1
      nu_true[as.logical(C_true[,k]), k] = alpha
      # Set to (1-a)/(K-1) when C=0
      nu_true[!as.logical(C_true[,k]),k] = 1-alpha
    }
  }
  
  #---------------------------------
  # Generate source specific clusterings
  #---------------------------------
  L_true <- list()
  for (m in 1:M){
    L_true[[m]]  <- matrix(0, nrow=N, ncol=K)
  }
  for (m in 1:M){
    for (i in 1:N){  # Generate source specific clusterings
      L_true[[m]][i,] <- rmultinom(1, 1, nu_true[i,])
    }
  }
  
  #---------------------------------
  # Generate the objects X_mn based on L_mn
  #---------------------------------
  X       <- list()
  lambdas <- vector(length = N)
  bp      <- vector(length = N)
  mus     <- vector(length = N)
  stds    <- vector(length = N)
  if (!is.null(w_init)){
    ws    <- matrix(0, ncol = NROW(w_init[[1]]), nrow = N)
  }else{
    ws <- matrix(0, ncol = 1, nrow = N)
  }
  
  for (m in 1:M){   # For each source
    for (k in 1:K){ # For each cluster
      # Get indices where L_mn = 1
      index <- which(L_true[[m]][,k] == 1)
      # Assign parameter values according to distribution and cluster
      if (identical(distr[m],"G")){
        mus[index]  <- mu_init[[G]][k]
        stds[index] <- std_init[[G]][k]
      }else if (identical(distr[m],"P")){
        lambdas[index] <- l_init[[P]][k]
      }else if (identical(distr[m],"B")){
        bp[index] <- p_init[[B]][k]
      }else if (identical(distr[m],"BPR")){
        for (i in 1:length(index)){
          ws[index[i], ]  <- w_init[[BPR]][, k]
        }
      }
    }
    # Generate the simulated data for each distribution
    X[[m]] <- list()
    if (identical(distr[m],"G")){
#       X[[m]]$x <- cbind(as.matrix(stats::rnorm(n=N, mean=mus, sd=stds)), 
#                         as.matrix(stats::rnorm(n=N, mean=mus, sd=stds)))
      X[[m]]$x <- as.matrix(stats::rnorm(n=N, mean=mus, sd=stds))
      X[[m]]$distr <- "G"
      G <- G + 1
    }else if (identical(distr[m],"P")){
      X[[m]]$x <- as.matrix(stats::rpois(n=N, lambda=lambdas))
      X[[m]]$distr <- "P"
      P <- P + 1
    }else if (identical(distr[m],"B")){
      X[[m]]$x <- as.matrix(stats::rbinom(n=N, size=r[,B], prob=bp))
      X[[m]]$distr <- "B"
      B <- B + 1
    }else if (identical(distr[m],"BPR")){
      X[[m]]$x <- list()
      X[[m]]$x[[1]] <- create_bpr_data(N = N, w = ws, basis = basis)
      X[[m]]$distr <- "BPR"
      BPR <- BPR + 1
    }
  }
  return(list(X=X, L_true=L_true))
}


#' Generate data from the Binomial Probit Regression (BPR) model.
#' The BPRMeth package (https://github.com/andreaskapou/BPRMeth) needs to be 
#' installed to generate these data.
create_bpr_data <- function(N = 300, w, basis, max_L = 35, xmin = -100, 
                            xmax=100){
  require(BPRMeth)
  # Create a list to store data for each methylation region
  X       <- list()
  
  # For each of the N objects
  for (i in 1:N){
    # L is the number of CpGs found in the ith region
    L <- stats::rbinom(n = 1, size = max_L, prob = .8)
    X[[i]] <- matrix(0, nrow = L, ncol = 3)
    # Randomly sample locations for the CpGs
    obs <- sort(sample(xmin:xmax, L))
    # Scale them, so the data lie in the (fmin, fmax) range
    X[[i]][, 1] <- BPRMeth::minmax_scaling(data = obs,
                                           xmin = xmin,
                                           xmax = xmax,
                                           fmin = -1,
                                           fmax = 1)
  }
  
  for (i in 1:N){
    H <- BPRMeth:::.design_matrix(basis, X[[i]][, 1])$H
    p_suc <- stats::pnorm(H %*% w[i, ]) + 
      stats::rnorm(NROW(H), mean = 0, sd = 0.05)
    
    p_suc[which(p_suc > (1 - 1e-10))] <- 1 - 1e-10
    p_suc[which(p_suc < 1e-10)] <- 1e-10
    
    # Total number of trials
    X[[i]][, 2] <- stats::rbinom(NROW(H), 25, 0.7)
    # Total number of successes
    X[[i]][, 3] <- stats::rbinom(NROW(H), X[[i]][, 2], p_suc)
  }
  return(X = X)
}
