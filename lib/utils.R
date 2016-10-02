#' Update the posterior probabilities of source-specific clusterings
#' 
#' \code{update_Lm} updates the posterior probabilities of source-specific 
#' clusterings Lm, that is it computes the conditional distribution of Lm, i.e. 
#' \deqn{L_{m}|X_{m},\theta,C,a \sim v(k, C_{n},a_{m}) * 
#' f_{m}(X_{mn}|\theta_{mk})} where v() is the dependence function and f() is 
#' the observation model (e.g. Normal, Poisson, etc.).
#' 
#' @param X List containing the M data sources X_{m}
#' @param distr Vector with distribution names for the corresponding sources.
#' @param D Vector with dimensions for each data source.
#' @param nu Source-specific probabilities P(L_{m}|C).
#' @param logF Log likelihood matrix for each data source.
#' @param logL Log likelihood posterior probs for each data source.
#' @param Poisson Poisson object in case X_{m} ~ Poisson.
#' @param Binom Binomial object in case X_{m} ~ Binomial.
#' @param Normal Normal object in case X_{m} ~ Normal.
#' @param PoisLL PoisLL object in case X_{m} ~ Poisson Log Linear.
#' @param BinPR BinPR object in case X_{m} ~ Binomial Probit Regression.
#'   
#' @return The updated posterior probabilities for the source-specific
#'   clusterings.
update_Lm <- function(X, distr, D, nu, logF, logL, Poisson, Binom, Normal, PoisLL, BinPR){
  # Obtain number of sources and clusters from the \nu array dimensions
  M <- dim(nu)[2]
  K <- dim(nu)[3]
  
  if(identical(X[[1]]$distr, "BPR")){
    N <- length(X[[1]]$x[[1]])
  }else{
    N <- NROW(X[[1]]$x)  # Sample size (number of objects)
  }
  G = B = P = PLL = BPR <- 1              # Initialize counts of each source
  for (m in 1:M){                         # For each source
    if (identical(distr[m], "P")){        # If source ~ Poisson
      if (D[m] == 1){                     # If univariate data
        for (k in 1:K)
          logF[[m]][, k] <- log(nu[, m, k]) + stats::dpois(X[[m]]$x, lambda = Poisson$l[[P]][k], 
                                                    log = TRUE)
      }else if (D[m] > 1){                # If multivariate data
        for (k in 1:K)
          logF[[m]][, k] <- log(nu[, m, k]) + rowSums(stats::dpois(X[[m]]$x, 
                                                            lambda = matrix(rep(Poisson$l[[P]][k, ], each = N), 
                                                                            ncol = D[m]), log = TRUE))
      }
      P <- P + 1
    }else if (identical(distr[m], "B")){  # If source ~ Binomial
      if (D[m] == 1){                     # If univariate data
        for (k in 1:K)
          logF[[m]][, k] <- log(nu[, m, k]) + stats::dbinom(X[[m]]$x, size = Binom$r[[B]], 
                                                     prob = Binom$p[[B]][k], log = TRUE)
      }else if (D[m] > 1){                # If multivariate data
        for (k in 1:K)
          logF[[m]][,k] <- log(nu[, m, k]) + rowSums(stats::dbinom(X[[m]]$x, size = Binom$r[[B]], 
                                                            prob = matrix(rep(Binom$p[[B]][k, ], each = N), 
                                                                          ncol = D[m]), log = TRUE))
      }
      B <- B + 1
    }else if (identical(distr[m], "G")){  # If source ~ Normal
      if (D[m] == 1){                     # If univariate data
        for (k in 1:K)
          logF[[m]][, k] <- log(nu[, m, k]) + stats::dnorm(X[[m]]$x, mean = Normal$mu[[G]][[k]], 
                                                    sd = 1 / sqrt(Normal$Tau[[G]][[k]]), log = TRUE)
      }else if (D[m] > 1){                # If multivariate data
        for (k in 1:K)
          logF[[m]][, k] <- log(nu[, m, k]) + mvtnorm::dmvnorm(X[[m]]$x, mean = Normal$mu[[G]][[k]], 
                                                               sigma = diag(1 / Normal$Tau[[G]][[k]]), 
                                                               log = TRUE)
      }
      G <- G + 1
    }else if (identical(distr[m], "PLL")){ # If source ~ Poisson Log Linear
      # Compute mean matrix using the estimated lambdas, normalization factors s and
      # the overall expression levels for each object w.
      mean.mat <- vector("list", K)        # List for holding the mean matrices l
      for (k in 1:K){
        lambda.mat <- matrix(rep(rep(PoisLL$l[[PLL]][, k], times = PoisLL$r[[PLL]]), each = N), 
                             nrow = N, ncol = PoisLL$D[[PLL]])
        mean.mat[[k]] <- PoisLL$w.mat[[PLL]] * PoisLL$s.mat[[PLL]] * lambda.mat
        
        logF[[m]][, k] <- log(nu[, m, k]) + rowSums(stats::dpois(X[[m]]$x, lambda = mean.mat[[k]], log = TRUE))
      }
      PLL <- PLL + 1
    }else if (identical(distr[m], "BPR")){
      for (k in 1:K){
        for (d in 1:D[m]){
          logF[[m]][, k] <- logF[[m]][, k] + vapply(X = 1:N,
                            FUN = function(y)
                              bpr_likelihood(w = BinPR$w[[BPR]][[d]][, k],
                                             H = BinPR$des_mat[[BPR]][[d]][[y]],
                                             data = X[[m]]$x[[d]][[y]][, 2:3],
                                             is_NLL = FALSE),
                            FUN.VALUE = numeric(1),
                            USE.NAMES = FALSE)
        }
        logF[[m]][, k] <- log(nu[, m, k]) + logF[[m]][, k]
      }
      BPR <- BPR + 1
    }      
    
    # Normalize log probability by using the logSumExp trick
    logL[[m]] <- logF[[m]] - apply(logF[[m]], 1, log_sum_exp)
    #print(logL)
  }
  
  # Exponentiate for normalized posterior source probs
  L_posterior = lapply(logL, exp)
  return(L_posterior)
}


#' Update the adherence parameter \alpha
#' 
#' \code{update_alpha} updates the adherence parameter \alpha, i.e. adherence of
#' data source m to the overall clustering C. We update it using a Truncated 
#' Beta distribution, i.e.: \deqn{a_{m}|C, L_{m} ~ TBeta(a+t_{m}, b+N-t_{m},
#' 1/K)}
#' 
#' @param C Overall clustering assignments.
#' @param L_m Source-specific clustering assignments.
#' @param indiv_alpha Equal or not equal adherence parameter alpha.
#' @param aBeta \alpha parameter of the TBeta distribution.
#' @param bBeta \beta parameter of the TBeta distribution.
#' 
#' @return Updated adherence parameter.
#' 
update_alpha <- function(C, L_m, indiv_alpha, aBeta, bBeta){
  # Obtain total number of sources, objects and clusters
  M <- length(L_m)
  N <- NROW(C)
  K <- NCOL(C)
  
  # Initialize \alpha
  alpha <- rep(1 / K, M)
  
  if (indiv_alpha){   # Assume not equal adherence to overall clustering
    for (m in 1:M){
      # Number of samples n satisfying L_mn == C_n
      num_eq <- sum(L_m[[m]] == 1 & C == 1)
      
      for(Count in 1:10){ # Try 10 times to sample from Truncated beta
        alpha_temp <- stats::rbeta(1, shape1 = aBeta + num_eq, 
                                   shape2 = bBeta + N - num_eq)
        if (alpha_temp > 1/K){
          # Update the adherence parameter for source m
          alpha[m] <- alpha_temp 
          break
        }
      }
    }
  }else{              # Assume equal adherence
    num_eq <- 0
    for (m in 1:M){   # Number of samples n satisfying L_mn == C_n
      num_eq <- num_eq + sum(L_m[[m]] == 1 & C == 1) 
    }
    for(Count in 1:10){ # Try 10 times to sample from Truncated beta
      alpha_temp <- stats::rbeta(1, shape1 = aBeta + num_eq, 
                                 shape2 = bBeta + N * M - num_eq)
      if (alpha_temp > 1/K){
        # Update the adherence parameter
        alpha[] <- alpha_temp 
        break
      }
    }
  }
  return(alpha)
}


#' Update posterior probabilities for overall clustering
#' 
#' \code{update_post_prob} updates posterior probabilities for the overall 
#' clustering C, i.e. compute the conditional distribution of C:
#' \deqn{C|L_{m}, P, a \sim \pi_{k} * \prod_{m=1}^{M} v(k, L_{mn},a_{m})}
#' where v() is the dependence function.
#' 
#' @param post_prob Overall clustering posterior probabilities (i.e.
#'   responsibilities).
#' @param L_m Source-specific clustering assignments.
#' @param pi_k Mixing proportions for overall clustering.
#' @param alpha Adherence parameter of data source m to overall clustering C.
#' 
#' @return Updated posterior probabilities
#' 
update_post_prob <- function(post_prob, L_m, pi_k, alpha){
  
  # Extract the total number of sources and clusters
  M <- length(alpha)
  K <- length(pi_k)
  
  for (k in 1:K){
    post_prob[, k] <- pi_k[k]
    # For each source ... 
    for (m in 1:M){
      # ... compute the posterior probabilities
      post_prob[,k]  <- post_prob[, k] * alpha[m] ^ (L_m[[m]][, k]) *
        ((1 - alpha[m]) / (K - 1)) ^ (1 - L_m[[m]][, k])
    }
  }
  # Obtain posterior probabilities by normalization
  post_prob  <- post_prob / rowSums(post_prob)
  return(post_prob)
}


#' Draw overall mixture components C
#' 
#' \code{update_C} draws the overall mixture components C using a Discrete 
#' distribution, i.e. multionomial from one sample.
#' 
#' @param post_prob Overall clustering posterior probabilities (i.e.
#'   responsibilities).
#' @param C Overall clustering assignments.
#' 
#' @return The overall clustering draws (i.e. assignments).
#' 
update_C <- function(post_prob, C){
  # Obtain number of objects
  N <- NROW(post_prob)
  for (i in 1:N){
    # Draw overall mixture components
    C[i, ] = stats::rmultinom(1, 1, post_prob[i, ])
  }
  return(C)
}


#' Update mixing proportions \pi
#' 
#' \code{update_pi} updates the mixing proportions \pi of overall clustering C, 
#' that is the probability that an object belongs to the overall cluster k. We 
#' update the mixing proportions by sampling from a Dirichlet distribution.
#' 
#' @param a_dirichlet Concentration parameter for the Dirichlet distribution.
#' @param C Overall clustering assignments.
#' 
#' @return Updates mixing proportions.
#' 
update_pi <- function(a_dirichlet, C){
  # Sample the mixing proportions from Dirichlet
  pi_k  <- MCMCpack::rdirichlet(n = 1, alpha = a_dirichlet + colSums(C))
  return(pi_k)
}


#' Update dependence function \nu
#' 
#' \code{update_nu} updates the dependence function \nu, which essentially 
#' needed for computing the source-specific clustering probabilities, i.e. 
#' \deqn{L_{m}|C ~ v(k, C_{n}, a_{m})}
#' 
#' @param nu Source-specific probabilities P(L_{m}|C).
#' @param C Overall clustering assignments.
#' @param alpha Adherence parameter of data source m to overall clustering C.
#'   
#' @return Updated source-specific probabilities.
#'   
update_nu <- function(nu, C, alpha){
  # Obtain number of sources and clusters from the \nu array dimensions
  M <- dim(nu)[2]
  K <- dim(nu)[3]
  # For each source m ...
  for(m in 1:M){
    # For each cluster k ...
    for(k in 1:K){
      # Set the rows to alpha when C = 1
      nu[as.logical(C[, k]), m, k] <- alpha[m]
      
      # Set the rows to (1-alpha)/(K-1) when C = 0
      nu[!as.logical(C[, k]), m, k] <- (1 - alpha[m]) / (K - 1)
    }
  }
  return(nu)
}


#' Compute the posterior of the adherence parameter.
#' 
#' \code{alpha_posterior} computes the posterior of the adherence parameter to 
#' the overall clustering by taking the average over the MCMC simulations and 
#' also computing the 95% credible intervals.
#' 
#' @param alpha_vec Adherence parameter value on each MCMC simulation.
#' @param M Number of data sources.
#' @param indiv_alpha Equal or not equal adherence parameter alpha.
#' @param mcmc_sims Number of MCMC simulations.
#'   
#' @return List with two elements, containing the posterior mean of alpha and
#'   the 95% credible intervals.
#'   
compute_alpha_posterior <- function(alpha_vec, M, indiv_alpha, mcmc_sims){
  # If not equal adherence to overall clustering
  if(indiv_alpha){
    alpha_cred_interv <- matrix(nrow = M, ncol = 2)
    for(m in 1:M){
      # Compute credible intervals
      alpha_cred_interv[m, 1] <- stats::quantile(alpha_vec[m, floor(mcmc_sims / 5):mcmc_sims], 0.025)
      alpha_cred_interv[m, 2] <- stats::quantile(alpha_vec[m, floor(mcmc_sims / 5):mcmc_sims], 0.975)
    }
    # Compute posterior mean
    alpha_post_mean <- rowMeans(alpha_vec[, floor(mcmc_sims / 5):mcmc_sims])
  }else{
    # Compute credible intervals
    alpha_cred_interv <- c()
    alpha_cred_interv[1] <- stats::quantile(alpha_vec[floor(mcmc_sims / 5):mcmc_sims], 0.025)
    alpha_cred_interv[2] <- stats::quantile(alpha_vec[floor(mcmc_sims / 5):mcmc_sims], 0.975)
    # Compute posterior mean
    alpha_post_mean <- mean(alpha_vec[floor(mcmc_sims / 5):mcmc_sims])
  }
  return(list(alpha_post_mean = alpha_post_mean, alpha_cred_interv = alpha_cred_interv))
}


#' Compute the library size normalization factors
#' 
#' \code{norm_factors} function computes the library size normalization factors 
#' using 4 different methods. These normalization factors 's' are estimated from
#' the data and are considered to be fixed parameters in the Poisson log linear 
#' mixture model.
#' 
#' Possible methods: 1. TC : simply normalize by total counts 2. MED : Compute
#' the median of the data and normalize 3. DESEQ : Median ration normalization
#' developed by DESeq package 4. TMM : Trimmed Mean of M-values normalization
#' used in edgeR package.
#' 
#' @param X Observation data for PoisLL model.
#' @param is_lib_size Logical, should we perform library size normalization.
#' @param lib_type Character, denoting the method for performing library size
#'   normalization.
#'
#' @return The library size normalization factor
norm_factors <- function(X, is_lib_size = TRUE, lib_type = "DESEQ"){
  # Number of variables
  q <- NCOL(X)
  # Normalized library size for each variable
  s <- rep(1, q)
  
  if (is_lib_size == TRUE) {
    if (lib_type == "TC"){              # 'TC' case
      s <- colSums(X) / sum(X)
    }else if (lib_type == "MED"){       # 'MED' case
      s <- apply(X, MARGIN = 2, FUN = median) / 
        sum(apply(X, MARGIN = 2, FUN = median))
    }else if (lib_type == "DESEQ"){     # 'DESEQ' case
      ## Code from DESeq function 'estimateSizeFactorsForMatrix'
      loggeomeans <- rowMeans(log(X))
      s <- apply(X, MARGIN = 2, FUN = function(x) 
        exp(median((log(x) - loggeomeans)[is.finite(loggeomeans)])))
      s <- s / sum(s)
    }else if (lib_type == "TMM"){       # 'TMM' case
      require(edgeR)
      f <- edgeR::calcNormFactors(X, method = "TMM")
      s <- colSums(X) * f / sum(colSums(X) * f)
    }
  }
  return(s)
}


#' Calculate lambdas for the Poisson Log Linear model
#' 
#' \code{calculate_lambdas} calculate the unknown parameters lambda of the 
#' Poisson Log Linear model, which correspond to the clustering parameters that 
#' define the profiles of the gene in cluster k across biological condition d. 
#' Thus, the dimension of the matrix is d x K, and is given by the following 
#' equation: \deqn{l_{jk} = (\sum_{i}w_{ik}X_{ij.}) /
#' (s_{j.}\sum_{i}t_{ik}X_{i..})}
#' 
#' @param X Observation data for PoisLL model.
#' @param w Overall expression for each object.
#' @param d Number of biological conditions.
#' @param s_dot Sum of s for all replicates on each condition d.
#' @param conds Different biological conditions.
#' @param post_prob Matrix of (hard) cluster assignments.
#' 
#' @return Calculated parameters lambda
#' 
calculate_lambdas <- function(X, w, d, s_dot, conds, post_prob){
  # Obtain total number of clusters
  K <- NCOL(post_prob)
  
  lambdas <- matrix(0, nrow = d, ncol = K)
  # Calculate RHS of denominator
  rhs_den <- colSums(post_prob * w)  
  for(j in 1:d) {
    # Multiply with LHS to get denominator value
    denom <- s_dot[j] * rhs_den         
    # Compute \sum_{i}\sum_{l} X_{ijl} = X_j
    X_j <- rowSums(as.matrix(X[, which(conds == (unique(conds))[j])]))
    # Compute the numerator
    num <- colSums(post_prob * matrix(rep(X_j, K), ncol = K))
    #Get the parameter value for condition d
    lambdas[j, ] <- num / denom
  }
  return(lambdas)
}


#' Compute stable Log-Sum-Exp
#'
#' \code{.log_sum_exp} computes the log sum exp trick for avoiding numeric
#' underflow and have numeric stability in computations of small numbers.
#'
#' @param x A vector of observations
#'
#' @return The logs-sum-exp value
#'
#' @references
#'  \url{https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/}
#'
log_sum_exp <- function(x) {
  # Computes log(sum(exp(x))
  offset <- max(x)
  return(log(sum(exp(x - offset))) + offset)
}


#' Compute source-specific clustering error
#' 
#' \code{compute_source_error} computes the source-specific clustering
#' assignment error, i.e. the average number of incorrect cluster assignments:
#' \deqn{SE = \sum_{m=1}^{M}(\sum_{n=1}^{N}(I(LT_{mn} \neq LP_{mn}))/M*N )}
#' 
#' @param L_true True source-specific cluster assignemnts.
#' @param L_post Posterior mean of predicted source-specific cluster 
#'   assignemnts.
#'   
#' @return The source-specific clustering assignment error
#'   
compute_source_error <- function(L_true, L_post){
  # Obtain the total number of sources and objects
  M <- length(L_post)
  N <- NROW(L_post[[1]])
  
  L_match <- 0
  for (m in 1:M){
    # Convert probs to hard clusterings
    L_post[[m]][L_post[[m]] < .5] <- 0 
    L_post[[m]][L_post[[m]] >= .5] <- 1
    
    # Align the cluster indices
    L_post[[m]]  <- align_clusters(L_true[[m]], L_post[[m]], type = "mat")
    # Count the number of correct assignments
    L_match     <- L_match + sum((L_post[[m]] == L_true[[m]])[, 1])
  }
  # Compute the source error
  source_error <- 1 - L_match / (M * N)
  return(source_error)
}


#' Compute overall clustering error
#' 
#' \code{compute_overall_error} computes the overall clustering assignment 
#' error, i.e. the average number of incorrect cluster assignments: \deqn{OE = 
#' \sum_{n=1}^{N}(I(LT_{mn} \neq LP_{mn})) / N}
#' 
#' @param C_true True overall cluster assignemnts.
#' @param C_post Posterior mean of predicted overall cluster assignemnts.
#'   
#' @return The overall clustering assignment error
#'   
compute_overall_error <- function(C_true, C_post){
  # Obtain the total number of objects
  N <- NROW(C_post)
  
  # Convert probs to hard clusterings
  C_post[C_post < .5] <- 0
  C_post[C_post >= .5] <- 1
  
  # Align cluster indices
  C_post   <- align_clusters(C_true, C_post, type = "mat")
  # Count number of correct assignments
  C_match <- sum((C_post == C_true)[, 1])
  
  # Compute the overall error
  overall_error <- 1 - C_match / N
  return(overall_error)
}


# Internal function to make all the appropriate type checks for BPR model.
do_bpr_checks <- function(x, K = 2, w = NULL, basis = NULL,
                          opt_method = "CG", opt_itnmax = 100,
                          is_parallel = TRUE, no_cores = NULL){
  # Number of dimensions
  D <- length(x)
  if (length(basis) == 0){
    basis <- BPRMeth::create_rbf_object(M = 3)
  }
  if (length(w) == 0){
    w <- list()
    for (d in 1:D){
      w_d <- rep(0.5, basis$M + 1)
      # Optimize the BPR function for each element in x
      out_opt <- BPRMeth::bpr_optim(x           = x[[d]],
                                    w           = w_d,
                                    basis       = basis,
                                    fit_feature = NULL,
                                    cpg_dens_feat = FALSE,
                                    method      = opt_method,
                                    itnmax      = opt_itnmax,
                                    is_parallel = is_parallel,
                                    no_cores    = no_cores)
      
      # Keep only the optimized coefficients
      W_opt <- out_opt$W_opt
      
      # Use Kmeans with random starts
      cl <- stats::kmeans(W_opt, K, nstart = 25)
      # Get the mixture components
      C.n <- cl$cluster
      # Mean for each cluster
      w[[d]] <- t(cl$centers)
    }
  }
  if (NROW(w[[1]]) != (basis$M + 1) ){
    stop("Coefficients vector should be M+1, M: number of basis functions!")
  }
  return(list(w = w, basis = basis, C.n = C.n))
}