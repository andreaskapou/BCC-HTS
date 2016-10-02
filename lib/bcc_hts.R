#' Run the BCC model for HTS data
#' 
#' \code{bcc_hts} implements the BCC model propose by Lock & Dunson (2013, 
#' Bioinformatics) and extends it with additional observation likelihoods, such 
#' as Binomial, Poisson and BPR model (https://github.com/andreaskapou).
#' 
#' @param X List containing the M data sources X_{m}. Each element in the list 
#'   contains the data source M and the name of the observation model that 
#'   source M follows. Currently, 'G' for Normal, 'B' for Binomial, 'P' for 
#'   Poisson, and 'BPR' for the BPR model. You can access the names as 
#'   \code{X[[ind]]$distr}.
#' @param K The total number of clusters. It is the same for source specific and
#'   overall clusterings.
#' @param pi_k Mixing proportions for overall clustering.
#' @param indiv_alpha Equal or not equal adherence parameter alpha.
#' @param hyper A list of the BCC model hyperparameters.
#' @param params A list of size M, containing the parameters of each data 
#'   source, i.e. initial values for the parameters of the prior, etc.
#' @param N_sims Number of MCMC simulations.
#' @param burnin Burn-in period, i.e how many of the initial draws we should 
#'   discard.
#' @param is_stephens Logical, whether to store all the MCMC draws, in order to 
#'   later use the Stephens relabelling algorithm for the label switching 
#'   problem. Currently, not implemented inside this function.
#' @param is_hard_cluster Logical, whether to produce a single (hard) clustering
#'   as proposed by Dahl (2006).
#' @param eq_mix_prop Logical, whether to impose constraint of equal mixing
#'   proportions.
#' @param is_parallel Logical, indicating if the code should be run in parallel.
#' @param no_cores Number of cores to be used.
#' 
#' @return A bcc object.
#'   
bcc_hts <- function(X, K, pi_k = NULL, indiv_alpha = FALSE, hyper = NA, 
                    params = NA, N_sims = 10000, burnin = 2000, 
                    is_stephens = FALSE, is_hard_cluster = FALSE, 
                    eq_mix_prop = FALSE, is_parallel = TRUE, no_cores = 3){
  
  ##------------------------------
  # Require libraries            #
  ##------------------------------
  require(MCMCpack)
  require(BPRMeth)
  require(mvtnorm)
  require(truncnorm)
  require(reshape2)
  
  ##------------------------------
  # Initialize variables         #
  ##------------------------------
  Normal  <- list()                 # Create a Normal object
  Binom   <- list()                 # Create a Binomial object
  Poisson <- list()                 # Create a Poisson object
  PoisLL  <- list()                 # Create a Poisson Log-Linear object
  BinPR   <- list()                 # Create a Binomial Probit Regression object
  
  if ("aBeta" %in% names(hyper)){   # alpha parameter of the TBeta distribution
    aBeta <- hyper$aBeta
  }else{
    aBeta <- 1
  }
  
  if ("bBeta" %in% names(hyper)){   # beta parameter of the TBeta distribution
    bBeta <- hyper$bBeta
  }else{
    bBeta <- 1
  }
  
  if ("aDir" %in% names(hyper)){    # Concentration parameter for Dirichlet 
    a_dirichlet <- rep(hyper$aDir, K)
  }else{
    a_dirichlet <- rep(1/K, K)
  }
  
  mcmc_sims <- N_sims - burnin      # Number of MCMC draws after burnin period
  M <- length(X)                    # Number of data sources
  if(identical(X[[1]]$distr,"BPR")){
    N <- length(X[[1]]$x[[1]])      # Sample size (number of objects)
  }else{
    N <- NROW(X[[1]]$x)     
  }
  
  # Vector of dimensions for each data source
  D <- vector(mode = "integer", length = M)
  # Vector of observation models for each data source
  distr <- vector(mode = "character", length = M)
  for (m in 1:M){
    if(identical(X[[m]]$distr,"BPR")){
      D[m] <- length(X[[m]]$x)   
    }else{
      D[m] <- NCOL(X[[m]]$x)
    }
    distr[m] <- X[[m]]$distr
  }
  
  L_m         <- list()             # Clustering for each data source
  LPost       <- list()             # Counts L_mn = k for each MCMC draw
  CPost       <- matrix(0, nrow = N, ncol = K) # Counts C_n = k for each MCMC draw
  logF        <- list()             # Log likelihood matrix for each data source
  logL        <- list()             # Log clustering probs for each data source
  if (is_hard_cluster){
    Llist     <- list()             # List of source-specific clusterings
    Clist     <- list()             # List of overall clusterings
  }
  
  if (is.null(pi_k)){
    pi_k <- rep(1/K,K)              # Cluster probabilites (mixing proportions)
  }
  
  # Overall clustering responsibilities
  post_prob <- matrix(0, nrow = N, ncol = K) 
  # Overall mixture components
  C <- matrix(t(stats::rmultinom(N, 1, rep(1/K, K))), nrow = N, ncol = K)
  # Size of each source-specific cluster
  N_m_k <- matrix(0, nrow = M, ncol = K)
  # Source-specific probs P(L_m|C)
  nu <- array(rep(1/K, N * M * K), dim = c(N, M, K)) 
  # Post resp for each MCMC run if Stephens method will be used
  if (is_stephens){
    post_resp_arr <- array(0, dim = c(N_sims - burnin, N, K))
  }
  
  ##-------------------------------------
  # Initialize adherence parameter      #
  ##-------------------------------------
  alpha_vec <- c()            # Vector of alpha (adherence) parameters
  alpha <- rep(0,M)           # Initialize alpha
  if (indiv_alpha){           # Not equal adherence to overall clustering
    for (m in 1:M)
      while(alpha[m] < 1/K)   # Sample from a Truncated Beta
        alpha[m] <- stats::rbeta(1, shape1=aBeta, shape2=bBeta)
  }else{
    while(alpha[1] < 1/K)
      alpha[] <- stats::rbeta(1, shape1=aBeta, shape2=bBeta)
  }
  
  pll.draws   <- list()   # l parameter of Poisson Log-Linear in each MCMC draw
  l.draws     <- list()   # Mean/variance of Poisson in each MCMC draw
  bp.draws    <- list()   # Success probability of Binomial in each MCMC draw
  mu.draws    <- list()   # Mean of Normal in each MCMC draw
  tau.draws   <- list()   # Precision (matrix) of Normal in each MCMC draw
  w.draws     <- list()   # Coefficients of BPR model in each MCMC draw
  pi.draws    <- matrix(0, nrow = N_sims-burnin, ncol = K) # Mixing proportions
  
  ##---------------------------------
  # Initialize Parameters           #
  ##---------------------------------
  G = B = P = PLL = BPR <- 1      # Initialize counts of each source
  
  for (m in 1:M){
    LPost[[m]]    <- matrix(0, nrow = N, ncol = K)
    L_m[[m]]      <- matrix(0, nrow = N, ncol = K)
    logF[[m]]     <- matrix(0, nrow = N, ncol = K)
    logL[[m]]     <- matrix(0, nrow = N, ncol = K)
    if (is_hard_cluster){
      Llist[[m]]    <- list()
    }
    
    if (identical(distr[m],"B")){
      if (!("r" %in% names(params[[m]]))){
        stop("'params' argument should contain value for the total number of trials 'r'")
      }
      Binom$r[[B]] <- as.matrix(params[[m]]$r)  # Total number of trials
      # Initialize probability of success for Binomial
      cl  <- stats::kmeans(X[[m]]$x/Binom$r[[B]], K, nstart=25) # Kmeans with random starts
      C.n <- cl$cluster # Get the mixture components
      
      # Initialize Beta prior hyperparameters
      if ( ("aBeta" %in% names(params[[m]])) && ("bBeta" %in% names(params[[m]])) ){
        Binom$Beta[[B]] <- list(a=params[[m]]$aBeta, b=params[[m]]$bBeta)
      }else{
        Binom$Beta[[B]] <- list(a=1, b=1)
      }
      
      if (D[m] == 1){                    # If univariate data
        Binom$p[[B]]  <- t(cl$centers)   # Success binomial probability
        bp.draws[[B]] <- matrix(0, nrow=N_sims-burnin, ncol=K)
      }else if (D[m] > 1){               # If multivariate data
        Binom$p[[B]]  <- cl$centers      # Success binomial probability
        bp.draws[[B]] <- array(0, dim = c(N_sims-burnin, K, D[m]))
      }
      B <- B + 1  # Increase count of Binomial sources
      
    }else if (identical(distr[m],"BPR")){
      if ("basis" %in% names(params[[m]])){
        BinPR$basis[[BPR]] <- params[[m]]$basis
      }else{
        BinPR$basis[[BPR]] <- list()
      }
      
      if ("w" %in% names(params[[m]])){
        BinPR$w[[BPR]] <- params[[m]]$w
      }else{
        BinPR$w[[BPR]] <- list()
      }
      if ("opt_method" %in% names(params[[m]])){
        BinPR$opt_method[[BPR]] <- params[[m]]$opt_method
      }else{
        BinPR$opt_method[[BPR]] <- "CG"
      }
      if ("opt_itnmax" %in% names(params[[m]])){
        BinPR$opt_itnmax[[BPR]] <- params[[m]]$opt_itnmax
      }else{
        BinPR$opt_itnmax[[BPR]] <- 100
      }
      res <- do_bpr_checks(x = X[[m]]$x, 
                           K = K, 
                           w = BinPR$w[[BPR]], 
                           basis = BinPR$basis[[BPR]],
                           opt_method = BinPR$opt_method[[BPR]], 
                           opt_itnmax = BinPR$opt_itnmax[[BPR]],
                           is_parallel = is_parallel, 
                           no_cores = no_cores)
      BinPR$basis[[BPR]] <- res$basis
      BinPR$w[[BPR]]     <- res$w
      
      # Get the mixture components
      C.n <- res$C.n
      
      # Initialize BPR prior hyperparameters
      if ( ("w_0_mean" %in% names(params[[m]])) && ("w_0_cov" %in% names(params[[m]])) ){
        BinPR$Norm[[BPR]] <- list(mu.0=params[[m]]$w_0_mean, s.0=params[[m]]$w_0_cov)
      }else{
        BinPR$Norm[[BPR]] <- list(mu.0=rep(0, BinPR$basis[[BPR]]$M+1), 
                                  s.0=diag(1, BinPR$basis[[BPR]]$M+1))
      }
      BinPR$prec_0[[BPR]] <- solve(BinPR$Norm[[BPR]]$s.0)
      BinPR$w_0_prec_0[[BPR]] <- BinPR$prec_0[[BPR]] %*% BinPR$Norm[[BPR]]$mu.0
      
      # List to keep w draws
      w.draws[[BPR]] <- list()
      
      BinPR$des_mat[[BPR]] < list()
      for (d in 1:D[m]){
        if (is_parallel){
          # Create design matrix for each observation
          BinPR$des_mat[[BPR]][[d]] <- parallel::mclapply(X = X[[m]]$x[[d]],
                              FUN = function(y)
                                BPRMeth:::.design_matrix(x = BinPR$basis[[BPR]], 
                                              obs = y[ ,1])$H,
                              mc.cores = no_cores)
        }else{
          # Create design matrix for each observation
          BinPR$des_mat[[BPR]][[d]] <- lapply(X = X[[m]]$x[[d]],
                            FUN = function(y)
                              BPRMeth:::.design_matrix(x = BinPR$basis[[BPR]], 
                                            obs = y[ ,1])$H)
        }
        
        ## ----------------------------------------------------------------------
        # Auxiliary variable model parameters
        BinPR$ext_des_mat[[BPR]][[d]] <- list()
        BinPR$data_y[[BPR]][[d]] <- list()
        # N1 in first column, N0 in second column
        BinPR$suc_fail_mat[[BPR]][[d]] <- matrix(NA_integer_, ncol = 2, nrow = N)
        
        # Iterate over each region
        for (i in 1:N){
          # Total number of reads for each CpG
          N_i <- X[[m]]$x[[d]][[i]][, 2]
          # Corresponding number of methylated reads for each CpG
          m_i <- X[[m]]$x[[d]][[i]][, 3]
          
          # Create extended vector y of length (J x 1)
          y <- vector(mode = "integer")
          for (j in 1:length(N_i)){
            y <- c(y, rep(1, m_i[j]), rep(0, N_i[j] - m_i[j]))
          }
          BinPR$data_y[[BPR]][[d]][[i]] <- y
          
          # Col1: Number of successes
          # Col2: Number of failures
          BinPR$suc_fail_mat[[BPR]][[d]][i, ] <- c(sum(y), sum(N_i) - sum(y))
          
          # TODO: Keep only one design matrix POSSIBLE MEMORY ISSUE
          # Create extended design matrix H of dimension (J x M)
          BinPR$ext_des_mat[[BPR]][[d]][[i]] <- 
            as.matrix(BinPR$des_mat[[BPR]][[d]][[i]][
              rep(1:NROW(BinPR$des_mat[[BPR]][[d]][[i]]), N_i), ])
        }
      }
      BPR <- BPR + 1  # Increase count of BPR sources
      
    }else if (identical(distr[m],"PLL")){
      if (!("conds" %in% names(params[[m]]))){
        stop("'params' argument should contain conditions experiments 'conds'")
      }
      
      PoisLL$w[[PLL]] <- rowSums(X[[m]]$x)      # Overall expression for each object
      PoisLL$conds[[PLL]] <- params[[m]]$conds  # Different biological conditions
      PoisLL$d[[PLL]] <- length(unique(params[[m]]$conds))    # Total number of conditions
      PoisLL$r[[PLL]] <- as.vector(table(params[[m]]$conds))  # Num of replicates in each cond
      PoisLL$D[[PLL]] <- D[m]                   # Total number of variables (dimensions)
      
      # Initialize libSize parameter
      if ("libSize" %in% names(params[[m]])){
        PoisLL$libSize[[PLL]] <- params[[m]]$libSize
      }else{
        PoisLL$libSize[[PLL]] <- TRUE
      }
      
      # Initialize libType parameter
      if ("libType" %in% names(params[[m]])){
        PoisLL$libType[[PLL]] <- params[[m]]$libType
      }else{
        PoisLL$libType[[PLL]] <- "DESEQ"
      }
      
      # Compute the library size normalization factors for each variable
      PoisLL$s[[PLL]] <- norm_factors(X = X[[m]]$x,
                                      is_lib_size = PoisLL$libSize[[PLL]],
                                      lib_type = PoisLL$libType[[PLL]])
      
      # Sum of s for all replicates l on each condition d
      PoisLL$s.dot[[PLL]] <- rep(NA, PoisLL$d[[PLL]])
      for (j in 1:PoisLL$d[[PLL]]){
        PoisLL$s.dot[[PLL]][j] <- sum( PoisLL$s[[PLL]][which(PoisLL$conds[[PLL]] == 
                                                               unique(PoisLL$conds[[PLL]])[j])] )
      }
      # Create matrices of dimension N x D, for faster computations
      PoisLL$w.mat[[PLL]] <- matrix(rep(PoisLL$w[[PLL]], times=PoisLL$D[[PLL]]), 
                                    nrow=N, ncol=PoisLL$D[[PLL]])
      PoisLL$s.mat[[PLL]] <- matrix(rep(PoisLL$s[[PLL]], each=N) , nrow=N, ncol=PoisLL$D[[PLL]])
      
      # Use Kmeans with random starts
      cl    <- stats::kmeans(X[[m]]$x/PoisLL$w[[PLL]], K, nstart = 25)
      C.n   <- cl$cluster                       # Get the mixture components
      C.mat <- matrix(0, nrow=N, ncol=K)        # Create matrix of cluster assignments
      for(i in 1:N){
        C.mat[i, C.n[i]] <- 1
      }
      # Call function for calculating the unknown parameters lambda
      PoisLL$l[[PLL]] <- calculate_lambdas(X = X[[m]]$x, 
                                           w = PoisLL$w[[PLL]], 
                                           d = PoisLL$d[[PLL]],
                                           s_dot = PoisLL$s.dot[[PLL]], 
                                           conds = PoisLL$conds[[PLL]], 
                                           post_prob = C.mat)
      
      # Initialize Gamma prior hyperparameters
      if ( ("aGamma" %in% names(params[[m]])) && ("bGamma" %in% names(params[[m]])) ){
        PoisLL$Gamma[[PLL]] <- list(a=params[[m]]$aGamma, b=params[[m]]$bGamma)
      }else{
        PoisLL$Gamma[[PLL]] <- list(a=1, b=1)
      }
      
      if (PoisLL$d[[PLL]] > 1){
        pll.draws[[PLL]] <- array(0, dim=c(N_sims-burnin, K, PoisLL$d[[PLL]]))
      }
      PLL <- PLL + 1 # Increase count of PoisLL sources
      
    }else if (identical(distr[m],"P")){
      # Initialize means\variances for Poisson
      cl      <- kmeans(X[[m]]$x, K, nstart=25) # Kmeans with random starts
      C.n     <- cl$cluster                   # Get the mixture components
      
      # Initialize Gamma prior hyperparameters
      if ( ("aGamma" %in% names(params[[m]])) && ("bGamma" %in% names(params[[m]])) ){
        Poisson$Gamma[[P]]  <- list(a=params[[m]]$aGamma, b=params[[m]]$bGamma)
      }
      else{
        Poisson$Gamma[[P]]  <- list(a=1, b=1) # Initialize Gamma hyperparameters
      }
      
      if (D[m] == 1){                         # If univariate data
        Poisson$l[[P]]  <- t(cl$centers)      # Poisson mean for each cluster
        l.draws[[P]]    <- matrix(0, nrow=N_sims-burnin, ncol=K)
      }else if (D[m] > 1){                    # If multivariate data
        Poisson$l[[P]]  <- cl$centers         # Poisson mean vector for each cluster
        l.draws[[P]]    <- array(0, dim=c(N_sims-burnin, K, D[m]))
      }
      P <- P + 1  # Increase count of Poisson sources
      
    }else if (identical(distr[m],"G")){
      # Kmeans with random starts
      cl  <- stats::kmeans(X[[m]]$x, K, nstart=25)
      C.n <- cl$cluster # Get the mixture components
      
      # Align Cluster indices since they are used later only for Gaussian data
      if (m==1){
        prevC.n <- C.n
      }else{ # Step to align cluster indices
        C.n     <- align_clusters(prevC.n, C.n) 
        prevC.n <- C.n
      }
      
      Normal$mu[[G]]  <- list()
      Normal$Tau[[G]] <- list()
      if (D[m] == 1){             # If univariate data
        mu.draws[[G]] <- matrix(0, nrow=N_sims - burnin, ncol = K)
        tau.draws[[G]] <- matrix(0, nrow=N_sims - burnin, ncol = K)
        for (k in 1:K){
          Normal$mu[[G]][[k]] <- cl$centers[k,] # Normal mean vector
          # Normal precision vector
          Normal$Tau[[G]][[k]] <- 1/var(X[[m]]$x[C.n == k])
        }
      }else if (D[m] > 1){                    # If Multivariate data
        mu.draws[[G]] <- list()
        tau.draws[[G]] <- list()
        for (k in 1:K){
          Normal$mu[[G]][[k]] <- cl$centers[k,] # Normal mean vector
          # Normal precision matrix
          Normal$Tau[[G]][[k]] <- 1/(apply(X[[m]]$x[C.n==k,],MARGIN=2,FUN='sd')^2)
        }
      }
      
      # Initialize Normal and Gamma hyperparameters
      if ( ("mu0" %in% names(params[[m]])) && ("tau0" %in% names(params[[m]])) ){
        Normal$Norm[[G]] <- list(mu.0=params[[m]]$mu0, tau.0=params[[m]]$tau0)
      }else{
        Normal$Norm[[G]] <- list(mu.0=colMeans(X[[m]]$x), tau.0=rep(1/100,D[m]))
      }
      
      if ( ("aGamma" %in% names(params[[m]])) && ("bGamma" %in% names(params[[m]])) ){
        Normal$Gamma[[G]] <- list(b=params[[m]]$aGamma, a=params[[m]]$bGamma)
      }else{
        Normal$Gamma[[G]] <- list(a=rep(1,D[m]), b=apply(X[[m]]$x, 2, sd)^2)
      }
      G <- G + 1
    }
    
    if (m==1){
      prevC.n <- C.n
    }else{    # Step to align cluster indices
      C.n     <- align_clusters(prevC.n, C.n) 
      prevC.n <- C.n
    }
  }
  
  ## ------------------------------------------------------------------------
  
  ## ------------------------------------------------------------------------
  
  cat("Starting MCMC simulations...\n")
  # Show progress bar
  pb <- txtProgressBar(min = 1, max = N_sims, style = 3)
  
  ##--------------------------------------------
  # Start the MCMC simulations                 #
  ##--------------------------------------------
  for (t in 1:N_sims){
    ##------------------------------------------------------------
    # Update posterior source-specific clustering probabilities  #
    ##------------------------------------------------------------
    Lp <- update_Lm(X = X, distr = distr, D = D, nu = nu, logF = logF, 
                    logL = logL, Poisson = Poisson, Binom = Binom, 
                    Normal = Normal, PoisLL = PoisLL, BinPR = BinPR)
    
    G = B = P = PLL = BPR <- 1  # Initialize counts of each source
    for (m in 1:M){
      ##-------------------------------------------
      # Generate source specific assignments L_m  #
      ##-------------------------------------------
      for(i in 1:N){
        L_m[[m]][i,] <- rmultinom(1, 1, Lp[[m]][i, ]) # Generate L_m from Lp
      }
      if (t > 1){ # Align cluster indices
        L_m[[m]] = align_clusters(C, L_m[[m]], type="mat") 
      }
      N_m_k[m,] <- colSums(L_m[[m]]) # Cluster component counts for source m
      
      ##-------------------------------------------
      #       Update Binomial parameters          #
      ##-------------------------------------------
      if (identical(distr[m],"B")){
        for (k in 1:K){
          if (N_m_k[m,k] == 0){ # If cluster empty, set posterior = prior
            alpha.n   <- Binom$Beta[[B]]$a
            beta.n    <- Binom$Beta[[B]]$b
          }else{
            if (D[m] == 1){
              # Calculate sample sum for each cluster for source m
              x.k.sum   <- sum(X[[m]]$x[L_m[[m]][,k]==1])
              # Calculate sample sum of differences for each cluster for source m
              x.k.diff  <- sum(Binom$r[[B]][L_m[[m]][,k]==1]-X[[m]]$x[L_m[[m]][,k]==1])
            }else if (D[m] > 1){
              # Calculate sample sum for each cluster for source m
              x.k.sum   <- colSums(X[[m]]$x[L_m[[m]][,k] == 1,])
              # Calculate sample sum of differences for each cluster for source m
              x.k.diff  <- colSums(Binom$r[[B]][L_m[[m]][,k]==1,]-X[[m]]$x[L_m[[m]][,k]==1,])              
            }
            # Update parameters
            alpha.n   <- Binom$Beta[[B]]$a + x.k.sum
            beta.n    <- Binom$Beta[[B]]$b + x.k.diff
          }
          #Sample from a Beta for the posterior probability parameter
          if (D[m] == 1)
            Binom$p[[B]][k]   <- rbeta(D[m], shape1=alpha.n, shape2=beta.n)
          else if (D[m] > 1)
            Binom$p[[B]][k,]  <- rbeta(D[m], shape1=alpha.n, shape2=beta.n)
        }
        B <- B + 1
      }
      
      ##-------------------------------------------
      #       Update Poisson parameters           #
      ##-------------------------------------------
      else if (identical(distr[m],"P")){
        for (k in 1:K){ 
          if (N_m_k[m,k] == 0){ # If cluster empty, set posterior = prior
            alpha.n   <- Poisson$Gamma[[P]]$a
            beta.n    <- Poisson$Gamma[[P]]$b
          }else{
            if (D[m] == 1){
              # Calculate sample sum for each cluster for source m
              x.k.sum   <-  sum(X[[m]]$x[L_m[[m]][,k] == 1])
              # Update beta parameter
              beta.n    <- Poisson$Gamma[[P]]$b + N_m_k[m,k]
            }else if (D[m] > 1){
              # Calculate sample sum for each cluster for source m
              x.k.sum   <- colSums(X[[m]]$x[L_m[[m]][,k] == 1,])
              # Update beta parameter
              beta.n    <- rep((Poisson$Gamma[[P]]$b + N_m_k[m,k]), D[m])
            }
            # Update alpha parameter
            alpha.n   <- Poisson$Gamma[[P]]$a + x.k.sum
          }
          #Sample from a Gamma for the posterior mean parameter
          if (D[m] == 1)
            Poisson$l[[P]][k]   <- rgamma(D[m], shape=alpha.n, rate=beta.n)
          else if (D[m] > 1)
            Poisson$l[[P]][k,]  <- rgamma(D[m], shape=alpha.n, rate=beta.n)
        }
        P <- P + 1
      }
      
      ##-------------------------------------------------
      #       Update Poisson Log Linear parameters      #
      ##-------------------------------------------------
      else if (identical(distr[m],"PLL")){
        X.k.sum     <- colSums(L_m[[m]] * PoisLL$w[[PLL]])  # Calculate sum of data points
        for (j in 1:PoisLL$d[[PLL]]){
          beta.n    <- PoisLL$Gamma[[PLL]]$b + PoisLL$s.dot[[PLL]][j]*X.k.sum
          X.j.      <- rowSums(as.matrix(X[[m]]$x[,which(PoisLL$conds[[PLL]] == 
                                                           (unique(PoisLL$conds[[PLL]]))[j])]))
          alpha.n   <- PoisLL$Gamma[[PLL]]$a + colSums(L_m[[m]] * matrix(rep(X.j., K), ncol=K))
          
          PoisLL$l[[PLL]][j,] <- rgamma(K, shape=alpha.n, rate=beta.n) # Sample from Gamma
        }
        PLL <- PLL + 1
      }
      
      ##--------------------------------------------
      #       Update Normal/Gaussian parameters    #
      ##--------------------------------------------
      else if (identical(distr[m],"G")){
        for (k in 1:K){ 
          if (N_m_k[m,k] == 0){ # If cluster empty, set posterior = prior
            alpha.n   <- Normal$Gamma[[G]]$a
            beta.n    <- Normal$Gamma[[G]]$b
            mu.n      <- Normal$Norm[[G]]$mu.0
            Tau.n     <- Normal$Norm[[G]]$tau.0
            # Sample from a Gamma for the posterior precision parameter
            Normal$Tau[[G]][[k]]  <- rgamma(D[m], alpha.n, beta.n)
          }else{
            if (D[m] == 1){
              x.k.bar <- mean(X[[m]]$x[L_m[[m]][,k] == 1])
              ssd.k   <- sum((X[[m]]$x[L_m[[m]][,k] == 1] - x.k.bar)^2)
            }else if (D[m] > 1){
              x.k.bar <- colMeans(X[[m]]$x[L_m[[m]][,k] == 1,])
              ssd.k   <- rowSums(apply(X[[m]]$x[L_m[[m]][,k]==1,],1,'-',x.k.bar)^2)
            }
            # Update Gamma hyperparameters
            alpha.n   <- Normal$Gamma[[G]]$a + N_m_k[m,k]/2
            beta.n    <- Normal$Gamma[[G]]$b + ssd.k/2 + ( 
              (N_m_k[m,k] * Normal$Norm[[G]]$tau.0 * 
                 ((x.k.bar - Normal$Norm[[G]]$mu.0)^2)) / 
                (2*(N_m_k[m,k] + Normal$Norm[[G]]$tau.0)) )
            # Sample from a Gamma for the posterior precision parameter
            Normal$Tau[[G]][[k]]  <- rgamma(D[m], alpha.n, beta.n)
            # Update Normal hyperparameters
            t.n.k     <- Normal$Norm[[G]]$tau.0 + N_m_k[m,k]
            mu.n      <- (Normal$Norm[[G]]$mu.0*Normal$Norm[[G]]$tau.0 + 
                            N_m_k[m,k] * x.k.bar ) / t.n.k
            Tau.n     <- Normal$Tau[[G]][[k]] * t.n.k
          }
          #Sample from a Normal for the posterior mean parameter
          Normal$mu[[G]][[k]] <- rnorm(D[m], mean=mu.n, sd=1/sqrt(Tau.n))
        }
        G <- G + 1
      }
      
      ##--------------------------------------
      #       Update BPR parameters          #
      ##--------------------------------------
      else if (identical(distr[m],"BPR")){
        for (k in 1:K){
          # Which regions are assigned to cluster k
          L_M_k_idx <- which(L_m[[m]][, k] == 1)
          for (d in 1:D[m]){
            # Concatenate data from all regions in cluster k
            H <- do.call(rbind, BinPR$ext_des_mat[[BPR]][[d]][L_M_k_idx])
            
            # Concatenate y from all regions in cluster k
            y <- do.call(c, BinPR$data_y[[BPR]][[d]][L_M_k_idx])
            
            # Add all successes and failures from all regions in cluster k
            N1_N0 <- colSums(BinPR$suc_fail_mat[[BPR]][[d]][L_M_k_idx, ])
            
            # Compute posterior variance of w
            V <- solve(BinPR$prec_0[[BPR]] + crossprod(H, H))
            
            # Update Mean of z
            mu_z <- H %*% BinPR$w[[BPR]][[d]][, k]
            # Draw latent variable z from its full conditional: z | \w, y, X
            z <- rep(NA_real_, sum(N1_N0))
            z[y == 1] <- rtruncnorm(N1_N0[1], mean = mu_z[y == 1], sd = 1,
                                    a = 0, b = Inf)
            z[y == 0] <- rtruncnorm(N1_N0[2], mean = mu_z[y == 0], sd = 1,
                                    a = -Inf, b = 0)
            
            # Compute posterior mean of w
            Mu <- V %*% (BinPR$w_0_prec_0[[BPR]] + crossprod(H, z))
            # Draw variable \w from its full conditional: \w | z, X
            BinPR$w[[BPR]][[d]][, k] <- c(rmvnorm(1, Mu, V))
          }
        }
        BPR <- BPR + 1
      }
    }
    ##------------------------------------------
    # Update alpha (i.e. adherence parameter)  #
    ##------------------------------------------
    if (t > 1){
      alpha <- update_alpha(C, L_m, indiv_alpha, aBeta, bBeta)
    }
    
    ##--------------------------------------------
    # Update overall posterior responsibilities  #
    ##--------------------------------------------
    post_prob <- update_post_prob(post_prob, L_m, pi_k, alpha)
    
    ##------------------------------------
    # Draw overall mixture components C  #
    ##------------------------------------
    C <- update_C(post_prob, C)
    
    ##----------------------------------------
    # Update pi_k (i.e. mixing proportions)  #
    ##----------------------------------------
    if (!eq_mix_prop){
      pi_k <- update_pi(a_dirichlet, C)
    }
    
    ##--------------------------------------------------
    # Update source-specific clustering probabilities  #
    ##--------------------------------------------------
    nu <- update_nu(nu, C, alpha)
    
    
    ##==============================
    # Store simulations parameters #
    ##==============================
    
    # Keep only the simulations after the burn-in period has passed
    if (t > burnin){
      if (is_stephens){
        post_resp_arr[t-burnin, , ] <- post_prob
      }
      if(indiv_alpha) # Keep the adherence parameters for each MCMC draw
        alpha_vec = cbind(alpha_vec,alpha)
      else
        alpha_vec[t-burnin] = alpha[1]
      
      for (m in 1:M)    # Keep counts of L_mn=k for each MCMC draw
        LPost[[m]] <- LPost[[m]] + L_m[[m]]
      CPost <- CPost + C # Keep counts of C_n=k for each MCMC draw
      
      if(is_hard_cluster){
        for (m in 1:M){  
          Llist[[m]][[t-burnin]] = L_m[[m]]
        }
        Clist[[t-burnin]] = C
      }
      
      pi.draws[t-burnin,] <- pi_k     # Keep mixing proportions
      
      G = B = P = PLL = BPR <- 1      # Initialize counts of each source
      for (m in 1:M){
        if (identical(distr[m],"G")){ # Keep mean and variance
          if (D[m] == 1){
            mu.draws[[G]][t-burnin,]   <- melt(Normal$mu[[G]])$value
            tau.draws[[G]][t-burnin,]  <- 1/(melt(Normal$Tau[[G]])$value)
          }else if (D[m] > 1){
            mu.draws[[G]][[t-burnin]]  <- Normal$mu[[G]]
            tau.draws[[G]][[t-burnin]] <- list()
            for (k in 1:K){
              tau.draws[[G]][[t-burnin]][[k]] <- 1/Normal$Tau[[G]][[k]]
            }
          }
          G <- G + 1
        }else if (identical(distr[m],"P")){ # Keep Poisson mean/variance
          if (D[m] == 1){
            l.draws[[P]][t-burnin,]    <- Poisson$l[[P]]
          }else if (D[m] > 1){
            l.draws[[P]][t-burnin, ,]  <- Poisson$l[[P]]
          }
          P <- P + 1
        }else if (identical(distr[m],"B")){ # Keep Binomial success prob
          if (D[m] == 1){
            bp.draws[[B]][t-burnin,]   <- Binom$p[[B]]
          }else if (D[m] > 1){
            bp.draws[[B]][t-burnin, ,] <- Binom$p[[B]]
          }
          B <- B + 1
        }else if (identical(distr[m],"PLL")){
          pll.draws[[PLL]][t-burnin, ,]   <- PoisLL$l[[PLL]]
          PLL <- PLL + 1
        }else if (identical(distr[m],"BPR")){
          w.draws[[BPR]][[t-burnin]] <- BinPR$w[[BPR]]
          BPR <- BPR + 1
        }
      }
    }
    setTxtProgressBar(pb,t)
  }
  close(pb)
  
  cat("Generating final statstics...\n")
  
  ##-----------------------------------------------------
  # Find posterior alpha by averaging over all MCMC     #
  # iterations and computing the 95% credible interval. #
  ##-----------------------------------------------------
  alpha_post <- compute_alpha_posterior(alpha_vec, M, indiv_alpha, mcmc_sims)
  
  # Transpose the matrix when we have not equal adherence parameter
  if (indiv_alpha)
    alpha_vec <- t(alpha_vec)
  
  # Calculate probabilities from counts of source specific clusterings and create
  # cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility.
  Llabels <- list()
  for (m in 1:M){
    LPost[[m]]    <- LPost[[m]] / mcmc_sims
    
    Llabels[[m]]  <- unlist(apply(LPost[[m]], 1, function(x) which(x == max(x, na.rm = TRUE))[1]))
  }
  
  # Calculate probabilities from counts of overall clusterings
  CPost   <- CPost / mcmc_sims
  Clabels <- unlist(apply(CPost, 1, function(x) which(x == max(x, na.rm = TRUE))[1]))
  
  
  # Object to keep input data
  dat         <- NULL
  #dat$X       <- X
  dat$K       <- K
  dat$N       <- N
  dat$M       <- M
  dat$N_sims  <- N_sims
  dat$burnin  <- burnin
  dat$hard_cluster <- is_hard_cluster
  
  # Object to hold all the MCMC draws
  draws       <- NULL
  draws$pi    <- pi.draws
  draws$mu    <- mu.draws
  draws$tau   <- tau.draws
  draws$l     <- l.draws
  draws$bp    <- bp.draws
  draws$pll   <- pll.draws
  draws$w     <- w.draws
  draws$alpha <- alpha_vec
  
  # Object to hold the summaries for the parameters
  summary       <- NULL
  summary$alpha <- alpha_post$alpha_post_mean
  summary$alpha_cred_interv <- alpha_post$alpha_cred_interv
  summary$CPost <- CPost
  summary$LPost <- LPost
  summary$Llab  <- Llabels
  summary$Clab  <- Clabels
  if (is_stephens) # Use Stephens algorithm for relabelling MCMC outputs
    summary$PR  <- post_resp_arr
  
  if (is_hard_cluster){
    clustAssign   <- hard_clustering(Llist, Clist, M, N, mcmc_sims)
    summary$Lbest <- clustAssign$Lbest
    summary$Cbest <- clustAssign$Cbest
  }
  
  # Create a BCC object
  BCC         <- list()
  BCC$dat     <- dat
  BCC$draws   <- draws
  BCC$summary <- summary
  class(BCC) <- "bcc"
  
  return(BCC)
}
