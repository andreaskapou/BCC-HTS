#' Hard clustering by least squares
#' 
#' \code{hard_clustering} functions aggregates over the MCMC iterations and 
#' produces a single (hard) clustering, proposed by Dahl (2006).
#' 
#' @param Llist List of source-specific clustering assignments from MCMC 
#'   simulations.
#' @param Clist List of overall clustering assignments from MCMC simulations.
#' @param M Number of sources.
#' @param N Number of objects.
#' @param mcmc_sims Number of MCMC simulations.
#' 
#' @return Hard source-specific and overall clustering assignments.
#' 
hard_clustering <- function(Llist, Clist, M, N, mcmc_sims){
  # List of best clustering for each data source
  Lbest <- list()
  
  # Choose best source clustering
  for(m in 1:M){
    Lkern <- Llist[[m]][[floor(mcmc_sims / 5)]] %*% t(Llist[[m]][[floor(mcmc_sims / 5)]])
    for(t in floor(mcmc_sims / 5 + 1):mcmc_sims){
      Lkern <- Lkern + Llist[[m]][[t]] %*% t(Llist[[m]][[t]])
    }
    Lkern <- Lkern / (mcmc_sims - floor(mcmc_sims / 5) + 1)
    CountLbest <- N ^ 2 + 1
    for(t in floor(mcmc_sims / 5):mcmc_sims){
      CountL <- norm(Lkern - Llist[[m]][[t]] %*% t(Llist[[m]][[t]]), 'F') ^ 2
      if(CountL < CountLbest){
        Lbest[[m]] <- Llist[[m]][[t]]
        CountLbest <- CountL
      }
    }
  }
  
  # Choose best overall clustering
  Ckern <- Clist[[floor(mcmc_sims / 5)]] %*% t(Clist[[floor(mcmc_sims / 5)]])
  for(t in floor(mcmc_sims / 5 + 1):mcmc_sims){
    Ckern <- Ckern + Clist[[t]] %*% t(Clist[[t]])
  }
  Ckern <- Ckern / (mcmc_sims - floor(mcmc_sims / 5) + 1)
  CountCbest <- N ^ 2 + 1
  for(t in floor(mcmc_sims / 5):mcmc_sims){
    CountC <- norm(Ckern - Clist[[t]] %*% t(Clist[[t]]), 'F')^2
    if(CountC < CountCbest){
      Cbest <- Clist[[t]]
      CountCbest <- CountC
    }
  }
  return(list(Lbest = Lbest, Cbest = Cbest))
}