#' Align cluster indices
#' 
#' \code{align_clusters} function helps to align cluster indices after each MCMC
#' simulation for each data source.
#' 
#' @param Z1 Previous cluster assignments.
#' @param Z2 Current cluster assignments
#' @param type Object type of the cluster assignemnts, either 'vector' or
#'   'matrix'.
#'   
#' @return The aligned indices of the current cluster assignments.
#' 
align_clusters <- function(Z1, Z2, type = "vec"){
  if(type == "vec"){
    for(k in 1:length(unique(Z1))){ # For each cluster k in previous Cluster
      # Find Max 
      Max <- sum(Z1==k & Z2==k)/(.01 + sum(Z1==k) + sum(Z2==k))
      for(tempk in  1:length(unique(Z2))){ # For each cluster k in current Cluster
        # Check if the proportions are higher than Max
        if( (sum(Z1==k & Z2==tempk)/(.01 + sum(Z1==k) + sum(Z2==tempk))) > Max){
          # Get the proportion that the two cluster indices are the same
          Max <- sum(Z1==k & Z2==tempk)/(.01 + sum(Z1==k) + sum(Z2==tempk))
          dummy <- (Z2==k)      # Keep indices that do not match
          Z2[Z2==tempk] <- k    # Swap the incorrect indices
          Z2[dummy] <- tempk    # Swap the incorrect indices
        }
      }
    }
  }else if(type == "mat"){
    for(k in 1:dim(Z1)[2]){         # For each cluster k in previous Cluster
      for(tempk in  1:dim(Z2)[2]){  # For each cluster k in current Cluster
        Max <- sum(Z1==Z2)          # Number of matches between the cluster indices
        Z2dummy <- Z2               # Keep the current indices in a dummy variable
        Z2dummy[,k] = Z2[,tempk]    # Swap the incorrect indices
        Z2dummy[,tempk] = Z2[,k]    # Swap the incorrect indices
        if(sum(Z1==Z2dummy) > Max){ # If the swaps make a better alignment, update indices
          Z2 <- Z2dummy
        }
      }
    }
  }
  return(Z2) # Return the aligned cluster indices
}