#' Plot BCC cluster outputs
#' 
#' Function for plotting BCC cluster assignments. Current implmentation assumes 
#' that we have data sources from gene expression data sources (GE) and DNA 
#' methylation data sources (Meth).
#' 
#' @param bcc_obj Output object of the 'bcc_hts' function.
#' @param original_data Original data sources, the same as the input to the 
#'   bcc_hts function.
#' @param annot Annotations for the plot x, y labels.
#' @param w The fitted coefficients after learning the shape of methylation 
#'   patterns for each genomic region of interest. These can be obtained using 
#'   the 'bpr_optim' function in the BPRMeth package 
#'   (https://github.com/andreaskapou/BPRMeth).
#' @param pca Logical, whether to project the w parameters in 2-D using
#'   Principal Component Analysis or t-SNE method.
#'   
plot_bcc_clusters <- function(bcc_obj, original_data, annot, w = NULL, pca = TRUE){
  # Keep summary data
  D <- bcc_obj$summary
  # Set the x and y labels in the plots
  sub.title <- NULL
  x.lab <- params$x.lab
  y.lab <- params$y.lab
  # Obtain the total number of clusters
  K <- bcc_obj$dat$K
  # For each data source
  for (m in 1:bcc_obj$dat$M){
    X_label <- TRUE
    if (length(D$alpha) == 1){
      adherence <- round(D$alpha[1], 2)
    }else{
      adherence <- round(D$alpha[m], 2)
    }
    # If source is Gaussian assume that we had gene expression data (GE)
    if (identical(original_data[[m]]$distr, "G")){
      X <- original_data[[m]]$x
      main.expr  <- bquote("K562: GE " ~ alpha == .(adherence))
    }else if (identical(original_data[[m]]$distr, "B")){
      # If source is Binomial assume that we had DNA methylation data (Meth)
      X <- original_data$X[[m]]
      X <- X / original_data$r
      main.expr  <- bquote("K562: Meth " ~ alpha == .(adherence))
    }else if (identical(original_data[[m]]$distr, "BPR")){
      # If source is Binomial assume that we had DNA methylation data (Meth)
      if (!is.null(w)){
        X <- w
        main.expr  <- bquote("K562: Meth " ~ alpha == .(adherence))
      }else{
        X_label <- FALSE
      }
    }
    
    if (X_label){
      if (NCOL(X) > 2){
        if (pca){
          X <- prcomp(X)$x
        }else{
          require(tsne)
          max_iter <- 300
          perplexity <- 100
          X <- tsne(X, max_iter = max_iter, perplexity = perplexity)
        }
      }
  
      ## Start Plot
      
      # Colours denote the overall clusterings and symbols denote 
      # the source-specific clusterings.
      cols <- c("blue2", "forestgreen", "red2", "goldenrod3")
      cex_t <- c(0.2, 0.5, 0.5, 0.5, 0.2, 0.5, 0.5, 0.5, 0.2)
      cex_t_k <- c(0.5, 0.5, 0.5, 0.2)
      
      i = j <- 1
      par(cex=1.15, mai=c(1.37,1.37,.7,.3))
      if (NCOL(X) > 1){
        plot(X[,1], X[,2], cex=0, xlab=NA, ylab=NA, main=NULL)
        if (pca){
          mtext(side = 1, "Methylation PC1", line = 2.3, cex = 1.5)
          mtext(side = 2, "Methylation PC2", line = 2.3, cex = 1.5)
        }else{
          mtext(side = 1, "t-SNE dim 1", line = 2.3, cex = 1.5)
          mtext(side = 2, "t-SNE dim 2", line = 2.3, cex = 1.5)
        }
        title(main=main.expr, sub = sub.title, line = 1, cex.main=1.3)
      }else{
        plot(X, seq(1:length(X)), cex=0, xlab=NA, ylab=NA, yaxt='n', main=NULL)
        mtext(side = 1, x.lab, line = 2.3, cex = 1.5)
        mtext(side = 2, y.lab, line = 1.5, cex = 1.5)
        title(main=main.expr, sub = sub.title, line = 1, cex.main=1.3)
      }
      if (bcc_obj$dat$hardCl){
        for (k in 1:K){
          if (NCOL(X) > 1){
            points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,1]==1,1], 
                   X[D$Cbest[,k]==1 & D$Lbest[[m]][,1]==1,2], 
                   col = cols[k], pch = 2, cex = cex_t[i])
            i <- i + 1
            points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,2]==1,1], 
                   X[D$Cbest[,k]==1 & D$Lbest[[m]][,2]==1,2], 
                   col = cols[k], pch = 8, cex = cex_t[i])
            i <- i + 1
            points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,3]==1,1], 
                   X[D$Cbest[,k]==1 & D$Lbest[[m]][,3]==1,2], 
                   col = cols[k], pch = 1, cex = cex_t[i])
            i <- i + 1
            if (bcc_obj$dat$K == 4){
              points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,4]==1,1], 
                     X[D$Cbest[,k]==1 & D$Lbest[[m]][,4]==1,2], 
                     col = cols[k], pch = 22, cex = cex_t_k[j])
              j <- j + 1
            }
          }else{
            points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,1]==1], 
                   which(D$Cbest[,k]==1 & D$Lbest[[m]][,1]==1), 
                   col = cols[k], pch = 2, cex = cex_t[i])
            i <- i + 1
            points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,2]==1], 
                   which(D$Cbest[,k]==1 & D$Lbest[[m]][,2]==1), 
                   col = cols[k], pch = 8, cex = cex_t[i])
            i <- i + 1
            points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,3]==1], 
                   which(D$Cbest[,k]==1 & D$Lbest[[m]][,3]==1), 
                   col = cols[k], pch = 1, cex = cex_t[i])
            i <- i + 1
            if (bcc_obj$dat$K == 4){
              points(X[D$Cbest[,k]==1 & D$Lbest[[m]][,4]==1], 
                     which(D$Cbest[,k]==1 & D$Lbest[[m]][,4]==1), 
                     col = cols[k], pch = 22, cex = cex_t_k[j])
              j <- j + 1
            }
          }
        }
      }else{
        for (k in 1:K){
          if (NCOL(X) > 1){
            points(X[D$Clab==k & D$Llab[[m]]==1,1], 
                   X[D$Clab==k & D$Llab[[m]]==1,2], 
                   col = cols[k], pch = 2, cex = cex_t[i])
            i <- i + 1
            points(X[D$Clab==k & D$Llab[[m]]==2,1], 
                   X[D$Clab==k & D$Llab[[m]]==2,2], 
                   col = cols[k], pch = 8, cex = cex_t[i])
            i <- i + 1
            points(X[D$Clab==k & D$Llab[[m]]==3,1], 
                   X[D$Clab==k & D$Llab[[m]]==3,2], 
                   col = cols[k], pch = 1, cex = cex_t[i])
            i <- i + 1
            if (K == 4){
              points(X[D$Clab==k & D$Llab[[m]]==4,1], 
                     X[D$Clab==k & D$Llab[[m]]==4,2], 
                     col = cols[k], pch = 22, cex = cex_t_k[j])
              j <- j + 1
            }
          }else{
            points(X[D$Clab==k & D$Llab[[m]]==1], 
                   which(D$Clab==k & D$Llab[[m]]==1), 
                   col = cols[k], pch = 2, cex = cex_t[i])
            i <- i + 1
            if (k == 2){
              XX <- X[D$Clab==k & D$Llab[[m]]==2]
              points(XX - 0.12, 
                     which(D$Clab==k & D$Llab[[m]]==2), 
                     col = cols[k], pch = 8, cex = cex_t[i])
            }else{
              points(X[D$Clab==k & D$Llab[[m]]==2], 
                     which(D$Clab==k & D$Llab[[m]]==2), 
                     col = cols[k], pch = 8, cex = cex_t[i])
            }
            i <- i + 1
            points(X[D$Clab==k & D$Llab[[m]]==3], 
                   which(D$Clab==k & D$Llab[[m]]==3), 
                   col = cols[k], pch = 1, cex = cex_t[i])
            i <- i + 1
            if (K == 4){
              points(X[D$Clab==k & D$Llab[[m]]==4], 
                     which(D$Clab==k & D$Llab[[m]]==4), 
                     col = cols[k], pch = 22, cex = cex_t_k[j])
              j <- j + 1
            }
          }
        } # End for loop 1:K
      } # End if hard_cluster
    } # End if X_label
  }# End for loop m = 1:M
}