#---------------------------------
# Set the working directory and load any required libraries 
#---------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}

require(coda)
require(R.utils)

# Source the 'lib' directory
R.utils::sourceDirectory("lib", modifiedOnly=FALSE) 

# Set seed for reproducible results
set.seed(12345)

##-----------------------------
# Initialize main variables   #
##-----------------------------
M      <- 3      # Number of sources
K      <- 2      # Number of clusters
N      <- 300    # Number of objects
N_sims <- 5000   # Set the number of simulations
burnin <- 2000   # Set how many samples should be burned in

# Model hyperparameters
hyper <- list(aBeta = 1, bBeta = 1, aDir = 1)


##------------------------------------
#       Generate synthetic data      #
##------------------------------------

# Vector with the source distribution list, i.e. we create data from 3 
# different sources, 1st source: Gaussian, 2nd source: Poisson and 
# 3rd source: Binomial
distr <- c("G", "P", "B")

# Generate the overall clusterings
C_true <- matrix(0, nrow = N, ncol = K) 
C_true[1:N/2, 1] <- 1        # Set half of them to overall cluster 1
C_true[((N/2)+1):N, 2] <- 1  # The other half to overall cluster 2

# Generate alpha ~ U[0.5,1]
a_true <- runif(1, min=.5, max=1)

# Give values to the parameters of the distributions to be generated
mu_init = std_init = l_init = p_init <- list()
mu_init[[1]] <- c(1.5, -1.5)
std_init[[1]] <- c(1, 1)
l_init[[1]] <- c(7, 22)
p_init[[1]] <- c(.2, .65)
# Total number of trials for source Binomial
r1 <- rbinom(n=N, size=50, prob=.8)  
r <- matrix(r1, ncol=1)

# Generate the data
synth_data <- gen_bcc_data(N = N, M = M, K = K, C_true = C_true, alpha = a_true, 
                           distr = distr, mu_init = mu_init, std_init = std_init, 
                           l_init = l_init, p_init = p_init, basis = basis, r = r)

# Obtain the true source-specific clusterings
L_true <- synth_data$L_true

##---------------------------------------------
# Initialize parameters for each data source  #
##---------------------------------------------
params      <- list()
params[[1]] <- list()
params[[2]] <- list()
params[[3]] <- list(r = r)


##----------------------------
#     Run BCC model          #
##----------------------------
message("Start Updated BCC")
bcc_out <- bcc_hts(X           = synth_data$X,
                   K           = K,
                   indiv_alpha = FALSE,
                   hyper       = hyper,
                   params      = params,
                   N_sims      = N_sims,
                   burnin      = burnin,
                   is_parallel = TRUE,
                   no_cores    = 3)


# Compute source-specific and overall clustering assignment errors
source_error <- compute_source_error(L_true = L_true, L_post = bcc_out$summary$LPost)
overall_error <- compute_overall_error(C_true = C_true, C_post = bcc_out$summary$CPost)


##---------------------------------------------
# Some plots to assess MCMC simulation and    #
# convergence to posterior distribution.      #
##---------------------------------------------

#-- Trace plot of adherence parameter alpha
alpha_draws <- as.data.frame(bcc_out$draws$alpha)
names(alpha_draws) <- "alpha"
plot(coda::mcmc(alpha_draws))

#-- Trace plots of mixing proportions \pi
pi_draws <-  as.data.frame(bcc_out$draws$pi)
names(pi_draws) <- c("pi1", "pi2")
plot(coda::mcmc(pi_draws))

#-- Trace plots of means of Gaussian source
mu_draws <- as.data.frame(bcc_out$draws$mu[[1]])
names(mu_draws) <- c("mu 1", "mu 2")
plot(mcmc(mu_draws))
