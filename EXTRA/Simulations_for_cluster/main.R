##### Main script to run the cluster

## Get the slurm id
#slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
slurm_id <- 4

## Import the libraries
library(foreach)
library(doParallel)
library(nleqslv)
library(coda)

## Import the functions
source("../R/DPMM_v2.R")
source("Simulations_functions.R")
sourceCpp("../src/useful_functions.cpp")

## Size of the simulated datasets
n_list <- c(100, 500, 1000, 10000)

## Create the data
set.seed(42)
n <- n_list[slurm_id]
data <- sample_data(size = n, sigma = 1, K = 2)

## hyperparameters of the function
# Try three different priors
y <- scale(data[,-3])

## Prior of the base measure
p <- ncol(y)
nu0 <- kappa0 <- p
S0 <- diag(p)
mu0 <- rep(0,p)

## Prior on the alpha parameter
a <- 2
b <- 1

R <- 10000
burnin <- 1000
prior_moments <- moments_A(a=a, b=b, size = n)
a_gamma <-  prior_moments[1]^2/prior_moments[2]
b_gamma <-  prior_moments[1]/prior_moments[2]


## Cores
n_cores <- 6
registerDoParallel(n_cores)


## Run the simulation in parallel now
set.seed(42, kind = "L'Ecuyer-CMRG")

out <- foreach(j = 1:6) %dopar% {
  if(j == 1) {
    ### My conjugate prior with a = 2, b = 1
    m <- DPMM(y = y, R = R, burnin = burnin, a = a, b = b, nu0 = nu0, mu0 = mu0,
         S0 = S0, kappa0 = kappa0)
  } else if (j == 2) {
    ### My conjugate prior with a = 1, b = .1
    DPMM(y = y, R = R, burnin = burnin, a = .2, b = .1, nu0 = nu0, mu0 = mu0,
         S0 = S0, kappa0 = kappa0)
  } else if (j == 3){
    ### My conjugate prior with a = 10, b = 1
    DPMM(y = y, R = R, burnin = burnin, a = 10, b = 1, nu0 = nu0, mu0 = mu0,
         S0 = S0, kappa0 = kappa0)
  } else if (j == 4) {
    # Escobar and West
    DPMM(y = y, R = R, burnin = burnin, a = 2, b = 4, nu0 = nu0, mu0 = mu0,
         S0 = S0, kappa0 = kappa0, Gamma_prior = T)
  } else if (j == 5) {
    alpha_K <- nleqslv(1, fn, n = n, K = a/b)$x
    DPMM(y = y, R = R, burnin = burnin, alpha = alpha_K, alpha_fixed = T,
                       b=b_gamma, nu0 = nu0, mu0 = mu0, S0 = S0, kappa0 = kappa0)
  } else if (j == 6) {
    DPMM(y = y, R = R, burnin = burnin, alpha = 1, alpha_fixed = T,
                       b=b_gamma, nu0 = nu0, mu0 = mu0, S0 = S0, kappa0 = kappa0)
  }

}

names(out) <- c("DY_2_1", "DY_.2_.1", "DY_10_1", "EW", "Fixed_K", "Fixed_1")

## Save the output now
out_name <- paste0("output/comp2/simulation_result_", slurm_id, ".rds.gzip")
saveRDS(out,file = out_name, compress = "gzip")
plot(table(out$DY_2_1$Kn))
plot(table(out$DY_.2_.1$Kn))
plot(table(out$DY_10_1$Kn))
plot(table(out$EW$Kn))
plot(table(out$Fixed_1$Kn))


