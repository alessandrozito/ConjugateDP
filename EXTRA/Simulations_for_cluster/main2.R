##### Main script to run the cluster

## Get the slurm id
slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#slurm_id <- 1

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
n_list <- rep(c(100, 500, 1000, 10000), 2)
K_list <- sort(rep(c(2,5),4))

# Prior hyperparameters for the distribution based on slurm id
a1_list <- c(2,2,2,2,5,5,5,5)

## Create the data
n <- n_list[slurm_id]
K <- K_list[slurm_id]


## Prior on the alpha parameter
a <- a1_list[slurm_id]
b <- 1

R <- 10000
burnin <- 1000

## Cores
n_cores <- 20
n_sims <- 20
registerDoParallel(n_cores)


## Run the simulation in parallel now
set.seed(42, kind = "L'Ecuyer-CMRG")
out <- foreach(j = 1:n_sims) %dopar% {
  ## Create the dataset
  data <- sample_data(size = n, sigma = 1, K = K)
  # Try three different priors
  y <- scale(data[,-3])

  ## Prior of the base measure
  p <- ncol(y)
  nu0 <- kappa0 <- p
  S0 <- diag(p)
  mu0 <- rep(0,p)

  ### My conjugate prior with a = 2, b = 1
    out1 <- DPMM(y = y, R = R, burnin = burnin, a = a, b = b, nu0 = nu0, mu0 = mu0,
              S0 = S0, kappa0 = kappa0)
  ### My conjugate prior with a = 1, b = .1
    out2 <- DPMM(y = y, R = R, burnin = burnin, a = 1, b = .1, nu0 = nu0, mu0 = mu0,
         S0 = S0, kappa0 = kappa0)
  ### My conjugate prior with a = 10, b = 1
    out3 <- DPMM(y = y, R = R, burnin = burnin, a = 10, b = 1, nu0 = nu0, mu0 = mu0,
         S0 = S0, kappa0 = kappa0)
  ### Escobar and West
    out4 <- DPMM(y = y, R = R, burnin = burnin, a = 2, b = 4, nu0 = nu0, mu0 = mu0,
         S0 = S0, kappa0 = kappa0, Gamma_prior = T)
  ### Fixed alpha = 1
    out5 <- DPMM(y = y, R = R, burnin = burnin, alpha = 1, alpha_fixed = T,
         b=b_gamma, nu0 = nu0, mu0 = mu0, S0 = S0, kappa0 = kappa0)

    list("data" = data, "DY1"= out1, "DY2" = out2, "DY3" = out3, "EW" = out4, "fixed" = out5)
}

#names(out) <- c("DY_2_1", "DY_.2_.1", "DY_10_1", "EW", "Fixed_K", "Fixed_1")

## Save the output now
out_name <- paste0("output/comp3/simulation_result_", slurm_id, ".rds.gzip")
saveRDS(out,file = out_name, compress = "gzip")


