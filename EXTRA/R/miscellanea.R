
out <- rStirgamma(nsamples = 1e4, a = a, b = b, m = m)
plot(density(rStirgamma(nsamples = 1e4, a = a, b = b, m = m)))

#f <-function(x) {- log_pdf_logSg(x, a, b, m) - log(x^2)}
#X <- seq(-5, 5,length.out = 1000)
#plot(X, sapply(X, function(x) -f(x)), type = "l")

#curve(f(x), from  = -5, to = 5)
#
# a <- 2; b <- 1; m <- 100
# bounds <- getBounds_uv(a, b, m)
# nsamples <- 1e4
# out <- c(rSg_logSg(nsamples, a, b, m, Mu= bounds$Mu, Mv_minus = bounds$Mv_minus, Mv_plus = bounds$Mv_plus))
#
# #density(exp(out)/log(m))
#
# #quantile(out)
# x <-log(rStirgamma(nsamples, a, b, m) * log(m))
# plot(sort(out), sort(x))
#
#
# plot(density(rStirgamma(1e4, a, b,  m)), col = "red")
# lines(density(rgamma(1e4, a - b, b * log(m))))


#
# x_star <- stats::nlminb(
#   start = 1,
#   function(x) - log_pdf_logSg(x, a, b, m),
#   gradient = function(x) {
#     - a + b * exp(x) * (digamma(exp(x) / log(m) + m) -
#                           digamma(exp(x) / log(m))) / log(m)
#   },
#   lower = 0,
#   control = list(abs.tol = 1e-12)
# )$par
#
# log_pdf_logSg(x_star, a, b, m)
#
#
# x2logSg_plus <- -stats::nlminb(
#   start = 1,
#   lower = 0,
#   control = list(abs.tol = 1e-13),
#     function(x) - log_pdf_logSg(x, a, b, m) - log(x^2),
#   gradient = function(x) {
#     - a + b * exp(x) * (digamma(exp(x) / log(m) + m) -
#                           digamma(exp(x) / log(m))) / log(m) - 2 / x
#   },
# )$objective
#
#
# x2logSg_minus <- -stats::nlminb(
#   start = -1,
#   upper = 0,
#   control = list(abs.tol = 1e-13),
#   function(x) - log_pdf_logSg(x, a, b, m) - log(x^2),
#   gradient = function(x) {
#     - a + b * exp(x) * (digamma(exp(x) / log(m) + m) -
#                           digamma(exp(x) / log(m))) / log(m) - 2 / x
#   },
# )$objective
#
# Mv_plus <- exp((x2logSg_plus - log_pdf_logSg(x_star, a, b, m)) / 2)
# Mv_minus <- - exp((x2logSg_minus - log_pdf_logSg(x_star, a, b, m)) / 2)
#
#


logf <- function(x, a,b,m){
  (a-1) * log(x) - b * (lgamma(x + m) - lgamma(x))
}

m <- 100
b <- 1
a <- 2
curve(from = 0, to = .1, expr = logf(x, a,b, m))

x <- rSg(1e5, 2.5, 1, 10)
y <- rSg(1e5, 3, 1, 10)
plot(density(x/(x+y)))


P <- function(x, a, b, m){
  (a-1) * log(x) - b* (lgamma(x + m) - lgamma(x))
}

h <- function(alpha, a, b, m){
  (a-1)/alpha - b * (digamma(alpha + m) - digamma(alpha))
}

hprime <- function(alpha, a, b, m){
  -(a-1)/(alpha^2) - b * (trigamma(alpha + m) - trigamma(alpha))
}


fsecond <- function(x, a, b, m){
  h(x, a,b, m)^2 + hprime(x,a,b, m)
}


m <- 100
b <- 1
a <- 2
curve(from = 1e-8, to = 5, expr = log(0.5) + 2 * log(x) + P(x,a,b,m) + log(fsecond(x, a,b, m)))

log(fsecond(0.001, a,b, m))

lbeta(a-b, m*b - a + b)

ConjugateDP:::ComputeStirgammaNormConst(a,b,m)


m <- 1000
b <- 0.7
a <- 1.2
x <- rbeta(1e6, a - b, (m-1)*b-a+b)
y <- rSg(1e6, a,b,m)
plot(density(y))
lines(density(exp(lgamma(m)/(m-1))* x/(1-x)), col = "red")


r <- 2
all(sapply(seq(0, 100000, length.out = 1e5), function(x) sum(log(x + 1:(m-1))) > (m-1) * log(x + r)))


x <- 2e6
sum(log(x + 1:(m-1))) > (m-1) * log(x + r)



S <- ConjugateDP:::ComputeStirgammaNormConst(a,b,m)
B <- (a-m*b)/(m-1) * lgamma(m) + lbeta(a -b, m * b - a)

1/exp(B - S)


f <- function(x, m) {
  xx <- sapply(x, function(i) sum(log(i + 1:(m-1))))
  exp(xx/(m-1)) - x
}
curve(from = 0, to = 10000, f(x, m), n = 1000)


nlminb(start =1,
       objective =function(x) f(x, m),lower = 0)



alpha <- 0.02
(m-1) * log(alpha + 3)
sum(log(alpha + (1:m-1)))


########## Let's try accept/reject from the BetaPrime
logA <- function(x, b, m) {
  -b * (lgamma(x + m) - lgamma(x + 1)) + ((m * b - b - 1)/(m-1)) * lgamma(m) +
    (m * b - 1) * log(1 + x / exp(lgamma(m)/(m-1)))
}


sampleSg <- function(nsamples, a, b, m){
  out <- rep(NA, nsamples)
  for(i in 1:nsamples) {
    accept <- 0
    while(accept == 0) {
      # Sample the generalized beta prime
      x <- rbeta(1, a - b, m * b - a)
      y <- exp(lgamma(m)/(m-1)) * x /(1 - x)
      # Sample from the uniform
      u <- runif(1)
      # Evaluate the acceptance ratio
      if(log(u) <= logA(y, b, m)){
        out[i] <- y
        accept <- 1
      }
    }
  }
  return(out)
}

a <- 0.5
b <- 0.2
m <- 10000
AccrejSg <- sampleSg(1000, a = a, b = b, m=m)
runifSg <- rStirgamma(20000, a = a, b = b, m = m)

plot(density(runifSg))
lines(density(AccrejSg), col = "red")


getBounds_uv <- function(a, b, m) {
  ###############################################################
  # Case 1 - a - b > 1.
  ###############################################################
  # Calculate the maximum value of the log density
  x_star <- stats::nlminb(
    start = 1,
    function(x) - log_pdf_logSg(x, a, b, m),
    gradient = function(x) {
      - a + b * exp(x) * (digamma(exp(x) / log(m) + m) -
                            digamma(exp(x) / log(m))) / log(m)
    },
    control = list(abs.tol = 1e-12)
  )$par
  Mu <- log_pdf_logSg(x_star, a, b, m)

  # Calculate the maximum value of the log density + log(x^2) when x > 0
  Mv_plus <- - stats::nlminb(
    start = 1,
    lower = 0,
    control = list(abs.tol = 1e-13),
    function(x) - log_pdf_logSg(x, a, b, m) - log(x^2),
    gradient = function(x) {
      - a + b * exp(x) * (digamma(exp(x) / log(m) + m) -
                            digamma(exp(x) / log(m))) / log(m) - 2 / x
    },
  )$objective

  #Calculate the maximum value of the log density + log(x^2) when x < 0
  Mv_minus <- -stats::nlminb(
    start = -1,
    upper = 0,
    control = list(abs.tol = 1e-13),
    function(x) - log_pdf_logSg(x, a, b, m) - log(x^2),
    gradient = function(x) {
      - a + b * exp(x) * (digamma(exp(x) / log(m) + m) -
                            digamma(exp(x) / log(m))) / log(m) - 2 / x
    },
  )$objective

  return(list(Mu = Mu, Mv_plus = Mv_plus, Mv_minus = Mv_minus))
}


# //
#   // // //[[Rcpp::export]]
# // // double log_pdf_logSg(double x, double a, double b, double m){
#   // //   double logf;
#   // //   logf = a * x - b * (lgamma(exp(x) / log(m) + m) - lgamma(exp(x) / log(m)));
#   // //   return logf;
#   // // }
# //
#   // //[[Rcpp::export]]
# // double log_pdf_logSg(double x, double a, double b, double m){
#   //   double logf;
#   //   logf = a * x - b * (lgamma(exp(x) + m) - lgamma(exp(x)));
#   //   return logf;
#   // }


logA <- function(x, b, m) {
  #  -b * (lgamma(x + m) - lgamma(x + 1)) + ((m * b - b - 1)/(m-1)) * lgamma(m) +
  #    (m * b - 1) * log(1 + x / exp(lgamma(m)/(m-1)))
  -b * (lgamma(x + m) - lgamma(x + 1)) + b * (m - 1) * log(x + exp(lgamma(m)/(m-1)))
}


sampleSg <- function(nsamples, a, b, m){
  out <- rep(NA, nsamples)
  for(i in 1:nsamples) {
    accept <- 0
    while(accept == 0) {
      # Sample the generalized beta prime
      x <- rbeta(1, a - b, m * b - a)
      y <- exp(lgamma(m)/(m-1)) * x /(1 - x)
      # Sample from the uniform
      u <- runif(1)
      # Evaluate the acceptance ratio
      if(log(u) <= logA(y, b, m)){
        out[i] <- y
        accept <- 1
      }
    }
  }
  return(out)
}


trySg <- function(ntrials, a, b, m){
  accept <- rep(0, ntrials)
  for(i in 1:ntrials) {
    # Sample the generalized beta prime
    x <- rbeta(1, a - b, m * b - a)
    y <- exp(lgamma(m)/(m-1)) * x /(1 - x)
    # Sample from the uniform
    u <- runif(1)
    # Evaluate the acceptance ratio
    if(log(u) <= logA(y, b, m)){
      #out[i] <- y
      accept[i] <- 1
    }
  }
  return(mean(accept))
}



# a <- 0.6
# b <- 0.1
# m <- 1000
# nsamples <- 100000
# sg <- rStirgamma(nsamples, a, b, m)
# y <- sampleSg(nsamples, a, b, m)
# plot(density(y))
# lines(density(sg), col ="red")
#
# x <- rbeta(nsamples, a-b, m * b - a)
# bp <- x/(1-x)
# r <- exp(lgamma(m)/(m-1))
# plot(density(sg))
# lines(density(r * bp), col = "red")
#
#
#
# S <- ConjugateDP:::ComputeStirgammaNormConst(a,b,m)
# B <- (a-m*b)/(m-1) * lgamma(m) + lbeta(a -b, m * b - a)
# 1/exp(B- S)
#
y1 <- c(rSg_beta_prime(10000,a,b,m))
y2 <- sampleSg(10000, a, b, m)
y3 <- rStirgamma(10000, a, b, m)
plot(density(y1))
lines(density(y2), col = "red")
lines(density(y3), col = "green")

# trySg(10000, a, b, m)
