##### In this file we do some analysis and study the distribution
library(ConjugateDP)
library(ggpubr)
library(copula)

######### Useful functions

## Marginal distribution for Kn
Kn_marginal <- function(a, b, size){
  Kn <- rep(0, size)
  stir <- log(abs(Stirling1.all(n = size)))
  C_den <- NormConst(a=a, b=b,size=size)
  for(k in 1:size){
    C_num <- NormConst(a=a+k, b=b+1,size=size)
    Kn[k] <- exp(C_num - C_den + stir[k])
  }
  return(Kn)
}

########## Code

################################################################## Shape of the distribution
x <- c(1e-8, 0.0001,0.005, 0.001,0.002, seq(0.01, 90, .005))
K0 <- 40
size <-100
seqs_b <- c(0.03, 0.1, 0.5, 1, 2)
data <- data.frame()
for(b in seqs_b){
  print(b)
  a <- b*K0
  df <- data.frame("alpha" = x,
             "density" = dA(x = x, a = a, b = b, size = size),
             "b" = rep(b, length(x)))
  data <- rbind(data, df)
}



P1 <- ggplot(data = data %>% mutate(b = as.factor(b))) +
  geom_line(aes(x = alpha, y = density, color = b, linetype = b)) +
  theme_bw() +
  scale_color_grey(start = 0, end = 0.4)+
  xlab(expression(alpha))+
  ylab("Density") +
  theme(text=element_text(family="Helvetica"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = .8,
        legend.position="right")+
  facet_grid(~"Prior")


################################################################## Shape of the marginal

## Marginal distribution for Kn
data_marginal <- data.frame()
size = 100
K0 = 40
seqs_b <- c(0.03, 0.1, 0.5, 1, 2)
for(b in seqs_b){
  print(b)
  a <- b*K0
  df <- data.frame("Kn" = seq(1:size),
                   "Probability" = Kn_marginal(a = a, b = b, size = size),
                   "b" = rep(b, size))
  data_marginal <- rbind(data_marginal, df)
}


P2 <- ggplot(data = data_marginal %>% mutate(b = as.factor(b))) +
  geom_point(aes(x = Kn, y = Probability, shape = b, color = b)) +
  theme_bw() +
  xlab(expression(K[n]))+
  ylab("Probability") +
  theme(text=element_text(family="Helvetica"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = .8,
        legend.position="right")+
  facet_grid(~"Prior predictive") +
  scale_shape_manual(values = c(20,2, 3,1, 18))+
  scale_color_grey(start = 0, end = 0.4)



ggarrange(P1, P2)
ggsave("ConjugatePrior_parameters.eps", width = 11, height = 4.4, dpi = 1000)


df <- data.frame("Kn" = seq(1:100),
                 "Probability" = Kn_marginal(a = .2, b = .1, size = 100))
data_marginal <- rbind(data_marginal, df)


### Miscellanea
C1 <- NormConst(a=50, b=1,size=100)
C2 <- NormConst(a=51, b=1,size=100)
C3 <- NormConst(a=52, b=1,size=100)

exp(C2-C1)
sqrt(exp(C3-C1) - (exp(C2-C1))^2)

mean(rA(1e4, a = 50, b =1, size =100))
mean(rA(1e4, a = 50, b =1, size =100)^2)
var(rA(1e4, a = 50, b =1, size =100))


library(ggplot2)


alphas <-rA(n = 1e6, a = 50, b = 1, size = 100)
mean(alphas*(digamma(alphas + 100) - digamma(alphas)))


set.seed(10)
DPsample <- table(simulate_PY(size = 1000, alpha = 5, sigma = 0))
length(DPsample)


alpha_star <- nlminb(
  start = 1,
  function(x) -ConjugateDP:::log_pdf_A(x, a=length(DPsample)+1, b=1, size=1000),
  gradient = function(x) {
    -(a - 1) / x + b * (digamma(x + size) -
                          digamma(x))
  },
  lower = 1e-7
)$par


K0 <- 25
b <- 2
a <- b*K0
size <- 50
Kn <- Kn_marginal(a,b, size)

plot(Kn)
EKn <- sum(Kn*c(1:size))
EKn2 <- sum(Kn*c(1:size)^2)
VarKn <- EKn2 - EKn^2

library(ConjugateDP)
n <- 100000
samples <- rA(1e5, a = 20, size = n)
samples_gamma <- rgamma(1e5, 20, log(n))

plot(density(samples))
lines(density(samples_gamma), col="red")




