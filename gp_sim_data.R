# Sampling from Gaussian process

# install.packages("rdetools")
library(MASS)
library(rdetools)
library(ggplot2)

N <- 3 # number of patients
T_final <- 60 # total units of time to plot
T_stop <- 10 + sample(c(1,3,5,10,15), size=N,
                      replace=TRUE, prob=c(0.1, 0.2, 0.2, 0.3, 0.2)) # sample treatment stopping times
cat(T_stop)

# Baseline covariate, L ~ N(mu_0,Sigma_0)
mu_0 <- rep(0, 2)
Sigma_0 <- diag(2)
L <- mvrnorm(n = N, mu_0, Sigma_0)

mu_1 <- 10.0 # starting mean
mu_2 <- 0.1 # slope of linear decrease
# mu_3 <- 10.0
mu_4 <- 0.1 # scale of exponential decrease

for(i in 1:N) {

  X <- c(1:T_final) # X is a vector of 60 time points each 1 unit apart

  ind_T_stop <- as.numeric((X >= T_stop[i])) # binary vector indicating time points >= treatment stop time
  
  mu_3 <- mu_1 - mu_2*T_stop[i] # mean at treatment stopping time
  
  # SOFA score, Z ~ N(mu, K)
  # Means for two separate case
  # linear decrease before treatment stop: mu_1 - mu_2*X, and
  # exponential decrease after treatment stop: mu_3*exp(-mu_4*(X - T_stop[i]))
  mu = (1-ind_T_stop) * (mu_1 - mu_2*X) + ind_T_stop * mu_3*exp(-mu_4*(X - T_stop[i]))
  
  K <- rbfkernel(as.matrix(X), sigma=1) # covariance matrix from RBF function
  
  # Generate 10 sample path with mean mu and covariance matrix K
  Z <- mvrnorm(n = 20, mu, K)
  
  # Plotting mean and standard deviation across 20 samples of function
  mean <- colMeans(Z)
  stdev <- apply(Z, 2, sd)
  dat <- data.frame("time" = X,
                    "mean" = mean,
                    "lower" = mean - stdev,
                    "upper" = mean + stdev)
  p <- ggplot(data=dat, aes(x=time, y=mean)) + geom_point() + geom_line() +
    geom_ribbon(aes(ymin=dat$lower, ymax=dat$upper), linetype=2, alpha=0.1) +
    ggtitle(paste("Patient",i)) +
    geom_vline(xintercept=T_stop[i]) +
    geom_text(aes(x=T_stop[i]+1, label="Stop", y=0), angle=90) +
    theme_bw()
  print(p)
}

