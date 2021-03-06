########Q and A learning doubly robust #########
################################################
rm(list = ls())
source('simulation.stop.R')

#Based on the method proposed by Zhang et al 2013 
#(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3843953/)
library(DynTxRegime) #tutorial here - shiny::runGitHub('DynTxRegimeTutorial','ShupingR')
library(rgenoud)

n = 1000
data <- as.data.frame(generate(n = n))
head(data)

# Define the propensity for treatment model and methods.
moPropen1 <- buildModelObj(model = ~ L1,
                          solver.method = 'glm',
                          solver.args = list('family'='binomial'),
                          predict.method = 'predict.glm',
                          predict.args = list(type='response'))
moPropen2 <- buildModelObj(model = ~ L2,
                           solver.method = 'glm',
                           solver.args = list('family'='binomial'),
                           predict.method = 'predict.glm',
                           predict.args = list(type='response'))
moPropen <- list(moPropen1, moPropen2)

# outcome model second stage
### specify the covariates of the main effects component of the outcome regression model
moMain2 <- buildModelObj(model = ~ L1 + A1 + A1:L1 + L2 + A2 + A2:L2 + A1:L2,
                         solver.method = 'lm')

### specify the covariates of the contrasts component of the outcome regression model
moCont2 <- buildModelObj(model = ~ A2,
                         solver.method = 'lm')

# outcome model first stage
moMain1 <- buildModelObj(model = ~ L1 + A1 + A1:L1,
                         solver.method = 'lm')
moCont1 <- buildModelObj(model = ~ A1,
                         solver.method = 'lm')

moMain <- list(moMain1, moMain2)
moCont <- list(moCont1, moCont2)

# regime function second stage
regime2 <- function(eta2, data) {
  tst <- {as.numeric(I(data$L2 > eta2 & data$A1 == 0))} #those with survival score at time point 2 is lower than the 2nd threshold given 1
  rec <- rep(1, nrow(x = data))
  rec[!tst] <- 0 
  return( rec )
}
# regime function first stage
regime1 <- function(eta1, data) {
  tst <- {as.numeric(I(data$L1 > eta1))}
  rec <- rep(1, nrow(x = data))
  rec[!tst] <- 0
  return( rec )
}
regimes <- list(regime1, regime2)

#### Analysis using AIPW
fit_AIPW <- optimalSeq(moPropen = moPropen,
                       moMain = moMain, moCont = NULL,
                       regimes = regimes,
                       data = data, response = data$Y, txName = c('A1', 'A2'),
                       Domains = cbind(rep(200, 2), rep(400, 2)),
                       pop.size = 1000, 
                       starting.values = rep(300, 2))

##Available methods
# Coefficients of the regression objects
coef(object = fit_AIPW)

# Description of method used to obtain object
DTRstep(object = fit_AIPW)
# Estimated value of the optimal treatment regime for training set
estimator(x = fit_AIPW)
# Value object returned by regression methods
fitObject(object = fit_AIPW)
# Retrieve the results of genetic algorithm
genetic(object = fit_AIPW)
# Estimated optimal treatment and decision functions for training data
optTx(x = fit_AIPW)
# Estimated optimal treatment and decision functions for new data
optTx(x = fit_AIPW, newdata = data)
# Value object returned by outcome regression method
outcome(object = fit_AIPW)
# Plots if defined by regression methods
dev.new()
par(mfrow = c(2,4))
plot(x = fit_AIPW)
plot(x = fit_AIPW, suppress = TRUE)
# Retrieve the value object returned by propensity regression method
propen(object = fit_AIPW)
# Show main results of method
show(object = fit_AIPW)
# Show summary results of method
summary(object = fit_AIPW)
#estimated values for eta 
regimeCoef(object = fit_AIPW)
