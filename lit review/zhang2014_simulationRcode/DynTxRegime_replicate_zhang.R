########Q and A learning doubly robust #########
################################################
rm(list = ls())
source('simulation_replicate_zhang.R')

#Based on the method proposed by Zhang et al 2013 
#(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3843953/)
library(DynTxRegime) #tutorial here - shiny::runGitHub('DynTxRegimeTutorial','ShupingR')
library(rgenoud)

n = 500
data <- as.data.frame(generate(n = n))
head(data)

# Define subsets of patients to limit available treatments 
fSet1 <- function(data){
  subsets <- list(list("s1",c(0L,1L)))
  
  txOpts <- rep(x = 's1', times = nrow(x = data))
  
  return(list("subsets" = subsets, "txOpts" = txOpts))
}

fSet2 <- function(data){
  subsets <- list(list("s1",1L),
                  list("s2",c(0L,1L)))
  
  txOpts <- rep(x = 's2', times = nrow(x = data))
  txOpts[data$A1 == 1L] <- "s1"
  
  return(list("subsets" = subsets, "txOpts" = txOpts))
}

# Define the propensity for treatment model and methods.
p1 <- modelObj::buildModelObj(model = ~ L1,
                              solver.method = 'glm',
                              solver.args = list(family='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))
p2 <- modelObj::buildModelObj(model = ~ L2,
                              solver.method = 'glm',
                              solver.args = list(family='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))


# outcome model second stage
### specify the covariates of the main effects component of the outcome regression model
q2Main <- modelObj::buildModelObj(model = ~ L1 + L2,
                                  solver.method = 'lm',
                                  predict.method = 'predict.lm')
### specify the covariates of the contrasts component of the outcome regression model
q2Cont <- modelObj::buildModelObj(model = ~ L2,
                                  solver.method = 'lm',
                                  predict.method = 'predict.lm')

# outcome model first stage
q1Main <- modelObj::buildModelObj(model = ~ L1,
                                  solver.method = 'lm',
                                  predict.method = 'predict.lm')
q1Cont <- modelObj::buildModelObj(model = ~ L1,
                                  solver.method = 'lm',
                                  predict.method = 'predict.lm')

# regime function second stage
regime2 <- function(eta2, data){return(data$A1 + {1L-data$A1}*{data$L2 < eta2}) }
# regime function first stage
regime1 <- function(eta1, data){return(as.integer(x = {data$L1 < eta1})) }

#### Analysis using AIPW
fit_AIPW <- optimalSeq(moPropen =list(p1, p2),
                       moMain = list(q1Main, q2Main), 
                       moCont = list(q1Cont, q2Cont),
                       fSet = list(fSet1,fSet2),
                       regimes = list(regime1, regime2),
                       data = data, response = data$Y, txName = c('A1', 'A2'),
                       Domains = rbind(c(200, 300), c(300,400)),
                       pop.size = 500, 
                       starting.values = c(250, 350))

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
