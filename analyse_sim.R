source('gaussian_sim_data.R')

library(DynTxRegime)

data = as.data.frame(generate(n = 500, seq.rand = 1, psi10 = 250, psi20 = 720, psi11=-1, psi21 = -2))

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

# regime function second stage
regime2 <- function(eta2, data){return(data$A1 + {1L-data$A1}*{data$L2 < eta2})}
# regime function first stage
regime1 <- function(eta1, data){return(as.integer(x = {data$L1 < eta1})) }

out = optimalSeq(moPropen =list(p1, p2),
           moMain = NULL, 
           moCont = NULL,
           fSet = list(fSet1,fSet2),
           regimes = list(regime1, regime2),
           data = data, response = data$Y, txName = c('A1', 'A2'),
           Domains = rbind(c(50, 500), c(50, 500)),
           pop.size = 2000, 
           starting.values = c(250, 350))
