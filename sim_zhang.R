#######Replication of simulations from Zhang at al ########
###https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3843953/###
###########################################################
library('DynTxRegime')
rm(list = ls())

##HIV patients randomised to 1 (start) or 0 (not start)
##Time points baseline, 6 months 

expit <- function(x) { exp(x) / (1 + exp(x)) } #inverse logit

n=500 #number of patients 

#at baseline 
C0 = rnorm(n, mean = 450, sd = 150)
A1 = expit(2-0.006*C0) #first decision depends on baseline CD4

#at 6 months 
C1 = rnorm (n, mean = 1.25 * C0, sd = 60) #CD4 improves with time ?? how is it dependent on A1
A2 = A1 + (1 - A1) * expit(0.8 - 0.004 * C1) # second decision depends on 6-month CD4

#at 1 year 
C2 = 400 + 1.6 * C0 - abs(250 - C0)*(A1 - ifelse((250-C0)>0,1,0))**2 - (1 - A1) * abs(720 - 2 * C1) * (A2 - ifelse((720 - 2*C1) > 0, 1, 0))**2
Y = rnorm(n, mean = C2, sd = 60)

### 400 + 1.6 * C0 -> CD4 improves from baseline after 1 year if subsequent terms are 0 
### first decision: 
### -abs(250 - C0)*(A1 - ifelse((250-C0)>0,1,0))**2  -> when C0 < 250, A1 = 1, whole term = 0 
#                                                                      A1= 0, whole term = -(250 − c1) 
#                                                       when C0 > 250, A1= 1, whole term = -(250 − c1)
#                                                                      A1= 0, whole term = 0
# At first decision, if start treatment when C0<250, no benefit, withholding treatment -(250 − c1) 
#                    if start treatment when C0>250, -(250 − c1), withholding treament no benefit 

###second decision: 
### -abs(720 - 2 * C1) * (A2 - ifelse((720 - 2*C1) > 0, 1, 0))**2 -> when C0 < 360, A1 = 1, whole term = 0 
#                                                                                   A1= 0, whole term = -(720 - 2*c1)
#                                                                    when C0 > 360, A1= 1, whole term = -(720 - 2*c1)
#                                                                                   A1= 0, whole term = 0
# At first decision, if start treatment when C0<360, no benefit, withholding treatment -(720 - 2*c1) 
#                    if start treatment when C0>360, -(720 - 2*c1), withholding treament no benefit 

#true Q contrast functions - 
## cont2(c0, c1, a1) = (1 − a1) (720 − 2*C1)
## cont1(c0) = 250 − c1 
# the optimal treatment regime gopt=(gopt1,gopt2) 
## gopt1(c1)= ifelse((250−c1>0), 1, 0)
## gopt2(c2,a1)=I{a1+(1−a1)(720−*c1)>0} and =I{a1+(1−a1)(360−c1)>0}
## E{Y*(gopt)} = 1120 ???

d = cbind.data.frame(id = 1:n, #ID 
                     C0 = C0,
                     C1 = C1, 
                     A1 = A1, 
                     A2 = A2, 
                     Y = Y)  #Outcome
d$A1=ifelse(d$A1==1, 'start', 'hold')
d$A2=ifelse(d$A2==1, 'start', 'hold')


#################################################################################
####################################Analysis#####################################
#################################################################################

# Scenario 1 - Q functions are misspecified 
###########################################

moPropen1 <- buildModelObj(model = ~ C1,
                           solver.method = 'glm',
                           solver.args = list('family'='binomial'),
                           predict.method = 'predict.glm',
                           predict.args = list(type='response'))
moPropen2 <- buildModelObj(model = ~ C0,
                           solver.method = 'glm',
                           solver.args = list('family'='binomial'),
                           predict.method = 'predict.glm',
                           predict.args = list(type='response'))
moPropen <- list(moPropen1, moPropen2)

# outcome model second stage
moMain2 <- buildModelObj(model = ~ C1,
                         solver.method = 'lm')
moCont2 <- buildModelObj(model = ~ C0,
                         solver.method = 'lm')

# outcome model first stage
moMain1 <- buildModelObj(model = ~ C1,
                         solver.method = 'lm')
moCont1 <- buildModelObj(model = ~ A1,
                         solver.method = 'lm')

moMain <- list(moMain1, moMain2)
moCont <- list(moCont1, moCont2)

# regime function second stage
regime2 <- function(eta2, data) {
  tst <- {data$C2 >= eta2}
  rec <- rep('start', nrow(x = data))
  rec[!tst] <- 'hold'
  return( rec )
}
# regime function first stage
regime1 <- function(eta1, data) {
  tst <-  {data$C1 >= eta1}
  rec <- rep('start', nrow(x = data))
  rec[!tst] <- 'hold'
  return( rec )
}
regimes <- list(regime1, regime2)

#### Analysis using AIPW
fit_AIPW <- optimalSeq(moPropen = moPropen,
                       moMain = moMain, moCont = moCont,
                       regimes = regimes,
                       data = d, response = Y, txName = c('A1','A2'),
                       Domains = cbind(rep(100,2),rep(800,2)),
                       pop.size = n, starting.values = rep(200,2))

#estimated values for eta 
regimeCoef(object = fit_AIPW)

# Scenario 2 - models for Q-contrast functions are misspecified 
###############################################################

moPropen1 <- buildModelObj(model = ~ C1,
                           solver.method = 'glm',
                           solver.args = list('family'='binomial'),
                           predict.method = 'predict.glm',
                           predict.args = list(type='response'))
moPropen2 <- buildModelObj(model = ~ C0,
                           solver.method = 'glm',
                           solver.args = list('family'='binomial'),
                           predict.method = 'predict.glm',
                           predict.args = list(type='response'))
moPropen <- list(moPropen1, moPropen2)

# outcome model second stage
moMain2 <- buildModelObj(model = ~ C1,
                         solver.method = 'lm')
moCont2 <- buildModelObj(model = ~ C0,
                         solver.method = 'lm')

# outcome model first stage
moMain1 <- buildModelObj(model = ~ C1,
                         solver.method = 'lm')
moCont1 <- buildModelObj(model = ~ A1,
                         solver.method = 'lm')

moMain <- list(moMain1, moMain2)
moCont <- list(moCont1, moCont2)

# regime function second stage
regime2 <- function(eta2, data) {
  tst <- {data$C2 >= eta2}
  rec <- rep('start', nrow(x = data))
  rec[!tst] <- 'hold'
  return( rec )
}
# regime function first stage
regime1 <- function(eta1, data) {
  tst <-  {data$C1 >= eta1}
  rec <- rep('start', nrow(x = data))
  rec[!tst] <- 'hold'
  return( rec )
}
regimes <- list(regime1, regime2)

#### Analysis using AIPW
fit_AIPW <- optimalSeq(moPropen = moPropen,
                       moMain = moMain, moCont = moCont,
                       regimes = regimes,
                       data = d, response = Y, txName = c('A1','A2'),
                       Domains = cbind(rep(100,2),rep(800,2)),
                       pop.size = n, starting.values = rep(200,2))

#estimated values for eta 
regimeCoef(object = fit_AIPW)
