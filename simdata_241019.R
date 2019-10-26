########Simulation of observational data for treatment duration of VAP#########
###############################################################################

# Based on the simulation setting used in Zhang et al. 2013 https://academic.oup.com/biomet/article-abstract/100/3/681/303040, Section 5

expit <- function(x) { exp(x) / (1 + exp(x)) } #inverse logit
set.seed(1234)

#Parameters
n = 1000    #number of patients

#.......................................................#
#....................Data generation....................#
#.......................................................#

##Baseline
A0 = rep(1, n)                                        # Start treatment for everyone
# L =  rbeta(n, shape1 = 2, shape2 = 2)                  # baseline characteristics: comobidities score
# C0 = rbeta(n, shape1 = 3, shape2 = 2) + 0.15 * L        # Time varying covariate: disease severity score
L =  rnorm(n, mean = 1, sd = 0.5)                  # baseline characteristics: comobidities score
C0 = rnorm(n, mean = 2, sd = 0.2) + 0.15 * L        # Time varying covariate: disease severity score

##Time varying, t=1 
# C1 depends on L, C0, A0
mu_C1 = C0 - 0.1 + 0.2 * L - 0.35 * A0
C1 = rnorm(n, mu_C1, 0.1)

# P(A1) depends on L, C0, C1
P_A1 = expit(0.3 + 0.2 * L + 0.3 * C0 + 0.7 * C1)
A1 = rbinom(n , 1, P_A1)

##Time varying, t=2
# C2 depends on L, C1
mu_C2 = C1 - 0.05 + 0.2 * L - 0 * A1 #does not depend on previous A1
C2 = rnorm(n, mu_C2, 0.1)

# P(A2) depends on L, C2
P_A2 = expit(0.2 * L + 0.4 * C2)
A2 = rbinom(n , 1, P_A2)
A2[which(A1 == 0)] = 0 #those who had antibiotics stopped at t=1 should have 0 in t=2

##Regime = benefits survival if at t=1, L>4 C1>6 continues antibiotics and 
#                               at t=2, L>4 C2>8 continues antibiotics

# R = rep(0, n)
# R[intersect(which(L> 0.4 & C1> 0.6 & A1 ==1), which(L> 0.4 & C2> 0.8 & A2 ==1))] = 1

##Outcome
# P_Y = expit(0.5 * L + 1 * R + 0.1 * C0 + 0.1 * C1 + 0.2 * C2) # A2 has no effect on outcome
# Y = rbinom(n, 1, P_Y)

## Note: Can the outcome be a continuous variable? Having a discrete variable complicates analysis. 
## Have to consider survival models otherwise.

## Note
#  Objective: Minimize Y and minimize time on treatment

# Assuming continuous Y, say, disease severity score at t=10
# Y is generated such that optimal regime is the following:
# If C1<=1.55 -> Stop treatment i.e. A1=0,A2=0,
# If C1>1.55 & C2<=2.55 -> Stop treatment i.e. A2=0, and
# If C1>1.55 & C2>2.55 -> Keep treatment
# Y depends on C1, A1, C2, A2
# Y = 2 + 0.6*C1 + (1-A1)*|C1 - 1.5|*{1(C1 - 1.5 > 0) - A1}^2 + A1*(1-A2)*|C2 - 2.5|*{1(C2 - 2.5 > 0) - A2}^2\
#     + 0.05*A1 + 0.05*A1*A2
mu_Y = 2 + 0.6 * C1 + (1 - A1) * abs(C1-1.5) * (ifelse((C1-1.5)>0,1,0) - A1)**2 + A1 * (1 - A2) * abs(C2-2.5) * (ifelse(C2-2.5>0,1,0) - A2)**2 + 0.1*A1 + 0.1*A1*A2
Y = rnorm(n, mu_Y, sd=0.1)

# Optimal regime
# g1 = 1(C1 - 1.5 > 0), g2 = 1(-1 + A1 + A1*(C2-2.5) > 0)

##Data
d = cbind.data.frame(id = 1:n, #ID 
                     L = L,    #Baseline covariates
                     C0 = C0,  #Time varying covariates
                     C1 = C1, 
                     C2 = C2, 
                     A0 = A0,  #Treatment decisions
                     A1 = A1, 
                     A2 = A2, 
                     Y = Y)  #Outcome

mean(Y[which(A1==1 & A2 ==1)]) #mortality in the regime 111
mean(Y[which(A1==1 & A2 ==0)]) #mortality in the regime 110
mean(Y[which(A1==0 & A2 ==0)]) #mortality in the regime 100

## Using R package
library(DynTxRegime)

# Outcome model second stage
moMain2 <- buildModelObj(model=~C1+C2+A1, solver.method='lm')
moCont2 <- buildModelObj(model=~C2, solver.method='lm')
moPropen <- buildModelObj(model=~1, solver.method='glm',
                          solver.args=list('family'='binomial'),
                          predict.method='predict.glm',
                          predict.args=list(type='response'))

# Outcome model first stage
moMain1 <- buildModelObj(model=~C1+L, solver.method='lm')
moCont1 <- buildModelObj(model=~C1, solver.method='lm')

moPropen <- list(moPropen, moPropen)
moMain <- list(moMain1, moMain2)
moCont <- list(moCont1, moCont2)

# Regime function second stage
regime2 <- function(eta1, eta2, data){
  tst <- {data$C1 > eta1} & {data$C2 <= eta2}
  rec <- rep('0', nrow(x=data))
  rec[!tst] <- '1'
  return(rec)
}

# Regime function first stage
regime1 <- function(eta1, data){
  tst <- {data$C1 <= eta1}
  rec <- rep('0', nrow(x=data))
  rec[!tst] <- '1'
  return(rec)
}

d$A1 <- as.factor(d$A1)
d$A2 <- as.factor(d$A2)

regimes <- list(regime1, regime2)

fit_AIPW <- optimalSeq(moPropen=moPropen,
                       moMain=moMain, moCont=moCont,
                       regimes=regimes,
                       data=d, response=Y, txName=c('A1', 'A2'),
                       Domains=cbind(rep(1,3), rep(4,3)),
                       pop.size=n, starting.values=rep(2.5,3))
