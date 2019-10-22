########Simulation of observational data for treatment duration of VAP#########
###############################################################################

expit <- function(x) { exp(x)/(1+exp(x)) } #inverse logit 

#Parameters
n = 500      #number of patients
t_observe=2  #number of points of data 

coef_age= 0.02
coef_comorb= 0.03
coef_icu= 0.015

coef_map = -0.02
coef_tw = 0.035
coef_fever = 0.04

#.......................................................#
#....................Data generation....................#
#.......................................................#

#patient id
id = seq(1:n)

#covariates
###baseline 
age = rnorm(n, mean = 65, sd = 5)                      # age (continuous)
comorb = round(rnorm(n, mean = 3, sd = 1))             # comobidities score (categorical)
icu = sample(rep(c('medical', 'surgical'), n/2))       # medical or surgical ICU (binary)
  
###time varying 
######## TW (continuous)

gen_tw=function(t_observe, n=n){ #t_observe refers to the total observational period 
  
  t_df=tw_df=matrix(NA,ncol=t_observe, nrow = n)
  
  # 1st test 
  t_df[,1] = round(rnorm(n, mean=0))
  
  for (i in 2:t_observe) { #matrix for time of tests with reference to VAP start date
    t_test = ceiling(rexp(n, 1/(i+1))) #time between tests are from an exponential distribution
    t_df[,i] =  t_df[, (i-1)] + t_test 
  }
  colnames(t_df)=paste0('day_tw_',1:t_observe)
  
  for (i in 1:t_observe) { #matris for TW 
    tw_df[,i] = abs(rnorm(n, mean = (15-0.2*i), sd=3)) #TW improves with time 
  }
  colnames(tw_df)=paste0('tw_',1:t_observe)
  
  df=cbind.data.frame(t_df,tw_df)

  
  return(df)
  
}

tw = gen_tw(t_observe = t_observe, n = n)

######## MAP (categorical)

gen_map=function(t_observe, n=n){ #t_observe refers to the total observational period 
  
  map_df=matrix(NA,ncol=t_observe, nrow = n)
  
  map_df[,1] = abs(round(rnorm(n, mean=70, sd=10)/10))
 
  for (i in 2:t_observe) { #matrix for time of tests with reference to VAP start date
    map_df[,i] =  map_df[, (i-1)] + abs(round(rnorm(n,mean=0.25))) #blood pressure improves with time 
  }
  colnames(map_df)=paste0('map_',1:t_observe)
  
  t_df=matrix(1:t_observe,ncol=t_observe, byrow = T, nrow = n)
  colnames(t_df)=paste0('day_map_',1:t_observe)
  
  df=cbind.data.frame(t_df, map_df)
  
  return(df)
  
}

map = gen_map(t_observe = t_observe, n = n)

######## presence of fever > 38.3°C  (binary)

gen_fever=function(t_observe, n=n){ #t_observe refers to the total observational period 
  
  fever_df=matrix(NA,ncol=t_observe, nrow = n)
  
  for (i in 1:t_observe) { #matrix for time of tests with reference to VAP start date
    prob= ifelse((0.8-i*0.05)>0, (0.8-i*0.05), 0) #probability of fever decreased with time 
    fever_df[,i] = rbinom(n, 1, prob)
  }
  colnames(fever_df) = paste0('fever_',1:t_observe)
  
  t_df=matrix(1:t_observe,ncol=t_observe, byrow = T, nrow = n)
  colnames(t_df)=paste0('day_fever_',1:t_observe)
  
  df=cbind.data.frame(t_df, fever_df)
  
  return(df)
  
}

fever=gen_fever(t_observe = t_observe, n = n)

#treatment - culture appropriate antibiotics 

#### A1, A2 - binary indicator for continuing antibiotics at day 3 and day 7 respectively 
####Probability to continue antibiotics at point 1 (day 3)
prob.continue.1 = expit(0 + coef_age * age + coef_comorb * comorb + coef_icu * as.numeric(as.factor(icu)) + 
                        coef_map * map[,3] + coef_tw * tw[,3] + coef_fever * fever[,3])

  
#outcome: mortality at day 60 
prob = 
y = rbinom(n, 1, prob)

data = cbind.data.frame(id, 
                        age, comorb, icu, 
                        tw, map, fever, 
                        dur, y)
  