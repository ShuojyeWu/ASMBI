### Inference for a simple step-stress model with progressively censored data from GD ###
rm(list = ls())
library(MASS)
library(nleqslv)
library(matlib)
#----------------------------------------------------
## Intilization of parameter
alpp = 0.1
bett = c(0.3, 0.9)
n = 60
m = 55
Tau = 6.5
b1 = 1.7
b2 = 0.5
b3 = 0.9
itrr = 10000
los = 0.05
c1 = -0.5
c2 = 0.5
a1 = b1 * alpp
a2 = b2 * bett[1]
a3 = b3 * bett[2]
N = 1001
#----------------------------------------------------
## Censoring Schemes
R = c(n - m, rep(0, m - 1))

# R=c(rep(1,n-m),rep(0, 2*m-n))

# if(m%%2 == 0) {
#   R=c(rep(0, m/2 -1), n-m, rep(0, m/2))
# }else{
#   R=c(rep(0, (m+1)/2 -1), n-m, rep(0, (m-1)/2))
# }
#----------------------------------------------------
## Progressive Sample generation
sample_progressive <- function(nn, mm, alp, bet, RR) {
  ww = runif(mm)
  vv = numeric(0)
  uu = numeric(0)
  for (i1 in 1:mm) {
    vv[i1] <- (ww[i1])^(1 / (i1 + sum(RR[(mm - i1 + 1):mm])))
  }
  for (i2 in 1:mm) {
    uu[i2] <- 1 - prod(vv[(mm - i2 + 1):mm])
    #These ui's are progressive Type-2 sample from the U[0,1] distribution
  }
  F_Tau <- 1 - exp(-alp * (exp(bet[1] * Tau) - 1))
  n1 <- sum(uu < F_Tau)
  n2 <- m - n1
  xx <- c(0)
  for (j1 in 1:n1) {
    xx[j1] <- 1 / bet[1] * log(1 - log(1 - uu[j1]) / alp)    # for 0 < t < tau
  }
  for (j2 in (n1 + 1):m) {
    xx[j2] <- 1 / bet[2] * log(1 - log(1 - uu[j2]) / alp) - Tau * (bet[1] /
                                                                     bet[2]) + Tau    # for t >= tau
  }
  return(list(xx, n1, n2))
}
ALT_sample <- sample_progressive(n, m, alpp, bett, R)
xsamp <- ALT_sample[[1]]
N1 <- ALT_sample[[2]]
N2 <- ALT_sample[[3]]
#------------------------------------------------------------
# Maximum likelihood estimates
ml_est <- function(xsrt, cc, NN1, NN2, MM, RR) {
  log_eq <- function(para) {
    beta1 <- para[1]
    beta2 <- para[2]
    if (beta1 < 0 || beta2 < 0) {
      beta1 <- xsrt[1]
      beta2 <- xsrt[2]
    }
    alp_trm1 <- 0
    alp_trm2 <- 0
    beta1_trm1 <- 0
    beta1_trm2 <- 0
    beta2_trm1 <- 0
    for (i in 1:NN1) {
      alp_trm1 <- alp_trm1 + (1 + RR[i]) * (exp(beta1 * cc[i]) - 1)
    }
    for (i in (NN1 + 1):MM) {
      alp_trm2 <- alp_trm2 + (1 + RR[i]) * (exp(beta2 * (cc[i] - Tau) + beta1 *
                                                  Tau) - 1)
    }
    alp_trm <- MM / (alp_trm1 + alp_trm2)
    
    for (k in 1:NN1) {
      beta1_trm1 <- beta1_trm1 + cc[k] * (1 - alp_trm * (1 + RR[k]) * exp(beta1 *
                                                                            cc[k]))
    }
    for (k in (NN1 + 1):MM) {
      beta1_trm2 <- beta1_trm2 + Tau * (1 - alp_trm * (1 + RR[k]) * exp(beta2 *
                                                                          (cc[k] - Tau) + beta1 * Tau))
    }
    for (k in (NN1 + 1):MM) {
      beta2_trm1 <- beta2_trm1 + (cc[k] - Tau) * (1 - alp_trm * (1 + RR[k]) *
                                                    exp(beta2 * (cc[k] - Tau) + beta1 * Tau))
    }
    eq <- numeric(2)
    eq[1] <- NN1 / beta1 + beta1_trm1 + beta1_trm2
    eq[2] <- NN2 / beta2 + beta2_trm1
    return(eq)
  }
  sol <- nleqslv(xsrt, log_eq, method = "Newton")$x
  ml_bet1 <- sol[1]
  ml_bet2 <- sol[2]
  ml_alp_trm1 <- 0
  ml_alp_trm2 <- 0
  for (i in 1:NN1) {
    ml_alp_trm1 <- ml_alp_trm1 + (1 + RR[i]) * (exp(ml_bet1 * cc[i]) - 1)
  }
  for (i in (NN1 + 1):MM) {
    ml_alp_trm2 <- ml_alp_trm2 + (1 + RR[i]) * (exp(ml_bet2 * (cc[i] - Tau) + ml_bet1 *
                                                      Tau) - 1)
  }
  ml_alp <- MM / (ml_alp_trm1 + ml_alp_trm2)
  
  return(c(ml_alp, ml_bet1, ml_bet2))
}
xstart <- c(bett[1] - .1, bett[2] - .1)
mle <- ml_est(xstart, xsamp, N1, N2, m, R)
#-----------------------------------------------------------------
inf_mat <- function(ml, ss, NN1, NN2, MM, RR) {
  a_22_trm1 <- 0
  a_22_trm2 <- 0
  a_33_trm1 <- 0
  a_12_trm1 <- 0
  a_12_trm2 <- 0
  a_13_trm1 <- 0
  a_23_trm1 <- 0
  for (i in 1:NN1) {
    a_22_trm1 <- a_22_trm1 + ml[1] * (1 + RR[i]) * (ss[i])^2 * exp(ml[2] * ss[i])
  }
  for (i in (NN1 + 1):MM) {
    a_22_trm2 <- a_22_trm2 + ml[1] * (1 + RR[i]) * Tau^2 * exp(ml[3] * (ss[i] -
                                                                          Tau) + ml[2] * Tau)
  }
  for (i in (NN1 + 1):MM) {
    a_33_trm1 <- a_33_trm1 + ml[1] * (1 + RR[i]) * (ss[i] - Tau)^2 * exp(ml[3] *
                                                                           (ss[i] - Tau) + ml[2] * Tau)
  }
  for (i in 1:NN1) {
    a_12_trm1 <- a_12_trm1 + (1 + RR[i]) * ss[i] * exp(ml[2] * ss[i])
  }
  for (i in (NN1 + 1):MM) {
    a_12_trm2 <- a_12_trm2 + Tau * (1 + RR[i]) * exp(ml[3] * (ss[i] - Tau) + ml[2] *
                                                       Tau)
  }
  for (i in (NN1 + 1):MM) {
    a_13_trm1 <- a_13_trm1 + (1 + RR[i]) * (ss[i] - Tau) * exp(ml[3] * (ss[i] -
                                                                          Tau) + ml[2] * Tau)
  }
  for (i in (NN1 + 1):MM) {
    a_23_trm1 <- a_23_trm1 + ml[1] * Tau * (1 + RR[i]) * (ss[i] - Tau) * exp(ml[3] *
                                                                               (ss[i] - Tau) + ml[2] * Tau)
  }
  a11 <-  MM / ml[1]^2
  a22 <- NN1 / ml[2]^2 + a_22_trm1 + a_22_trm2
  a33 <- NN2 / ml[3]^2 + a_33_trm1
  a12 <- a_12_trm1 + a_12_trm2
  a21 <- a_12_trm1 + a_12_trm2
  a13 <- a_13_trm1
  a31 <- a_13_trm1
  a23 <- a_23_trm1
  a32 <- a_23_trm1
  II <- matrix(
    c(a11, a12, a13, a21, a22, a23, a31, a32, a33),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )
  return(II)
}
var_cov_matrix <- inf_mat(mle, xsamp, N1, N2, m, R)
#----------------------------------------------------------------------------------
#MH_algorithm
MH_ieed <- function(cv_mt, ss, mll, itr, NN1, NN2, RR, MM, a11, b11, a22, b22, a33, b33, cc1, cc2){
  alpha <- betta1 <- betta2 <- c(0);
  alpha[1] <- mll[1]; betta1[1] <- mll[2]; betta2[1] <- mll[3];
  for(k in 1:itr){
    DD1 <- sum((1+RR[1:NN1])*(exp(betta1[k]*ss[1:NN1])-1))
    DD2 <- sum((1+RR[(NN1+1):MM])*(exp(betta2[k]*(ss[(NN1+1):MM]-Tau) + betta1[k]*Tau)-1))
    alp_temp <- rgamma(n = 1, shape = MM+a11, rate = b11 + DD1 + DD2 )
    repeat{
      bet1_temp <- rnorm(n = 1, mean = betta1[k] , sd = sqrt(cv_mt[2,2]))
      if(bet1_temp > 0 && bet1_temp < bett[1]*2.5){break}
    } 
    repeat{
      bet2_temp <- rnorm(n = 1, mean = betta2[k] , sd = sqrt(cv_mt[3,3]))
      if(bet2_temp > 0 && bet2_temp < bett[2]*2.5){break}
    }
    alpha[k+1] <- alp_temp
    ## Conditional posterior distribution of beta1
    post_distribution_bt1 <- function(xbet1){
      post_fun1 <-  xbet1^(NN1+a22-1)*exp(-b22*xbet1 + sum(xbet1*ss[1:NN1] -alp_temp*(1+RR[1:NN1])*exp(xbet1*ss[1:NN1]))
                                          + sum(xbet1*Tau - alp_temp*(1+RR[(NN1+1):MM])*exp(betta2[k]*(ss[(NN1+1):MM]-Tau) +xbet1*Tau))  )
      return(post_fun1)
    }
    ratio1 <- post_distribution_bt1(bet1_temp)/post_distribution_bt1(betta1[k])
    ratio1[is.nan(ratio1)] <- 0.5
    h1 <- min(1, ratio1)
    u1 <- runif(n = 1, min = 0, max = 1)
    ## Accept or reject the candidate
    if(u1 <= h1){
      betta1[k+1] <- bet1_temp
    } else{
      betta1[k+1] <- betta1[k]
    }
    ## Conditional posterior distribution of beta2
    post_distribution_bt2 <- function(xbet2){
      post_fun2 <- xbet2^(NN2 + a33 -1)*exp(-b33*xbet2   + sum(xbet2*(ss[(NN1+1):MM]-Tau) -alp_temp*(1+RR[(NN1+1):MM])*exp(xbet2*(ss[(NN1+1):MM]-Tau) +betta1[k+1]*Tau)))
      return(post_fun2)
    }
    ratio2 <- post_distribution_bt2(bet2_temp)/post_distribution_bt2(betta2[k])	
    ratio2[is.nan(ratio2)] <- 0.5
    h2 <- min(1, ratio2)
    u2 <- runif(n = 1, min = 0, max = 1)
    ## Accept or reject the candidate
    if(u2 <= h2){
      betta2[k+1] <- bet2_temp
    } else{
      betta2[k+1] <- betta2[k]
    }
  }
  alpha_new  <- alpha[2001:10000]            ##Burn the first 2000 sample from the generated posterior sample of alpha
  betta1_new <- betta1[2001:10000]					 ##Burn the first 2000 sample from the generated posterior sample of betta1
  betta2_new <- betta2[2001:10000]					 ##Burn the first 2000 sample from the generated posterior sample of betta2
  ## Estimation under Squared error loss function
  alpha_SE  <- mean(alpha_new)
  betta1_SE <- mean(betta1_new)
  betta2_SE <- mean(betta2_new) 
  # ## Estimation under LINEX loss function with positive value of parameter c
  alp_LI1  <- -(1/cc1)*log(mean(exp(-cc1*alpha_new)))
  bet1_LI1 <- -(1/cc1)*log(mean(exp(-cc1*betta1_new)))
  bet2_LI1 <- -(1/cc1)*log(mean(exp(-cc1*betta2_new)))
  # ## Estimation under LINEX loss function with negative value of parameter c
  alp_LI2  <- -(1/cc2)*log(mean(exp(-cc2*alpha_new)))
  bet1_LI2 <- -(1/cc2)*log(mean(exp(-cc2*betta1_new)))
  bet2_LI2 <- -(1/cc2)*log(mean(exp(-cc2*betta2_new)))
  return(c(alpha_SE, betta1_SE, betta2_SE, alp_LI1, bet1_LI1, bet2_LI1, alp_LI2, bet1_LI2, bet2_LI2))
}
bayes_IP   <- MH_ieed(var_cov_matrix, xsamp, mle, itrr, N1, N2, R, m, a1, b1, a2, b2, a3, b3, c1, c2)
bayes_IP_Se <- bayes_IP[1:3]
bayes_IP_LT1 <- bayes_IP[4:6]
bayes_IP_LT2 <- bayes_IP[7:9]
mle
bayes_IP_Se
bayes_IP_LT1
bayes_IP_LT2