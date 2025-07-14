rm(list=ls())
library(MASS)
library(nleqslv)
library(mcmcplots)
library(bayesplot)

#Only the failures are considered
Data1 <- ## Stage 1 data
Data2 <- ## Stage 2 data

Data <- c(Data1, Data2)
n <- length(Data)
m <- 47
Tau <- 96


#Scheme1
R <- c(n-m, rep(0, m-1))
D_temp1 <- Data[-1]
store1 <- Data[1]
removed_sample <- sample(D_temp1, size = R[1],replace = FALSE)
censored_data <- D_temp1[! D_temp1  %in% c(removed_sample)]
xsamp <- c(store1, censored_data)
N1 <- length(xsamp[which(xsamp <Tau)])
N2 <- m-N1


# #Scheme2
# R=c(rep(1,n-m),rep(0, 2*m-n))
# D_temp2 <- Data[-1]
# store2 <- Data[1]
# removed_sample <-c(0)
# for (i in 1:(n-m)) {
#   removed_sample[i] <- sample(D_temp2, size = R[i],replace = FALSE)
#   xx <- D_temp2[! D_temp2 %in% removed_sample[i]]
#   store2 <- c(store2,xx[1])
#   D_temp2 <- xx[-1]
# }
# xsamp <- sort(c(store2,D_temp2))
# N1 <- length(xsamp[which(xsamp <Tau)])
# N2 <- m-N1


# ##Scheme3
# if(m%%2 == 0) {
#   R=c(rep(0, m/2 -1), n-m, rep(0, m/2))
# }else{
#   R=c(rep(0, (m+1)/2 -1), n-m, rep(0, (m-1)/2))
# }
# D_temp3 <- Data[-c(1:((m-1)/2))]
# store3 <- Data[1:((m-1)/2)]
# removed_sample <- sample(D_temp3, size = R[((m+1)/2)],replace = FALSE)
# censored_data <- D_temp3[! D_temp3  %in% c(removed_sample)]
# xsamp <- sort(c(store3,censored_data))
# N1 <- length(xsamp[which(xsamp <Tau)])
# N2 <- m-N1


#................Function for MLE................................
ml_est <- function(xsrt, cc, NN1, NN2, MM, RR){
  log_eq <- function(para){
    beta1 <- para[1]
    beta2 <- para[2]
    if(beta1 < 0 || beta2 < 0){
      beta1 <- xsrt[1]
      beta2 <- xsrt[2]
    }
    alp_trm1 <- 0; alp_trm2 <- 0;
    beta1_trm1 <- 0; beta1_trm2 <- 0;
    beta2_trm1 <- 0;
    
    for (i in 1:NN1) {
      alp_trm1 <- alp_trm1 + (1+RR[i])*(exp(beta1*cc[i])-1)
    }
    for (i in (NN1+1):MM) {
      alp_trm2 <- alp_trm2 + (1+RR[i]) * (exp(beta2*(cc[i]-Tau) + beta1*Tau)-1)
    }
    alp_trm <- MM / ( alp_trm1 + alp_trm2)
    
    for (k in 1:NN1) {
      beta1_trm1 <- beta1_trm1 + cc[k]* (1 - alp_trm*(1+RR[k])*exp(beta1*cc[k]))
    }
    for (k in (NN1+1):MM) {
      beta1_trm2 <- beta1_trm2 + Tau*(1 - alp_trm*(1+RR[k])*exp(beta2*(cc[k]-Tau) + beta1*Tau))
    }
    for (k in (NN1+1):MM) {
      beta2_trm1 <- beta2_trm1 + (cc[k]-Tau)*(1 - alp_trm*(1+RR[k])*exp(beta2*(cc[k]-Tau) + beta1*Tau))
    }
    
    eq <- numeric(2)
    eq[1] <- NN1/beta1 + beta1_trm1 + beta1_trm2
    eq[2] <- NN2/beta2 + beta2_trm1
    return(eq)
  } 
  sol <- nleqslv(xsrt, log_eq, method = "Newton")$x
  ml_bet1 <- sol[1]
  ml_bet2 <- sol[2]
  
  ml_alp_trm1 <- 0; ml_alp_trm2 <- 0;
  
  for (i in 1:NN1) {
    ml_alp_trm1 <- ml_alp_trm1 + (1+RR[i])*(exp(ml_bet1*cc[i])-1)
  }
  for (i in (NN1+1):MM) {
    ml_alp_trm2 <- ml_alp_trm2 + (1+RR[i]) * (exp(ml_bet2*(cc[i]-Tau) + ml_bet1*Tau)-1)
  }
  ml_alp <- MM / ( ml_alp_trm1 + ml_alp_trm2)
  
  return(c(ml_alp, ml_bet1, ml_bet2))
}
##Initial guess
xstart <- c(0.01,0.02)                     
mle <- ml_est(xstart, xsamp, N1, N2, m, R)

##.................Function for fisher information matrix....................
inf_mat <- function(ml, ss, NN1, NN2, MM, RR){
  a_22_trm1 <- 0; a_22_trm2 <- 0;
  a_33_trm1 <- 0;
  a_12_trm1 <- 0; a_12_trm2 <- 0;
  a_13_trm1 <- 0;
  a_23_trm1 <- 0
  
  for (i in 1:NN1) {
    a_22_trm1 <- a_22_trm1 + ml[1]*(1+RR[i])* (ss[i])^2 *exp(ml[2]*ss[i])
  }
  for (i in (NN1+1):MM) {
    a_22_trm2 <- a_22_trm2 + ml[1]*(1+RR[i])* Tau^2 *exp(ml[3]*(ss[i]-Tau) + ml[2]*Tau)
  }
  for (i in (NN1+1):MM) {
    a_33_trm1 <- a_33_trm1 + ml[1]*(1+RR[i])* (ss[i]-Tau)^2 * exp(ml[3]*(ss[i]-Tau) + ml[2]*Tau)
  }
  for (i in 1:NN1) {
    a_12_trm1 <- a_12_trm1 + (1+RR[i])*ss[i]*exp(ml[2]*ss[i])
  }
  for (i in (NN1+1):MM) {
    a_12_trm2 <- a_12_trm2 + Tau*(1+RR[i])*exp(ml[3]*(ss[i]-Tau) + ml[2]*Tau)
  }
  for (i in (NN1+1):MM) {
    a_13_trm1 <- a_13_trm1 + (1+RR[i])*(ss[i]-Tau)*exp(ml[3]*(ss[i]-Tau) + ml[2]*Tau)
  }
  for (i in (NN1+1):MM) {
    a_23_trm1 <- a_23_trm1 + ml[1]*Tau*(1+RR[i])*(ss[i]-Tau)*exp(ml[3]*(ss[i]-Tau) + ml[2]*Tau)
  }
  
  a11 <-  MM/ml[1]^2
  a22 <- NN1/ml[2]^2 + a_22_trm1 + a_22_trm2
  a33 <- NN2/ml[3]^2 + a_33_trm1
  a12 <- a_12_trm1 + a_12_trm2
  a21 <- a_12_trm1 + a_12_trm2  
  a13 <- a_13_trm1
  a31 <- a_13_trm1
  a23 <- a_23_trm1
  a32 <- a_23_trm1
  
  II <- matrix(c(a11,a12,a13,a21,a22,a23,a31,a32,a33), nrow = 3, ncol = 3, byrow = TRUE)
  return(II)
}
fish_inf_matrix <- inf_mat(mle, xsamp, N1, N2, m, R)





#....................Function for ACIs.................
aci_est <- function(ml, ss, vv_mat, NN1, NN2, MM, RR){
  
  #Inverse of the matrix
  var_mat <- solve(vv_mat)
  #95% confidence 
  alp_low  <- max(0, ml[1]-(1.96*sqrt(var_mat[1,1])))
  alp_upr  <- ml[1]+(1.96*sqrt(var_mat[1,1]))
  bet1_low <- max(0, ml[2]-(1.96*sqrt(var_mat[2,2])))
  bet1_upr <- ml[2]+(1.96*sqrt(var_mat[2,2]))
  bet2_low <- max(0, ml[3]-(1.96*sqrt(var_mat[3,3])))
  bet2_upr <- ml[3]+(1.96*sqrt(var_mat[3,3]))
  
  return(list(alp_low, alp_upr, bet1_low, bet1_upr, bet2_low, bet2_upr,var_mat))
}
aci <- aci_est(mle, xsamp, fish_inf_matrix, N1, N2, m, R)
aci
var_cov_matrix <- aci[[7]]
l1 <- aci[[2]]-aci[[1]]
l2 <- aci[[4]]-aci[[3]]
l3 <- aci[[6]]-aci[[5]]





#..............................function for MH_algorithm............................
MH_ieed <- function(cv_mt, ss, mll, itr, NN1, NN2, RR, MM, a11, b11, a22, b22, a33, b33, cc1, cc2){
  
  alpha <- betta1 <- betta2 <- c(0);
  alpha[1] <- mll[1]; betta1[1] <- mll[2]; betta2[1] <- mll[3];
  
  sample_gen <- function(itr){
    for(k in 1:itr){
      DD1 <- sum((1+RR[1:NN1])*(exp(betta1[k]*ss[1:NN1])-1))
      DD2 <- sum((1+RR[(NN1+1):MM])*(exp(betta2[k]*(ss[(NN1+1):MM]-Tau) + betta1[k]*Tau)-1))
      repeat{
        alp_temp <- rgamma(n = 1, shape = MM+a11, rate = b11 + DD1 + DD2 )
        if(alp_temp <  mle[1]*3){break}
      }
      repeat{
        bet1_temp <- rnorm(n = 1, mean = betta1[k] , sd = sqrt(cv_mt[2,2]))
        if(bet1_temp > 0 && bet1_temp < mle[2]*3 ){break}
      } 
      repeat{
        bet2_temp <- rnorm(n = 1, mean = betta2[k] , sd = sqrt(cv_mt[3,3]))
        if(bet2_temp > 0 && bet2_temp < mle[3]*3 ){break}
      }
      
      alpha[k+1] <- alp_temp
      ## Conditional posterior distribution of beta1
      post_distribution_bt1 <- function(xbet1){
        post_fun1 <-  xbet1^(NN1+a22-1)*exp(-b22*xbet1 + sum(xbet1*ss[1:NN1] -alp_temp*(1+RR[1:NN1])*exp(xbet1*ss[1:NN1])) + sum(xbet1*Tau - alp_temp*(1+RR[(NN1+1):MM])*exp(betta2[k]*(ss[(NN1+1):MM]-Tau) +xbet1*Tau))  )
        return(post_fun1)
      }
      ratio1 <- post_distribution_bt1(bet1_temp)/post_distribution_bt1(betta1[k])
      # ratio1[is.nan(ratio1)] <- 0.5
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
      # ratio2[is.nan(ratio2)] <- 0.5
      h2 <- min(1, ratio2)
      u2 <- runif(n = 1, min = 0, max = 1)
      ## Accept or reject the candidate
      if(u2 <= h2){
        betta2[k+1] <- bet2_temp
      } else{
        betta2[k+1] <- betta2[k]
      }
    }
    alpha_new  <- alpha            
    betta1_new <- betta1			
    betta2_new <- betta2 	
    return(list(alpha_new, betta1_new, betta2_new))
  }
  
  sample1 <- sample_gen(itr)
  sample2 <- sample_gen(itr)
  sample3 <- sample_gen(itr)
  samp1 <- cbind(sample1[[2]],sample1[[3]])
  samp2 <- cbind(sample2[[2]],sample2[[3]])
  samp3 <- cbind(sample3[[2]],sample3[[3]])
  colnames(samp1) <- paste(c('\u03b21', '\u03b22' ))
  colnames(samp2) <- paste(c('\u03b21', '\u03b22' ))
  colnames(samp3) <- paste(c('\u03b21', '\u03b22' ))
  samplesList <- list(chain1 = samp1, chain2 = samp2, chain3 = samp3)
  
  jpeg(filename = "hist plot data set 3.jpeg", width = 6, height = 2, units = "in", quality = 100, res = 650)
  color_scheme_set("viridis")
  print(mcmc_hist(samplesList))
  dev.off()
  
  jpeg(filename = "trace plot data set 3.jpeg", width = 6, height = 2, units = "in", quality = 100, res = 650)
  print(mcmc_trace(samplesList))
  dev.off()
  
  
  alpha_new  <- sample1[[1]][5001:itr]
  betta1_new <- sample1[[2]][5001:itr]
  betta2_new <- sample1[[3]][5001:itr]
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
  
  ##Next, find 95% Credible Interval
  alp_sort  <- sort(alpha_new)
  bet1_sort <- sort(betta1_new)
  bet2_sort <- sort(betta2_new)
  
  
  bay_alp_low  <- bay_alp_upp  <- bay_alp_len <- c(0)
  bay_bet1_low <- bay_bet1_upp <- bay_bet1_len <- c(0)
  bay_bet2_low <- bay_bet2_upp <- bay_bet2_len <- c(0)
  
  for(j in 1:(los*(length(alpha_new)))){
    bay_alp_low[j]  <- alp_sort[j]
    bay_alp_upp[j]  <- alp_sort[j + length(alpha_new)*(1-los)]
    bay_alp_len[j]  <- bay_alp_upp[j] - bay_alp_low[j]
    bay_bet1_low[j] <- bet1_sort[j]
    bay_bet1_upp[j] <- bet1_sort[j + length(betta1_new)*(1-los)]
    bay_bet1_len[j] <- bay_bet1_upp[j] - bay_bet1_low[j]
    bay_bet2_low[j] <- bet2_sort[j]
    bay_bet2_upp[j] <- bet2_sort[j + length(betta2_new)*(1-los)]
    bay_bet2_len[j] <- bay_bet2_upp[j] - bay_bet2_low[j]
  }
  min_diff_alp  <- which(bay_alp_len  == min(bay_alp_len))[1]
  min_diff_bet1 <- which(bay_bet1_len == min(bay_bet1_len))[1]
  min_diff_bet2 <- which(bay_bet2_len == min(bay_bet2_len))[1]
  
  halp_low  <- bay_alp_low[min_diff_alp]
  halp_upp  <- bay_alp_upp[min_diff_alp]
  hbet1_low <- bay_bet1_low[min_diff_bet1]
  hbet1_upp <- bay_bet1_upp[min_diff_bet1]
  hbet2_low <- bay_bet2_low[min_diff_bet2]
  hbet2_upp <- bay_bet2_upp[min_diff_bet2]
  #return(c(lam_SE, alp1_SE, alp2_SE, lam_LI1, alp1_LI1, alp2_LI1, lam_LI2, alp1_LI2, alp2_LI2, hlam_low, hlam_upp, halp1_low, halp1_upp, halp2_low, halp2_upp))  
  return(c(alpha_SE, betta1_SE, betta2_SE, alp_LI1, bet1_LI1, bet2_LI1, alp_LI2, bet1_LI2, bet2_LI2, halp_low, halp_upp, hbet1_low, hbet1_upp, hbet2_low,hbet2_upp))
}
a1= 0.001; a2 = 0.001; a3 = 0.001; b1= 0.001; b2 = 0.001; b3 = 0.001;
itrr = 30000; 
los = 0.05; c1 = -0.5; c2 = 0.5;
bayes_NIP   <- MH_ieed(var_cov_matrix, xsamp, mle, itrr, N1, N2, R, m, a1, b1, a2, b2, a3, b3, c1, c2)
bayes_est <- bayes_NIP[1:9]
bayes_est_intrvl <- bayes_NIP[10:15]
b_ap_int <- bayes_est_intrvl[2]-bayes_est_intrvl[1]
b_bt1_int <- bayes_est_intrvl[4]-bayes_est_intrvl[3]
b_bt2_int <- bayes_est_intrvl[6]-bayes_est_intrvl[5]
ml_int <- round(c(l1,l2,l3), 5)
bs_int <-round(c(b_ap_int,b_bt1_int,b_bt2_int), 5)
round(mle,5)
round(bayes_est,5)
ml_int
bs_int



