rm(list=ls(all=TRUE));
library(nleqslv)

#Only the failures are considered
Data1 <- sort(c(12.07, 19.50, 22.10, 23.11, 24.00, 25.10, 26.90,
                36.64, 44.10, 46.30, 54.00, 58.09, 64.17, 72.25,
                86.90, 90.09, 91.22, 14.00, 17.95, 24.46, 26.46,
                26.58, 28.06, 34.00, 36.13, 40.85, 41.11, 42.63,
                52.51, 62.68, 72.13, 83.63, 91.56, 94.38))

Data2 <- sort(c(102.10, 105.10, 109.20, 114.90, 117.90, 121.90, 122.50, 123.60, 126.50, 
                130.10, 97.71, 101.13, 105.11, 112.11, 119.58, 120.20, 126.95, 129.25, 136.31))
Tau <- 96



Data <- Data1
n <- length(Data)
##mle complete data under normal stress
fn <- function(bet){
  term1 <- n / sum(exp(bet*Data) - 1)
  lam_fn <- (n /bet) + sum(Data) - term1 *sum(Data*exp(bet*Data))
  return(lam_fn)
}
bet_val <- nleqslv(1, fn)
bet_mle <- bet_val$x
alpp_mle <- n / sum(exp(bet_mle*Data) - 1)
## Let us define the CDF of GD under normal stress
cdf_gd1 <- function(x){
  cdf <- 1 - exp(-alpp_mle*(exp(bet_mle*x) -1))
  return(cdf)
}
# Apply KS test
ks1 = ks.test(Data,cdf_gd1)
ks1





# Data <- Data2
# n <- length(Data)
# ##mle complete data under high stress
# fn <- function(para){
#   beta1 <- para[1]
#   beta2 <- para[2]
#   term1 <- n / sum(exp((Data-Tau)*beta2 + Tau*beta1) -1)
#   eq <- numeric(2)
#   eq1 <- sum(Tau*beta1 -term1*Tau*exp((Data-Tau)*beta2 + Tau*beta1) )
#   eq2 <- n/beta2 + sum((Data-Tau) - term1*(Data-Tau)*exp((Data-Tau)*beta2 + Tau*beta1) )
#   return(c(eq1, eq2))
# }
# 
# bet_val <- nleqslv(c(1,1), fn)
# bet_mle <- bet_val$x
# alpp_mle <-  n / sum(exp((Data-Tau)*bet_mle[2] + Tau*bet_mle[1]) -1)
# ## Let us define the CDF of GD under high stress
# cdf_gd2 <- function(x){
#   cdf <- 1 - exp(-alpp_mle*(exp((x-Tau)*bet_mle[2] + bet_mle[1]*Tau)-1))
#   return(cdf)
# }
# 
# 
# # Apply KS test
# ks2 = ks.test(Data,cdf_gd2)
# ks2








par(mfrow = c(1,3))
# Empirical CDF
emp_CDF<- numeric(0)
for(i in 1:n){
  count=0
  for(j in 1:n){
    if(Data[j] <= Data[i]){
      count= count+ 1
      emp_CDF[i] <- count/n
    }
  }
}

# Theoritical CDF
theoretical_CDF <- cdf_gd1(Data);
# theoretical_CDF <- cdf_gd2(Data);
emp_CDF; 
theoretical_CDF; 

## Apply CDF CDF plot technique
plot(Data,emp_CDF, type="s", col="chartreuse3",ylab="F(x)",xlab= "Data", lwd = 2, main = "Under stress level 2.25 V")
lines(Data,theoretical_CDF, type="l", col="tomato",lwd = 2, lty = 1) 
legend("topleft", legend=c("Empirical_cdf", "Theoritical_cdf"),col=c("chartreuse3", "tomato"),lty=1:1, lwd = 2:2) 

##3) Apply QQ-Plot technique
sample_quantiles <- sort(Data)
theoretical_quantiles <- numeric(0)
proportion <- (1:n)/n
for(i in 1:length(proportion)){
  theoretical_quantiles[i] <- 1/bet_mle*(log(1 - 1/alpp_mle*log(1-proportion[i])))
  # theoretical_quantiles[i] <- 1/bet_mle[2] * log(1- log(1-proportion[i])/alpp_mle) - Tau*(bet_mle[1]/bet_mle[2]) + Tau 
}
plot(sample_quantiles,theoretical_quantiles,col="tomato", main = "Under stress level 2.25 V",pch = 19)
abline(0,1,col = "chartreuse3")                  # Add a diagonal line in q-q plot

##4) Apply PP-Plot technique
probDist <- cdf_gd1(Data)       #get probability distribution for data points
# probDist <- cdf_gd2(Data)       #get probability distribution for data points
plot(ppoints(length(Data)), sort(probDist), pch = 19, xlab = "Observed Probability", ylab = "Expected Probability",col="tomato",main = "Under stress level 2.25 V")
abline(0,1, col = "chartreuse3")

