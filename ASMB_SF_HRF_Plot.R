par(mfrow = c(1, 2))
S0 = 2
S1 = 2.25
S2 = 2.44
ml_alp = mle[1]
ml_beta1 = mle[2]
ml_beta2 = mle[3]
bs_alp = bayes_est[1]
bs_beta1 =  bayes_est[2]
bs_beta2 =  bayes_est[3]

ml_b = (log(ml_beta1) - log(ml_beta2)) / (log(S1) - log(S2))
ml_a = 1 / 2 * (log(ml_beta1) + log(ml_beta2) - ml_b * (log(S1) + log(S2)))
ml_beta0 = exp(ml_a + ml_b * log(S0))

bs_b = (log(bs_beta1) - log(bs_beta2)) / (log(S1) - log(S2))
bs_a = 1 / 2 * (log(bs_beta1) + log(bs_beta2) - bs_b * (log(S1) + log(S2)))
bs_beta0 = exp(bs_a + bs_b * log(S0))


##....MTTF and Median....
mttf <- function(x) {
  exp(ml_alp) / ml_beta0 * exp(-x) / x
}
MTTF <- integrate(mttf, lower = ml_alp, upper = Inf)$value

qt <- function(p) {
  1 / ml_beta0 * log(1 - 1 / ml_alp * log(1 - p))
}

quantle <- qt(0.5)


ml_R0 <- function(t1) {
  r <-  exp(-ml_alp * (exp(ml_beta0 * t1) - 1))
  return(r)
}
ml_h0 <- function(t2) {
  h <- ml_alp * ml_beta0 * exp(ml_beta0 * t2)
  return(h)
}

rsol = integrate(ml_R0, lower = 0, upper = Inf)



bs_R0 <- function(t1) {
  r <-  exp(-bs_alp * (exp(bs_beta0 * t1) - 1))
  return(r)
}
bs_h0 <- function(t2) {
  h <- bs_alp * bs_beta0 * exp(bs_beta0 * t2)
  return(h)
}


# Define new attractive colors
color_MLE <- "#FF5733"  # Vibrant orange for MLE
color_BE <- "#1F77B4"   # Bold blue for BE

# Plot 1: R0 curves with new colors
curve(
  ml_R0,
  from = 0,
  to = 600,
  col = color_MLE,
  lwd = 2,
  lty = 1,
  xlab = expression(t[0]),
  ylab = expression(R[0](t[0]))
)

curve(
  bs_R0,
  from = 0,
  to = 600,
  col = color_BE,
  lwd = 2,
  lty = 1,
  add = TRUE
)
legend(
  "topright",
  legend = c("MLE", "BE"),
  col = c(color_MLE, color_BE),
  lwd = c(2, 2),
  lty = c(1, 1)
)

# Plot 2: h0 curves with new colors
curve(
  ml_h0,
  from = 0,
  to = 1000,
  col = color_MLE,
  lwd = 2,
  lty = 1,
  xlab = expression(t[0]),
  ylab = expression(h[0](t[0]))
)

curve(
  bs_h0,
  from = 0,
  to = 1000,
  col = color_BE,
  lwd = 2,
  lty = 1,
  add = TRUE
)
legend(
  "topleft",
  legend = c("MLE", "BE"),
  col = c(color_MLE, color_BE),
  lwd = c(2, 2)
)