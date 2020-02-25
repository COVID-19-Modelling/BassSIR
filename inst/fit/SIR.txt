model {
  for (i in 1:n_t){
    s[i] ~ dnorm(mu[i], tau)
    mu[i] <- 2 * (alpha1[i] * x1[i] + alpha2 * x2[i])

    alpha1[i] = (beta * (1 - A[i] / m) - nu) / 2
  }

  alpha2 = (- beta / m)
  tau ~ dgamma(1E-3, 1E-3)

  m ~ dunif(mx, mx * 10)
  beta ~ dunif(0, 1)
}
