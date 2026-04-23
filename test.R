# This script is AI generated using GPT-5.2-Codex with the following prompts:
# - "I want to write a JAGS model to implement N-mixture model for count data. Write me the JAGS model and related R code to run it, include a small simulated data example too."
# - "Have you considered that JAGS expects the precision instead of variance for the Normal priors? Can you make sure you get that right?"

# ----------------------

## N-mixture model example (Poisson abundance + Binomial detection)
## Simulates repeated-count data and fits the model in JAGS via rjags.

set.seed(123)

if (!requireNamespace("rjags", quietly = TRUE)) {
  stop("Package 'rjags' is required. Install with install.packages('rjags').")
}
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Package 'coda' is required. Install with install.packages('coda').")
}

library(rjags)
library(coda)

## -----------------------------
## 1) Simulate repeated-count data
## -----------------------------
nsite <- 100
nvisit <- 3

## Site-level covariate for abundance
x <- rnorm(nsite, mean = 0, sd = 1)

## True parameters used for simulation
beta0_true <- 0.4
beta1_true <- 0.8
alpha0_true <- -0.2

lambda_true <- exp(beta0_true + beta1_true * x)
N_true <- rpois(nsite, lambda_true)
p_true <- plogis(alpha0_true)

## Observation matrix y[i, j]
y <- matrix(NA_integer_, nrow = nsite, ncol = nvisit)
for (i in seq_len(nsite)) {
  for (j in seq_len(nvisit)) {
    y[i, j] <- rbinom(1, size = N_true[i], prob = p_true)
  }
}

## -----------------------------
## 2) JAGS model definition
## -----------------------------
nmix_model <- "
model {
  ## Priors
  beta0 ~ dnorm(0, tau_beta0)
  beta1 ~ dnorm(0, tau_beta1)
  alpha0 ~ dnorm(0, tau_alpha0)

  ## Likelihood
  for (i in 1:nsite) {
    log(lambda[i]) <- beta0 + beta1 * x[i]
    N[i] ~ dpois(lambda[i])

    logit(p[i]) <- alpha0

    for (j in 1:nvisit) {
      y[i, j] ~ dbin(p[i], N[i])
    }
  }

  ## Derived quantities
  Ntotal <- sum(N[])
  lambda_mean <- mean(lambda[])
}
"

## Write model to a temporary file for JAGS
model_file <- tempfile(fileext = ".jags")
writeLines(nmix_model, con = model_file)

## -----------------------------
## 3) Prepare data and run MCMC
## -----------------------------
jags_data <- list(
  y = y,
  x = x,
  nsite = nsite,
  nvisit = nvisit,
  ## JAGS dnorm() uses precision tau = 1 / sd^2
  tau_beta0 = 1 / 10^2,
  tau_beta1 = 1 / 10^2,
  tau_alpha0 = 1 / 10^2
)

## Reasonable latent state initials: >= max observed count at each site
inits_fn <- function() {
  list(
    beta0 = rnorm(1, 0, 1),
    beta1 = rnorm(1, 0, 1),
    alpha0 = rnorm(1, 0, 1),
    N = apply(y, 1, max) + rpois(nsite, 1)
  )
}

params <- c("beta0", "beta1", "alpha0", "Ntotal", "lambda_mean")

nc <- 3
ni <- 8000
nb <- 3000
nt <- 5

jm <- jags.model(
  file = model_file,
  data = jags_data,
  inits = inits_fn,
  n.chains = nc,
  n.adapt = 1000
)

update(jm, n.iter = nb)

post <- coda.samples(
  model = jm,
  variable.names = params,
  n.iter = ni,
  thin = nt
)

## -----------------------------
## 4) Summarize posterior results
## -----------------------------
print(summary(post))
cat("\nGelman-Rubin diagnostics (close to 1 is good):\n")
print(gelman.diag(post, autoburnin = FALSE))

## Compare posterior means vs truth used for simulation
sm <- summary(post)$statistics
cat("\nPosterior mean vs truth:\n")
cat(sprintf(
  "beta0:  est = %0.3f, true = %0.3f\n",
  sm["beta0", "Mean"],
  beta0_true
))
cat(sprintf(
  "beta1:  est = %0.3f, true = %0.3f\n",
  sm["beta1", "Mean"],
  beta1_true
))
cat(sprintf(
  "alpha0: est = %0.3f, true = %0.3f\n",
  sm["alpha0", "Mean"],
  alpha0_true
))
cat(sprintf(
  "p:      est = %0.3f, true = %0.3f\n",
  plogis(sm["alpha0", "Mean"]),
  p_true
))
cat(sprintf(
  "Ntotal: est = %0.1f, true = %0.1f\n",
  sm["Ntotal", "Mean"],
  sum(N_true)
))

## Optional traceplots
## plot(post)
