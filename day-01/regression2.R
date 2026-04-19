## Multivariate analysis of spatial ecological data using R (MASE01)
## Northwest Atlantic Fisheries Centre
## St. John's, 19 June 2017 -- 23 June 2017
## Subhash Lele & Peter Solymos
## Module 2: generalized linear models

## set root directory
ROOT <- "~/Dropbox/courses/st-johns-2017"

## load packages and other required code
library(visreg)
library(rpart)
library(ResourceSelection)
library(mgcv)
library(intrval)
library(MASS)
library(lattice)
library(mefa4)
library(pROC)
source(file.path(ROOT, "R", "diagnostics-functions.R"))

## read in the data file
x <- read.csv(file.path(ROOT, "data", "ABbirds", "YellowrumpedWarbler.csv"))

## groups of variables for easier column subsetting
## land cover composition
cn_lc <- c(
  "pforest",
  "pconif",
  "pwet",
  "poldfor",
  "Cult",
  "UrbInd",
  "SoftLin",
  "HardLin",
  "CC",
  "Alien",
  "Succ",
  "THF"
)
## bioclimatic variables
cn_cl <- c("AHM", "PET", "FFP", "MAP", "MAT", "MCMT", "MWMT")
## response y and the two combined
cn <- c("y", cn_cl, cn_lc)

## look at the counts: peak at 0
table(x$y)
plot(table(x$y))

d <- readRDS("data/ab_birds_site_year.rds")
spp <- "yellowrumped_warbler"

x <- data.frame(y = d$counts[, spp], d$vars)

x$log_y <- log(x$y + 0.5)

## Data summaries
## prints the first few rows of tables
head(x) # try also tail(x)
## structure of objects: numeric/factor/etc types
str(x)
## summary stats and number of NAs
summary(x)
## have a look at correlations of the variables
round(cor(x[, 9:15]), 2)

## base graphics: univariate exploration
plot(table(x$y)) # barplot
hist(x$log_y) # histogram
## summarize individual variable (range, mean, median, quartiles)
summary(x$log_y)

## base graphics: bivariate exploration
plot(log_y ~ pforest, x) # scatterplot
boxplot(y ~ nrname, x, horizontal = TRUE, las = 1) # boxplot

## base graphics: multivariate exploration
with(
  x,
  plot(
    latitude,
    longitude, # plot the species detections in space
    pch = ifelse(x$y > 0, 19, 3),
    asp = 1
  )
)
## bivariate plot with bubbles
with(x, plot(alien, succ, cex = 0.5 + pforest))
with(x, plot(pforest, pconif, cex = 0.5 + x$y / 20))
## scatterplot matrix
pairs(x[, 9:15], pch = ".")
## look at y vs. multiple x variables
plot(log_y ~ ., data = x[, cl])

## lattice graphics: facetting
library(lattice)
histogram(~ THF | NRNAME, x)
xyplot(Succ ~ Alien | NRNAME, x)
xyplot(y ~ pforest | NRNAME, x)
xyplot(y ~ pconif | NRNAME, x)
splom(x[, 9:15])

## fit a linear model, intercept only, no covariates

## lm(log_y ~ 1)
fit0 <- lm(log_y ~ 1, data = x)
summary(fit0) # summary
summary(fit0)$sigma # residual SD
coef(fit0) # coefficient
confint(fit0) # confidence interval

## prediction: 95% CI and PI for *new* observations
## mu=beta0
xnew <- x
xnew <- xnew[order(xnew$log_y), ] # sorted by y
## confidence interval
pr_ci <- predict(fit0, newdata = xnew, interval = "confidence", level = 0.95)
## prediction interval
pr_pi <- predict(fit0, newdata = xnew, interval = "prediction", level = 0.95)
head(pr_ci) # show the head of the tables
head(pr_pi)

## plot log_y, CI, PI -- clearly need more covariates
plot(xnew$log_y, type = "l", ylim = range(pr_ci, pr_pi))
abline(h = pr_ci[1, "fit"], col = 2, lty = 1)
abline(h = pr_ci[1, c("lwr", "upr")], col = 2, lty = 2)
abline(h = pr_pi[1, c("lwr", "upr")], col = 4, lty = 2)
legend(
  "topleft",
  col = c(1, 2, 2, 4),
  lty = c(1, 1, 2, 2),
  lwd = 1,
  legend = c("log_y", "fit", "CI", "PI"),
  bty = "n"
)

cn_lc <- c(
  "pforest",
  "pconif",
  "pwet",
  "poldfor",
  "cult",
  "urbind",
  "softlin",
  "hardlin",
  "cc",
  "alien",
  "succ",
  "thf"
)
## bioclimatic variables
cn_cl <- c("ahm", "pet", "ffp", "map", "mat", "mcmt", "mwmt")
## response log_y and the two combined
cn <- c("log_y", cn_cl, cn_lc)

ct <- partykit::ctree(log_y ~ ., data = x[, cn])
plot(ct)


## semiparametric exploration with GAM
fit_gam <- mgcv::gam(
  log_y ~ s(pforest) + s(pconif) + s(poldfor) + s(alien) + s(ahm) + s(map),
  data = x,
  family = gaussian
)
summary(fit_gam)

## mu=beta0 + combination of splines
pred_gam <- fitted(fit_gam) # there is also predict function for GAM

## fitted splines with residuals, all on 1 page
plot(fit_gam, residuals = TRUE, pages = 1, scheme = 1)
## bivariate fitted surfaces with detections and nondetection on top
## try plot.type="persp" as well
mgcv::vis.gam(fit_gam, c("pforest", "ahm"), plot.type = "contour") # "contour" or "persp"
#with(x[x$log_y < mean(x$log_y),], points(pforest, ahm, cex=0.6, col="turquoise"))
#with(x[x$log_y >= mean(x$log_y),], points(pforest, ahm, cex=0.6, col="blue"))

## lm with covariates

## lm(log_y ~ x)

## note: here we use the untransformed counts as response!
fit1y <- lm(y ~ pconif, data = x)
plot(fit1y, 1) # not constant variance of residuals

## log_y as response, compare residuals to the previous one
fit1 <- lm(log_y ~ pconif, data = x)
plot(fit1, 1) # more homoskedastic now, probably nonlinear
summary(fit1)
anova(fit1) # F-test
summary(fit1)$sigma # residual standard deviation
summary(fit1)$r.squared # R^2
summary(fit1)$adj.r.squared # adjusted R^2
coef(fit1) # estimated coefficients
confint(fit1) # confidence intervals (use levels if other than 95% is needed)
summary(fitted(fit1)) # summary of the fitted values
summary(x$log_y) # data mean should match with fitted mean

## prediction
## mu[i]=beta0 + beta1*x[i]
xnew <- x
xnew <- xnew[order(xnew$pconif), ] # order by pconif
pr_ci <- predict(fit1, newdata = xnew, interval = "confidence", level = 0.95)
pr_pi <- predict(fit1, newdata = xnew, interval = "prediction", level = 0.95)

## plot the pbserved, fitted values and the intervals
with(xnew, plot(pconif, pr_ci[, "fit"], type = "l", ylim = range(pr_ci, pr_pi)))
polygon(
  c(xnew$pconif, rev(xnew$pconif)),
  c(pr_pi[, "lwr"], rev(pr_pi[, "upr"])),
  col = rgb(0, 0, 1, 0.1),
  border = NA
)
with(xnew, lines(pconif, pr_pi[, "upr"], col = 4, lty = 2))
with(xnew, lines(pconif, pr_pi[, "lwr"], col = 4, lty = 2))
polygon(
  c(xnew$pconif, rev(xnew$pconif)),
  c(pr_ci[, "lwr"], rev(pr_ci[, "upr"])),
  col = rgb(1, 0, 0, 0.3),
  border = NA
)
with(xnew, lines(pconif, pr_ci[, "upr"], col = 2, lty = 2))
with(xnew, lines(pconif, pr_ci[, "lwr"], col = 2, lty = 2))
with(xnew, points(pconif, log_y, cex = 0.2, col = "turquoise"))

## almost nominal coverage
library(intrval)
table(xnew$log_y %[]% pr_pi[, c("lwr", "upr")]) / nrow(xnew)

## diagnostics
plot(fit1, which = 1:6)
x[c(132, 192, 441), ] # high leverage points

## making a QQ plot
res <- residuals(fit1)
qqnorm(res)
qqline(res)


## lm(log_y ~ x + z)

## use pconif and heat moisture index (AHM)
fit2 <- lm(log_y ~ pconif + ahm, data = x)
summary(fit2)
#anova(fit2)
#summary(fit2)$sigma
#coef(fit2)
#confint(fit2)

## mu[i]=beta0 + beta1*x1[i] + beta2*x2[i]
pr_ci <- predict(fit2, newdata = xnew, interval = "confidence", level = 0.95)
pr_pi <- predict(fit2, newdata = xnew, interval = "prediction", level = 0.95)
## coverage
table(xnew$log_y %[]% pr_pi[, c("lwr", "upr")]) / nrow(xnew) # good

## 2D fitted surface using expand.grid and predict
xnew <- expand.grid(
  pconif = seq(0, 1, length.out = 50),
  ahm = seq(min(x$ahm), max(x$ahm), length.out = 50)
)
str(xnew)
xnew$log_y_pred <- predict(fit2, newdata = xnew) # only fit, no CI/PI
pr_list <- list(
  x = seq(0, 1, length.out = 50),
  y = seq(min(x$ahm), max(x$ahm), length.out = 50),
  z = matrix(xnew$log_y_pred, 50, 50)
)
image(pr_list, xlab = "pconif", ylab = "AHM")
contour(pr_list, add = TRUE)
#with(x[x$log_y < mean(x$log_y),], points(pforest, AHM, cex=0.6, col="turquoise"))
#with(x[x$log_y >= mean(x$log_y),], points(pforest, AHM, cex=0.6, col="blue"))
## PI now is the interval between 2 surfaces

## let us now overfit
fit3 <- lm(
  log_y ~ (pforest + pconif + I(pconif^2) + ahm)^3,
  data = x,
  na.action = na.fail
)
summary(fit3)

## backwards selection: start from fit3 and drop terms one at a time
fit4 <- step(fit3, direction = "backward")
## forward selection: start with the null model and add terms one at a time
fit5 <- step(lm(log_y ~ 1, x), scope = formula(fit3), direction = "forward")
## different 'best' models depending on direction
formula(fit4) # backwards
formula(fit5) # forwards
## does it matter for prediction?
plot(fitted(fit4), fitted(fit5))
abline(0, 1, lty = 2, col = 2)

aic <- AIC(fit0, fit1, fit2, fit3, fit4, fit5)
aic$dAIC <- aic$AIC - min(aic$AIC)
aic

## put all out models in a list
mList <- list(
  fit0 = fit0,
  fit1 = fit1,
  fit2 = fit2,
  fit3 = fit3,
  fit4 = fit4,
  fit5 = fit5
)


## select the best model based on AIC
fitBest <- mList[[which.min(aic$AIC)]]

## run diagnostics
plot(fitBest, which = 1:6)

## conditional plots: effect of x on E[Y] holding all
visreg::visreg(fitBest) # linear conditional responses
visreg::visreg2d(fitBest, "pconif", "ahm", plot.type = "image") # try persp or rgl
visreg::visreg2d(fitBest, "pconif", "pforest", plot.type = "image") # try persp or rgl
visreg::visreg2d(fitBest, "pforest", "ahm", plot.type = "image") # try persp or rgl


predict_sim <-
  function(
    object,
    newdata = NULL,
    interval = c("none", "confidence", "prediction"),
    type = c("asymp", "pboot", "npboot"),
    level = 0.95,
    B = 99,
    ...
  ) {
    interval <- match.arg(interval)
    type <- match.arg(type)
    if (is.null(newdata)) {
      x <- model.frame(object)
      X <- model.matrix(object)
    } else {
      x <- model.frame(delete.response(terms(object)), newdata)
      X <- model.matrix(attr(x, "terms"), x)
    }
    n <- nrow(x)
    fun <- switch(
      family(object)$family,
      "gaussian" = function(x) rnorm(length(x), x, summary(object)$sigma),
      "poisson" = function(x) rpois(length(x), x),
      "binomial" = function(x) rbinom(length(x), 1, x)
    )
    if (interval == "none") {
      return(predict(object, newdata))
    }
    if (B < 2) {
      stop("Are you kidding? B must be > 1")
    }
    if (type == "asymp") {
      cm <- rbind(
        coef(object),
        MASS::mvrnorm(B, coef(object), vcov(object))
      )
      fm <- apply(cm, 1, function(z) X %*% z)
    }
    if (type == "boot") {
      cm <- matrix(0, B + 1, length(coef(object)))
      cm[1, ] <- coef(object)
      xx <- model.frame(object)
      for (i in 2:B) {
        j <- sample.int(n, n, replace = TRUE)
        cm[i, ] <- coef(update(object, data = xx[j, ]))
      }
      fm <- apply(cm, 1, function(z) X %*% z)
    }
    if (type == "npboot") {
      cm <- matrix(0, B + 1, length(coef(object)))
      cm[1, ] <- coef(object)
      xx <- model.frame(object)
      j <- attr(attr(xx, "terms"), "response")
      f <- fitted(object)
      for (i in 2:B) {
        xx[, j] <- fun(f)
        cm[i, ] <- coef(update(object, data = xx))
      }
      fm <- apply(cm, 1, function(z) X %*% z)
    }
    fm <- family(fit0)$linkinv(fm)
    if (interval == "prediction") {
      y <- matrix(fun(fm), n, B + 1)
    } else {
      y <- fm
    }
    rownames(y) <- rownames(x)
    p <- c(0.5, (1 - level) / 2, 1 - (1 - level) / 2)
    stat_fun <- function(x) {
      c(mean(x), sd(x), quantile(x, p))
    }
    out <- cbind(fm[, 1], t(apply(y, 1, stat_fun)))
    colnames(out) <- c("fit", "mean", "se", "median", "lwr", "upr")
    out[, c("fit", "lwr", "upr", "mean", "median", "se")]
  }

.r2_fun <-
  function(
    observed,
    fitted,
    distr = c("binomial", "poisson"),
    size = 1,
    null = NULL,
    p = 0
  ) {
    distr <- match.arg(distr)
    if (distr == "poisson") {
      if (is.null(null)) {
        null <- mean(observed)
      }
      ll0 <- sum(dpois(observed, null, log = TRUE))
      lls <- sum(dpois(observed, observed, log = TRUE))
      llf <- sum(dpois(observed, fitted, log = TRUE))
    } else {
      if (is.null(null)) {
        null <- mean(observed / size)
      }
      ll0 <- sum(dbinom(observed, size, null, log = TRUE))
      lls <- sum(dbinom(observed, size, observed / size, log = TRUE))
      llf <- sum(dbinom(observed, size, fitted, log = TRUE))
    }
    n <- length(observed)
    R2 <- 1 - (lls - llf) / (lls - ll0)
    R2adj <- 1 - (1 - R2) * ((n - 1) / (n - (p + 1)))
    D0 <- -2 * (ll0 - lls)
    DR <- -2 * (llf - lls)
    p_value <- 1 - pchisq(D0 - DR, p)
    c(
      R2 = R2,
      R2adj = R2adj,
      Deviance = D0 - DR,
      Dev0 = D0,
      DevR = DR,
      df = p,
      p_value = p_value
    )
  }

R2dev <-
  function(object, ...) {
    y <- model.response(model.frame(object), "numeric")
    f <- fitted(object)
    .r2_fun(
      y,
      fitted(object),
      distr = family(object)$family,
      size = 1,
      null = NULL,
      p = length(coef(object)) - 1
    )
  }

## ----------- POISSON --------------

## glm(y ~ 1)
fit0 <- glm(y ~ 1, data = x, family = poisson)
summary(fit0) # summary
coef(fit0) # coefficient
confint(fit0) # confidence interval

## prediction: 95% CI
## lambda=exp(beta0)
xnew <- x
xnew <- xnew[order(xnew$y), ]
pr_ci1 <- predict(fit0, newdata = xnew, type = "link", se.fit = TRUE)
pr_ci2 <- predict(fit0, newdata = xnew, type = "response", se.fit = TRUE)
str(pr_ci1) # beta0
str(pr_ci2) # exp(beta0)

plot(xnew$y, type = "l", ylim = range(xnew$y, pr_ci2$fit))
abline(h = pr_ci2$fit[1], col = 2, lty = 1)

## how do we get PI? (*new* observation)
pr_ci <- predict_sim(
  fit0,
  newdata = xnew,
  interval = "confidence",
  level = 0.95,
  B = 999
)
pr_pi <- predict_sim(
  fit0,
  newdata = xnew,
  interval = "prediction",
  level = 0.95,
  B = 999
)
head(pr_ci)
head(pr_ci2$fit)
head(pr_ci2$se.fit)
head(pr_pi)

## plot y, CI, PI -- clearly need more covariates
plot(xnew$y, type = "l", ylim = range(x$y, pr_ci, pr_pi))
abline(h = pr_ci[1, "fit"], col = 2, lty = 1)
abline(h = pr_ci[1, c("lwr", "upr")], col = 2, lty = 2)
abline(h = pr_pi[1, c("lwr", "upr")], col = 4, lty = 2)
legend(
  "topleft",
  col = c(1, 2, 2, 4),
  lty = c(1, 1, 2, 2),
  lwd = 1,
  legend = c("y", "fit", "CI", "PI"),
  bty = "n"
)

## much smaller than nominal: probably violating distr assumptions
## zero- and over-dispersion
table(xnew$y %[]% pr_pi[1, c("lwr", "upr")]) / nrow(xnew)

## diagnostics
plot(fit0)

## pairwise plots
plot(y ~ ., data = x[, cn])
## non-parametric exploration based on loess fit
## to the scatter
siplot(y ~ ., x[, cn])

## nonparametric exploration with CART
fit_rpart <- rpart(y ~ ., x[, cn], method = "poisson")
fit_rpart # look at most important variables

## display the tree
op <- par(xpd = TRUE) # all plotting clipped
plot(fit_rpart, compress = TRUE) # tree
text(fit_rpart, use.n = TRUE, cex = 0.5) # text
par(op)

## fitted/predicted values
## lambda=mean of obs in bucket according to splits
pred_rpart <- predict(fit_rpart, newdata = x)

## semiparametric exploration with GAM
fit_gam <- gam(
  y ~ s(pconif) +
    s(poldfor) +
    s(MCMT) +
    s(Alien) +
    s(AHM) +
    s(MWMT) +
    s(MAP) +
    s(FFP),
  data = x,
  family = poisson
)
summary(fit_gam)

## lambda=exp(beta0 + combination of splines)
pred_gam <- fitted(fit_gam)

## fitted splines with residuals, all on 1 page
plot(fit_gam, scheme = 1)
## bivariate fitted surfaces with detections and nondetection on top
## try plot.type="persp" as well
vis.gam(fit_gam, c("Alien", "MWMT"), plot.type = "contour") # "contour" or "persp"
#with(x[x$y == 0,], points(Alien, MWMT, cex=0.6, col="turquoise"))
#with(x[x$y > 0,], points(Alien, MWMT, cex=0.6, col="blue"))

## glm with covariates

## glm(y ~ x)

fit1 <- glm(y ~ pconif, data = x, family = poisson)
summary(fit1)
coef(fit1)
confint(fit1)
anova(fit1, test = "Chi")

## lambda=exp(beta0 + beta1*x[i])
xnew <- x
xnew <- xnew[order(xnew$pconif), ]
pr_ci <- predict(fit1, newdata = xnew, type = "response")
with(xnew, plot(pconif, pr_ci, type = "l"))

plot(fit1) # diagnostics
meptest(fit1) # check the form

fit2 <- glm(y ~ pconif + I(pconif^2), data = x, family = poisson)
summary(fit2)
coef(fit2)
confint(fit2)
anova(fit2, test = "Chi")

## lambda=exp(beta0 + beta1*x[i] + beta2*x[i]^2)
pr_ci <- predict(fit2, newdata = xnew, type = "response")
with(xnew, plot(pconif, pr_ci, type = "l"))

plot(fit2) # diagnostocs
meptest(fit2) # check the form

## added variable plot
avp(fit2, x[, cn_cl])
avp(fit2, x[, cn_lc])

## ultra alles inklusive
fit3 <- glm(
  y ~ pconif + I(pconif^2) + poldfor + Alien + AHM + MWMT + MCMT + MAP + FFP,
  x,
  family = poisson
)

visreg(fit3, scale = "response")
visreg2d(fit3, "poldfor", "pconif", scale = "linear", plot.type = "image") # try persp and rgl

mep(fit3)
meptest(fit3)

fit4 <- glm(
  y ~ pconif +
    I(pconif^2) +
    poldfor +
    I(poldfor^2) +
    Alien +
    AHM +
    MWMT +
    MCMT +
    MAP +
    FFP,
  x,
  family = poisson
)
meptest(fit4)

fit5 <- step(fit4)
formula(fit4)
formula(fit5)

aic <- AIC(fit0, fit1, fit2, fit3, fit4, fit5)
aic$dAIC <- aic$AIC - min(aic$AIC)
aic

## H0: the null model fits the data at least as well as the alternative model
## table with comulative deviance explained
(at <- anova(fit0, fit1, fit2, fit3, fit4, fit5, test = "Chi"))

## pseudo R^2: based on full, saturated, and null model
R2dev(fit5)
(R2 <- data.frame(t(sapply(
  list(
    fit0 = fit0,
    fit1 = fit1,
    fit2 = fit2,
    fit3 = fit3,
    fit4 = fit4,
    fit5 = fit5
  ),
  R2dev
))))

## ------------ BINOMIAL ------------

x$y01 <- ifelse(x$y > 0, 1, 0)

## glm(y ~ 1)
fit0 <- glm(y01 ~ 1, data = x, family = binomial)
summary(fit0) # summary
coef(fit0) # coefficient
confint(fit0) # confidence interval

## prediction: 95% CI
## p=exp(beta0) / (1 + exp(beta0)) = 1/(1+exp(-exp(beta0)))
xnew <- x
xnew <- xnew[order(xnew$y), ]
pr_ci1 <- predict(fit0, newdata = xnew, type = "link", se.fit = TRUE)
pr_ci2 <- predict(fit0, newdata = xnew, type = "response", se.fit = TRUE)
str(pr_ci1) # beta0
str(pr_ci2) # plogis(beta0) = exp(beta0) / (1 + exp(beta0))

plot(xnew$y01, type = "l", ylim = range(0, 1))
abline(h = pr_ci2$fit[1], col = 2, lty = 1)

## how do we get PI? (*new* observation)
pr_ci <- predict_sim(
  fit0,
  newdata = xnew,
  interval = "confidence",
  level = 0.95,
  B = 999
)
pr_pi <- predict_sim(
  fit0,
  newdata = xnew,
  interval = "prediction",
  level = 0.95,
  B = 999
)
head(pr_ci)
head(pr_ci2$fit)
head(pr_ci2$se.fit)
head(pr_pi)

## plot y, CI, PI -- clearly need more covariates
plot(xnew$y01, type = "l", ylim = range(0, 1))
abline(h = pr_ci[1, "fit"], col = 2, lty = 1)
abline(h = pr_ci[1, c("lwr", "upr")], col = 2, lty = 2)
abline(h = pr_pi[1, c("lwr", "upr")], col = 4, lty = 2)
legend(
  "bottomright",
  col = c(1, 2, 2, 4),
  lty = c(1, 1, 2, 2),
  lwd = 1,
  legend = c("y01", "fit", "CI", "PI"),
  bty = "n"
)

## coverage is 1, not really helpful
table(xnew$y01 %[]% pr_pi[1, c("lwr", "upr")]) / nrow(xnew)

## diagnostics
plot(fit0)

## non-parametric exploration based on loess
siplot(x$y01, x[, cn_lc])

## nonparametric exploration with CART
x$yfact <- as.factor(c("0", "1")[x$y01 + 1])
table(x$y01, x$yfact)
fit_rpart <- rpart(yfact ~ ., x[, c("yfact", cn[-1])], method = "class")
fit_rpart # look at most important variables

## fitted/predicted values
pred_rpart <- predict(fit_rpart, newdata = x)
head(pred_rpart)

## display the tree
op <- par(xpd = TRUE) # all plotting clipped to the figure region
plot(fit_rpart, compress = TRUE)
text(fit_rpart, use.n = FALSE)
par(op)

## semiparametric exploration with GAM
fit_gam <- gam(
  y01 ~ s(pconif) + s(poldfor) + s(AHM),
  data = x,
  family = binomial
)
summary(fit_gam)

## lambda=exp(beta0 + combination of splines)
pred_gam <- fitted(fit_gam)

## fitted splines with residuals, all on 1 page
plot(fit_gam)

## glm(y ~ x)

## the use of update
#fit1 <- glm(y01 ~ pconif, data=x, family=binomial)
fit1 <- update(fit0, . ~ . + pconif)
summary(fit1)
coef(fit1)
confint(fit1)
anova(fit1, test = "Chi")

## lambda=exp(beta0 + beta1*x[i])
xnew <- x
xnew <- xnew[order(xnew$pconif), ]
pr_ci <- predict(fit1, newdata = xnew, type = "response")
with(xnew, plot(pconif, pr_ci, type = "l"))

plot(fit1) # diagnostocs
mep(fit1)
meptest(fit1) # check the form

#fit2 <- glm(y01 ~ pconif + I(pconif^2), data=x,
#    family=binomial)
fit2 <- update(fit1, . ~ . + I(pconif^2))
summary(fit2)
coef(fit2)
confint(fit2)
anova(fit2, test = "Chi")

## lambda=exp(beta0 + beta1*x[i] + beta2*x[i]^2)
pr_ci <- predict(fit2, newdata = xnew, type = "response")
with(xnew, plot(pconif, pr_ci, type = "l"))

plot(fit2) # diagnostocs
meptest(fit2) # check the form

## added variable plot
avp(fit2, x[, cn_cl])
avp(fit2, x[, cn_lc])

## ultra alles inklusive
fit3 <- glm(
  y01 ~ pconif +
    I(pconif^2) +
    poldfor +
    Alien +
    AHM +
    MWMT +
    MCMT +
    MAP +
    FFP,
  x,
  family = binomial
)

visreg(fit3, scale = "response")

mep(fit3)
meptest(fit3)

fit4 <- step(fit3)
formula(fit3)
formula(fit4)

aic <- AIC(fit0, fit1, fit2, fit3, fit4)
aic$dAIC <- aic$AIC - min(aic$AIC)
aic

## H0: the null model fits the data at least as well as the alternative model
(at <- anova(fit0, fit1, fit2, fit3, fit4, test = "Chi"))

## pseudo R^2: based on full, saturated, and null model
(R2 <- data.frame(t(sapply(
  list(fit0 = fit0, fit1 = fit1, fit2 = fit2, fit3 = fit3, fit4 = fit4),
  R2dev
))))

## Hosmer-Lemeshow GoF test (p<alpha indicates poor fit)
hoslem.test(x$y01, fitted(fit4)) # p>0.05: no evidence for poor fit

## ROC/AUC
roc0 <- roc(x$y01, fitted(fit0))
roc1 <- roc(x$y01, fitted(fit1))
roc2 <- roc(x$y01, fitted(fit2))
roc3 <- roc(x$y01, fitted(fit3))
roc4 <- roc(x$y01, fitted(fit4))
roc4
plot(roc0, lty = 2)
lines(roc1, col = 4)
lines(roc4, col = 2)
auc(roc0)
auc(roc1)
auc(roc2)
auc(roc3)
auc(roc4)

## What should the link function be?
fit41 <- update(fit4, family = binomial("probit"))
fit42 <- update(fit4, family = binomial("cloglog"))
AIC(fit4, fit41, fit42)

xval <- seq(-10, 10, by = 0.1)
linkinv <- family(fit4)$linkinv
linkinv1 <- family(fit41)$linkinv
linkinv2 <- family(fit42)$linkinv
curve(linkinv, -6, 6)
curve(linkinv1, add = TRUE, col = 2)
curve(linkinv2, add = TRUE, col = 4)
curve(exp, add = TRUE, col = 4, lty = 2)
legend(
  "topleft",
  bty = "n",
  col = c(1, 2, 4, 4),
  lty = c(1, 1, 1, 2),
  legend = c("loglog", "probit", "cloglig", "log")
)

## offsets and transformation of response (Poisson case)

xC <- read.csv(file.path(ROOT, "data", "ABbirds", "YellowrumpedWarblerC.csv"))
xC$npt <- ifelse(xC$Method == "PT4", 4, 9) # sampling area
xC$off <- log(xC$npt / 9) # area as offset
xC$ytr <- round(xC$y / xC$npt) # transformed

boxplot(y ~ Method, xC)
boxplot(ytr ~ Method, xC)

fit_off <- list(
  fit11 = glm(y ~ pconif + I(pconif^2) + AHM, x, family = poisson), # equal effort (9 PT)
  fit12 = glm(y ~ pconif + I(pconif^2) + AHM, xC, family = poisson), # variable effort (4 or 9 PT)
  fit13 = glm(y ~ pconif + I(pconif^2) + AHM + Method, xC, family = poisson), # estimated
  fit14 = glm(
    y ~ pconif + I(pconif^2) + AHM + offset(off),
    xC,
    family = poisson
  ), # offset
  fit15 = glm(ytr ~ pconif + I(pconif^2) + AHM, xC, family = poisson)
) # transformed

lapply(fit_off, coef) # fit14 and fit15 are NOT equivalent

xnewC <- xC
xnewC$off <- 0
xnewC$Method <- "PT9"

ff <- sapply(fit_off, function(z) {
  predict(z, newdata = xnewC, type = "response")
})
boxplot(ff, range = 0, ylab = "Prediction")
matlines(1:ncol(ff), t(ff), col = rgb(0, 0, 0, 0.1), lty = 1)
boxplot(ff, range = 0, add = TRUE, col = rgb(1, 0.5, 0.5, 0.4))

## see also NB, ZIP, ZINB, we'll cover some related GLMMs as well
