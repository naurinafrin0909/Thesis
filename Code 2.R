library(data.table)
library(stringr)
library(dplyr)
library(mgcv)
library(splines)

# ==============================================================================
# 1. Generalized Double GLM (DGLM) Function
# ==============================================================================
DGLM_function <- function(mean.formula, dispersion.pred, data, method = "gam") {
  epsilon <- 1e-07
  maxit   <- 100
  iter    <- 0
  dev     <- 0
  devdiff <- 1
  
  # Extract response vector dynamically from the formula
  y <- model.response(model.frame(mean.formula, data = data))
  
  data$W_t <- 1
  ini.model <- glm(formula = mean.formula, family = quasipoisson(link = "log"), data = data)
  
  while (devdiff > epsilon && iter < maxit) {
    # Update mean model with observation weights
    ED_glm <- update(ini.model, weights = W_t, data = data)
    mu <- ED_glm$fitted.values
    
    # Calculate Poisson unit deviance: d_i = 2 * [y * log(y / mu) - (y - mu)]
    d <- numeric(length(y))
    nonzero <- y > 0
    d[nonzero]  <- 2 * (y[nonzero] * log(y[nonzero] / mu[nonzero]) - (y[nonzero] - mu[nonzero]))
    d[!nonzero] <- 2 * mu[!nonzero]
    
    # Guard against non-positive deviances (which crash the Gamma family)
    d <- pmax(d, 1e-6)
    data$d <- d
    
    dispersion.formula <- reformulate(dispersion.pred, response = "d")
    
    # Fit deviance via Gamma model with log link
    if (method == "gam") {
      d_glm <- gam(dispersion.formula, family = Gamma(link = "log"), data = data, select = TRUE, method = "REML")
    } else if (method == "glm") {
      d_glm <- glm(dispersion.formula, family = Gamma(link = "log"), data = data)
    }
    
    # Var(Y_i) = phi_i * mu_i, where phi_i is fitted dispersion
    var_y <- d_glm$fitted.values * mu
    
    # For log link, g'(mu) = 1/mu
    d_gmu <- 1 / mu
    
    # Observation weight: W_i = 1 / [Var(Y_i) * (g'(mu))^2] = mu_i / phi_i
    data$W_t <- 1 / (var_y * (d_gmu)^2)
    
    dev_t <- sum(d)
    devdiff <- abs(dev_t - dev)
    dev <- dev_t
    iter <- iter + 1
  }
  
  output <- list(ED_glm = ED_glm, d_glm = d_glm, dev = dev, iter = iter)
  return(output)
}

# ==============================================================================
# 2. Data Loading & Preprocessing
# ==============================================================================
# Read features and labels from files (or URLs)
features <- read.csv("https://s3.amazonaws.com/drivendata/data/44/public/dengue_features_train.csv")
labels   <- read.csv("https://s3.amazonaws.com/drivendata/data/44/public/dengue_labels_train.csv")

# Merge features and response (total_cases)
dat <- merge(features, labels, by = c("city", "year", "weekofyear"))

# Filter to San Juan ('sj') to maintain an uninterrupted weekly time-series
dat <- dat %>%
  filter(city == "sj") %>%
  arrange(as.Date(week_start_date))

# Date, calendar, and sequence transformations
dat$week_start_date <- as.Date(dat$week_start_date)
dat$month           <- as.factor(format(dat$week_start_date, "%m"))
dat$year            <- as.factor(dat$year)
dat$week_num        <- as.numeric(dat$weekofyear)
dat$time_idx        <- seq_len(nrow(dat))

# Impute median values for missing weather sensor records
continuous_vars <- c(
  "station_avg_temp_c",
  "reanalysis_relative_humidity_percent",
  "station_precip_mm"
)

for (col in continuous_vars) {
  dat[[col]][is.na(dat[[col]])] <- median(dat[[col]], na.rm = TRUE)
}

# ==============================================================================
# 3. Model Specification & Execution
# ==============================================================================
# Covariates for the Mean Model:
# - ns(station_avg_temp_c, 3): Flexible air temperature curve
# - ns(reanalysis_relative_humidity_percent, 3): Humidity effect
# - ns(station_precip_mm, 3): Precipitation effect
# - ns(time_idx, df = 18 * 4): Multi-year trend spline (~4 df per year)
mean_variables <- c(
  "ns(station_avg_temp_c, 3)",
  "ns(reanalysis_relative_humidity_percent, 3)",
  "ns(station_precip_mm, 3)",
  "ns(time_idx, 18 * 4)"
)

model.mean <- reformulate(mean_variables, response = "total_cases")

# Fit the DGLM
# Note: Using month (or cyclic spline on week_num) for dispersion avoids factor rank-deficiency
fit <- DGLM_function(
  mean.formula    = model.mean,
  dispersion.pred = "month",
  data            = dat,
  method          = "gam"
)

# ==============================================================================
# 4. Summaries
# ==============================================================================
summary(fit$ED_glm)  # Quasi-poisson mean regression
summary(fit$d_glm)   # Gamma dispersion regression
# 1. Fit the standard Quasi-Poisson (constant dispersion)
fit_quasipoisson <- glm(
  formula = model.mean, 
  family  = quasipoisson(link = "log"), 
  data    = dat
)

summary(fit$d_glm)
# ==============================================================================
# 1. Fit the Comparison Models
# ==============================================================================

# Baseline standard Poisson (has an exact AIC, but assumes phi = 1)
fit_pois <- glm(formula = model.mean, family = poisson(link = "log"), data = dat)

# Standard Quasi-Poisson (constant overdispersion phi)
fit_qp   <- glm(formula = model.mean, family = quasipoisson(link = "log"), data = dat)

# DGLM (already fitted; varying overdispersion phi_i by month)
# fit$ED_glm (mean model), fit$d_glm (dispersion model)

# ==============================================================================
# 2. Extract Key Quantities
# ==============================================================================
y <- dat$total_cases
n <- length(y)

# Degrees of freedom (number of estimated parameters)
# Mean model parameters (k_mean)
k_mean <- length(coef(fit_qp))

# Dispersion model parameters:
# For Quasi-Poisson: 1 constant phi parameter
# For DGLM: coefficients of the dispersion model
k_disp_qp   <- 1
k_disp_dglm <- length(coef(fit$d_glm))

k_qp   <- k_mean + k_disp_qp
k_dglm <- k_mean + k_disp_dglm

# Fitted means
mu_qp   <- fit_qp$fitted.values
mu_dglm <- fit$ED_glm$fitted.values

# Fitted dispersions (phi)
phi_qp   <- summary(fit_qp)$dispersion
phi_dglm <- fit$d_glm$fitted.values

# ==============================================================================
# 3. Method A: QAIC (Quasi-Akaike Information Criterion)
# Formula: QAIC = -2 * (logLik_poisson / c_hat) + 2 * K
# ==============================================================================
# Log-likelihood under standard Poisson
loglik_poisson_qp   <- sum(dpois(y, lambda = mu_qp, log = TRUE))
loglik_poisson_dglm <- sum(dpois(y, lambda = mu_dglm, log = TRUE))

qaic_qp   <- -2 * (loglik_poisson_qp / phi_qp) + 2 * k_qp
qaic_dglm <- -2 * sum(dpois(y, lambda = mu_dglm, log = TRUE) / phi_dglm) + 2 * k_dglm

# ==============================================================================
# 4. Method B: Exact DGLM Extended Quasi-Likelihood AIC (Smyth, 1989)
# Formula: -2 * logL(mu, phi) + 2 * (k_mean + k_disp)
# ==============================================================================
# Unit deviances
dev_qp   <- 2 * ifelse(y > 0, y * log(y / mu_qp) - (y - mu_qp), mu_qp)
dev_dglm <- fit$ED_glm$data$d

# Joint log-likelihood approximation (Saddlepoint / Extended Quasi-Likelihood)
eql_qp   <- -0.5 * sum(log(2 * pi * phi_qp * y + 1e-6) + dev_qp / phi_qp)
eql_dglm <- -0.5 * sum(log(2 * pi * phi_dglm * y + 1e-6) + dev_dglm / phi_dglm)

aic_eql_qp   <- -2 * eql_qp + 2 * k_qp
aic_eql_dglm <- -2 * eql_dglm + 2 * k_dglm

# ==============================================================================
# 5. Output Summary Table
# ==============================================================================
comparison_aic <- data.frame(
  Model           = c("Standard Poisson (phi = 1)", "Quasi-Poisson (constant phi)", "DGLM (varying phi)"),
  Parameters_K    = c(k_mean, k_qp, k_dglm),
  QAIC            = c(AIC(fit_pois), qaic_qp, qaic_dglm),
  Extended_QL_AIC = c(NA, aic_eql_qp, aic_eql_dglm)
)

print(comparison_aic)
head(features)
head(labels)