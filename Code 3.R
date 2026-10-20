# ============================================================
# DENGUE DGLM ANALYSIS
# Weekly dengue case counts
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(mgcv)
library(splines)


# ============================================================
# 2. READ DATA
# ============================================================

features <- read.csv(
  "https://s3.amazonaws.com/drivendata/data/44/public/dengue_features_train.csv"
)

labels <- read.csv(
  "https://s3.amazonaws.com/drivendata/data/44/public/dengue_labels_train.csv"
)


# ============================================================
# 3. MERGE FEATURES AND LABELS
# ============================================================




# ============================================================
# 4. PREPARE VARIABLES
# ============================================================

dat <- dat %>%
  mutate(
    week_start_date = as.Date(week_start_date),
    
    # numerical measure of time
    time = as.numeric(
      week_start_date - min(week_start_date)
    ),
    
    # month as categorical variable
    month = factor(
      format(week_start_date, "%m")
    ),
    
    # city as categorical variable
    city = factor(city)
  )


# ============================================================
# 5. SELECT VARIABLES FOR ANALYSIS
# ============================================================

dat_model <- dat %>%
  select(
    total_cases,
    city,
    year,
    weekofyear,
    week_start_date,
    time,
    month,
    precipitation_amt_mm,
    reanalysis_avg_temp_k,
    reanalysis_dew_point_temp_k,
    reanalysis_relative_humidity_percent
  ) %>%
  na.omit()
dim(dat_model)
names(dat_model)

# Check the data
dim(dat_model)

head(dat_model)

summary(dat_model$total_cases)


# ============================================================
# 6. DGLM FUNCTION
# ============================================================

DGLM_function <- function(mean.formula,
                          dispersion.formula,
                          data,
                          maxit = 100,
                          epsilon = 1e-7) {
  
  # Make local copy
 
  
  # -----------------------------------------
  # Initial values
  # -----------------------------------------
  
  data_work$W_t <- 1
  
  dev_old <- Inf
  
  devdiff <- Inf
  
  iter <- 0
  
  
  # -----------------------------------------
  # Iterative DGLM algorithm
  # -----------------------------------------
  
  while (devdiff > epsilon && iter < maxit) {
    
    iter <- iter + 1
    
    
    # =======================================
    # STEP 1:
    # Fit mean model
    # =======================================
    
    ED_glm <- glm(
      formula = mean.formula,
      family = quasipoisson(link = "log"),
      data = data_work,
      weights = W_t
    )
    
    
    # Fitted mean
    mu <- fitted(ED_glm)
    
    
    # Observed dengue counts
    y <- data_work$total_cases
    
    
    # =======================================
    # STEP 2:
    # Calculate Poisson deviance
    # contribution for each observation
    # =======================================
    
    d <- ifelse(
      y == 0,
      2 * mu,
      2 * (
        y * log(y / mu) -
          (y - mu)
      )
    )
    
    
    # Gamma regression requires positive values
    d <- pmax(d, 1e-8)
    
    
    data_work$d <- d
    
    
    # =======================================
    # STEP 3:
    # Fit dispersion model
    # =======================================
    
    d_glm <- gam(
      formula = dispersion.formula,
      family = Gamma(link = "log"),
      data = data_work,
      method = "REML"
    )
    
    
    # Predicted dispersion phi_i
    phi <- fitted(d_glm)
    
    
    # Prevent numerical problems
    phi <- pmax(phi, 1e-8)
    
    
    # =======================================
    # STEP 4:
    # Variance
    #
    # Var(Y_i) = phi_i * mu_i
    # =======================================
    
    var_y <- phi * mu
    
    
    # =======================================
    # STEP 5:
    # Derivative of log link
    #
    # g(mu) = log(mu)
    #
    # g'(mu) = 1 / mu
    # =======================================
    
    d_gmu <- 1 / mu
    
    
    # =======================================
    # STEP 6:
    # Update weights
    # =======================================
    
    new_weights <- 1 / (
      var_y * (d_gmu^2)
    )
    
    
    # equivalent to approximately mu / phi
    data_work$W_t <- new_weights
    
    
    # =======================================
    # STEP 7:
    # Calculate total deviance
    # =======================================
    
    dev_new <- sum(d)
    
    
    # =======================================
    # STEP 8:
    # Check convergence
    # =======================================
    
    if (is.finite(dev_old)) {
      
      devdiff <- abs(
        dev_new - dev_old
      )
      
    } else {
      
      devdiff <- Inf
    }
    
    
    dev_old <- dev_new
    
    
    # Show iteration information
    cat(
      "Iteration:",
      iter,
      "| Deviance:",
      round(dev_new, 4),
      "| Difference:",
      round(devdiff, 8),
      "\n"
    )
    
  }
  
  
  # ==========================================================
  # Final model refit using final weights
  # ==========================================================
  
  ED_glm <- glm(
    formula = mean.formula,
    family = quasipoisson(link = "log"),
    data = data_work,
    weights = W_t
  )
  
  
  mu <- fitted(ED_glm)
  
  y <- data_work$total_cases
  
  
  d <- ifelse(
    y == 0,
    2 * mu,
    2 * (
      y * log(y / mu) -
        (y - mu)
    )
  )
  
  
  d <- pmax(d, 1e-8)
  
  data_work$d <- d
  
  
  d_glm <- gam(
    formula = dispersion.formula,
    family = Gamma(link = "log"),
    data = data_work,
    method = "REML"
  )
  
  
  # ==========================================================
  # Return results
  # ==========================================================
  
  output <- list(
    
    ED_glm = ED_glm,
    
    d_glm = d_glm,
    
    dev = sum(d),
    
    iter = iter,
    
    data = data_work
    
  )
  
  
  return(output)
  
}


# ============================================================
# 7. DEFINE MEAN MODEL
# ============================================================

model.mean <- total_cases ~
  
  ns(precipitation_amt_mm) +
  
  ns(
    reanalysis_avg_temp_k,
    df = 6
  ) +
  
  ns(
    reanalysis_dew_point_temp_k,
    df = 6
  ) +
  
  ns(
    reanalysis_relative_humidity_percent,
    df = 6
  ) +
  
  ns(
    time,
    df = 12
  ) +
  
  city +
  
  month


# Look at formula
model.mean


# ============================================================
# 8. DEFINE DISPERSION MODEL
# ============================================================

model.dispersion <- d ~
  
  city +
  
  s(
    time,
    k = 12
  )


# Look at formula
model.dispersion


# ============================================================
# 9. FIT DGLM
# ============================================================

fit_dglm <- DGLM_function(
  mean.formula = model.mean,
  dispersion.formula = model.dispersion,
  data = dat_model

  
)


# ============================================================
# 10. MEAN MODEL RESULTS
# ============================================================

summary(
  fit_dglm$ED_glm
)


# ============================================================
# 11. DISPERSION MODEL RESULTS
# ============================================================

summary(
  fit_dglm$d_glm
)


# ============================================================
# 12. NUMBER OF ITERATIONS
# ============================================================

fit_dglm$iter


# ============================================================
# 13. FINAL DEVIANCE
# ============================================================

fit_dglm$dev


# ============================================================
# 14. CREATE RESULTS DATASET
# ============================================================

dat_results <- fit_dglm$data


dat_results$fitted_cases <-
  fitted(
    fit_dglm$ED_glm
  )


dat_results$fitted_dispersion <-
  fitted(
    fit_dglm$d_glm
  )


head(
  dat_results
)


# ============================================================
# 15. SUMMARY OF ESTIMATED DISPERSION
# ============================================================



# ============================================================
# 16. PLOT DGLM DISPERSION OVER TIME
# ============================================================

plot(
  dat_results$week_start_date,
  dat_results$fitted_dispersion,
  type = "l",
  xlab = "Date",
  ylab = "Estimated dispersion",
  main = "Time-varying Dispersion of Dengue Cases"
)


# ============================================================
# 17. OBSERVED VS FITTED DENGUE CASES
# ============================================================

plot(
  dat_results$total_cases,
  dat_results$fitted_cases,
  xlab = "Observed Dengue Cases",
  ylab = "Fitted Dengue Cases",
  main = "Observed vs Fitted Dengue Cases"
)


abline(
  0,
  1,
  lty = 2
)


# ============================================================
# 18. ORDINARY QUASI-POISSON MODEL
#     FOR COMPARISON
# ============================================================

fit_qp <- glm(
  formula = model.mean,
  family = quasipoisson(link = "log"),
  data = dat_model
)


summary(
  fit_qp
)


# ============================================================
# 19. QUASI-POISSON CONSTANT DISPERSION
# ============================================================

qp_dispersion <- summary(
  fit_qp
)$dispersion


qp_dispersion


# ============================================================
# 20. DGLM DISPERSION SUMMARY
# ============================================================

dglm_dispersion <- fitted(
  fit_dglm$d_glm
)


summary(
  dglm_dispersion
)


# ============================================================
# 21. COMPARE DISPERSION
# ============================================================

cat(
  "\nQuasi-Poisson constant dispersion =",
  round(qp_dispersion, 3),
  "\n"
)


cat(
  "DGLM mean estimated dispersion =",
  round(mean(dglm_dispersion), 3),
  "\n"
)


cat(
  "DGLM minimum dispersion =",
  round(min(dglm_dispersion), 3),
  "\n"
)


cat(
  "DGLM maximum dispersion =",
  round(max(dglm_dispersion), 3),
  "\n"
)


# ============================================================
# 22. PLOT:
# QUASI-POISSON VS DGLM DISPERSION
# ============================================================

plot(
  dat_results$week_start_date,
  dat_results$fitted_dispersion,
  type = "l",
  xlab = "Date",
  ylab = "Estimated Dispersion",
  main = "Quasi-Poisson vs DGLM Dispersion"
)


abline(
  h = qp_dispersion,
  lty = 2
)


legend(
  "topright",
  legend = c(
    "DGLM time-varying dispersion",
    "Quasi-Poisson constant dispersion"
  ),
  lty = c(1, 2),
  bty = "n"
)

model.dispersion <- d ~
  city +
  s(time, k = 12)

summary(fit_dglm$d_glm)