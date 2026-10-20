# ============================================================
# DENGUE DGLM ANALYSIS
# Time-varying dispersion model
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(mgcv)
library(splines)
library(lubridate)



# ------------------------------------------------------------
# 2. Read the dengue datasets
# ------------------------------------------------------------

features <- read.csv(
  "https://s3.amazonaws.com/drivendata/data/44/public/dengue_features_train.csv"
)

labels <- read.csv(
  "https://s3.amazonaws.com/drivendata/data/44/public/dengue_labels_train.csv"
)




# ------------------------------------------------------------
# 3. Merge features with dengue case counts
# ------------------------------------------------------------

dat <- features %>%
  left_join(
    labels,
    by = c("city", "year", "weekofyear")
  )



# ------------------------------------------------------------
# 4. Select San Juan
# ------------------------------------------------------------

dat <- dat %>%
  filter(city == "sj") %>%
  arrange(week_start_date)



# ------------------------------------------------------------
# 5. Create time variables
# ------------------------------------------------------------

dat <- dat %>%
  mutate(
    
    # Convert date variable to R Date format
    week_start_date = as.Date(week_start_date),
    
    # Month: 1 = January, ..., 12 = December
    month = factor(month(week_start_date)),
    
    # Sequential week number
    time_index = row_number()
  )



# ------------------------------------------------------------
# 6. Keep variables needed for the model
#    and remove observations with missing values
# ------------------------------------------------------------

dat_model <- dat %>%
  drop_na(
    total_cases,
    precipitation_amt_mm,
    reanalysis_avg_temp_k,
    reanalysis_dew_point_temp_k,
    month,
    time_index
  )



# ------------------------------------------------------------
# 7. Check the data
# ------------------------------------------------------------

glimpse(dat_model)

summary(dat_model$total_cases)

nrow(dat_model)



# ============================================================
# 8. CREATE DGLM FUNCTION
# ============================================================

DGLM_function <- function(mean.formula,
                          dispersion.pred,
                          data,
                          method = "gam") {
  
  
  # ----------------------------------------------------------
  # Convergence settings
  # ----------------------------------------------------------
  
  # Small number used to determine convergence
  epsilon <- 1e-07
  
  # Maximum number of iterations
  maxit <- 100
  
  
  # ----------------------------------------------------------
  # Initial values
  # ----------------------------------------------------------
  
  iter <- 0
  
  dev <- 0
  
  devdiff <- 1
  
  # Initial weights
  data$W_t <- 1
  
  
  
  # ----------------------------------------------------------
  # Initial mean model
  # ----------------------------------------------------------
  
  ini.model <- glm(
    formula = mean.formula,
    family = quasipoisson(link = "log"),
    data = data
  )
  
  
  
  # ==========================================================
  # Iterative DGLM procedure
  # ==========================================================
  
  while (devdiff > epsilon && iter < maxit) {
    
    
    # --------------------------------------------------------
    # A. Fit the mean model
    # --------------------------------------------------------
    
    ED_glm <- update(
      ini.model,
      weights = W_t,
      data = data
    )
    
    
    # Observed dengue cases
    y <- data$total_cases
    
    
    # Estimated mean dengue cases
    mu <- ED_glm$fitted.values
    
    
    
    # --------------------------------------------------------
    # B. Calculate Poisson deviance contribution
    # --------------------------------------------------------
    
    # Start with the expression for observations where y = 0
    d <- 2 * mu
    
    
    # For observations where dengue cases > 0
    positive <- y > 0
  
    
    d[positive] <- 2 * (
      
      y[positive] *
        log(
          y[positive] /
            mu[positive]
        ) -
        
        (
          y[positive] -
            mu[positive]
        )
    )
    
    
    # Gamma regression requires positive response values
    d <- pmax(d, 1e-08)
    
    
    # Add deviance contribution to dataset
    data$d <- d
    
    
    
    # --------------------------------------------------------
    # C. Create dispersion model formula
    # --------------------------------------------------------
    
    dispersion.formula <- reformulate(
      dispersion.pred,
      response = "d"
    )
    
    
    
    # --------------------------------------------------------
    # D. Fit dispersion model
    # --------------------------------------------------------
    
    if (method == "gam") {
      
      
      d_glm <- gam(
        dispersion.formula,
        family = Gamma(link = "log"),
        data = data,
        select = TRUE,
        method = "REML"
      )
      
      
    } else if (method == "glm") {
      
      
      d_glm <- glm(
        dispersion.formula,
        family = Gamma(link = "log"),
        data = data
      )
      
      
    } else {
      
      stop("method must be either 'gam' or 'glm'")
      
    }
    
    
    
    # --------------------------------------------------------
    # E. Estimate variance
    #
    # Var(Y_i) = phi_i * mu_i
    #
    # phi_i = estimated dispersion
    # mu_i  = estimated mean
    # --------------------------------------------------------
    
    var_y <- d_glm$fitted.values *
      ED_glm$fitted.values
    
    
    
    # --------------------------------------------------------
    # F. Derivative of link function
    #
    # g(mu) = log(mu)
    #
    # g'(mu) = 1/mu
    # --------------------------------------------------------
    
    d_gmu <- 1 /
      ED_glm$fitted.values
    
    
    
    # --------------------------------------------------------
    # G. Update weights
    #
    # W_i =
    # 1 / [Var(Y_i) * {g'(mu_i)}^2]
    # --------------------------------------------------------
    
    data$W_t <- 1 /
      (
        var_y *
          (d_gmu^2)
      )
    
    
    
    # --------------------------------------------------------
    # H. Calculate total deviance
    # --------------------------------------------------------
    
    dev_components <- 2 *
      ED_glm$fitted.values
    
    
    dev_components[positive] <- 2 * (
      
      y[positive] *
        log(
          y[positive] /
            ED_glm$fitted.values[positive]
        ) -
        
        (
          y[positive] -
            ED_glm$fitted.values[positive]
        )
    )
    
    
    dev_t <- sum(dev_components)
    
    
    
    # --------------------------------------------------------
    # I. Check convergence
    # --------------------------------------------------------
    
    devdiff <- abs(
      dev_t - dev
    )
    
    
    dev <- dev_t
    
    
    iter <- iter + 1
    
  }
  
  
  
  # ----------------------------------------------------------
  # Save final results
  # ----------------------------------------------------------
  
  output <- list(
    
    # Mean model
    ED_glm = ED_glm,
    
    # Dispersion model
    d_glm = d_glm,
    
    # Final deviance
    dev = dev,
    
    # Number of iterations
    iter = iter,
    
    # Dataset containing final weights/deviance
    data = data
    
  )
  
  
  return(output)
}



# ============================================================
# 9. SPECIFY THE MEAN MODEL
# ============================================================


# Number of years in San Juan data
n_years <- n_distinct(dat_model$year)

n_years


# Degrees of freedom for long-term time trend
# Starting choice: 4 df per year
time_df <- 4 * n_years

time_df



# Variables in the mean model
variables <- c(
  
  # Rainfall
  "precipitation_amt_mm",
  
  # Nonlinear temperature effect
  "ns(reanalysis_avg_temp_k, 6)",
  
  # Nonlinear dew-point/humidity effect
  "ns(reanalysis_dew_point_temp_k, 6)",
  
  # Long-term temporal trend
  paste0(
    "ns(time_index, ",
    time_df,
    ")"
  )
)



# Create mean model formula
model.mean <- reformulate(
  variables,
  response = "total_cases"
)


# Display formula
model.mean



# ============================================================
# 10. FIT DGLM
# ============================================================

fit <- DGLM_function(
  
  mean.formula = model.mean,
  
  # Allow dispersion to change according to month
  dispersion.pred = "month",
  
  data = dat_model,
  
  method = "gam"
)



# ============================================================
# 11. CHECK CONVERGENCE
# ============================================================

fit$iter

fit$dev



# ============================================================
# 12. MEAN MODEL RESULTS
# ============================================================

summary(fit$ED_glm)



# ============================================================
# 13. DISPERSION MODEL RESULTS
# ============================================================

summary(fit$d_glm)



# ============================================================
# 14. SAVE FITTED MEAN AND DISPERSION
# ============================================================

dat_results <- fit$data %>%
  mutate(
    
    # Estimated number of dengue cases
    fitted_cases = fitted(fit$ED_glm),
    
    # Estimated dispersion
    fitted_dispersion = fitted(fit$d_glm)
  )



# Look at results
dat_results %>%
  select(
    week_start_date,
    year,
    weekofyear,
    total_cases,
    fitted_cases,
    fitted_dispersion
  ) %>%
  head(20)



# ============================================================
# 15. PLOT OBSERVED AND FITTED DENGUE CASES
# ============================================================

ggplot(
  dat_results,
  aes(x = week_start_date)
) +
  
  geom_line(
    aes(
      y = total_cases,
      linetype = "Observed cases"
    )
  ) +
  
  geom_line(
    aes(
      y = fitted_cases,
      linetype = "Fitted cases"
    )
  ) +
  
  labs(
    x = "Date",
    y = "Weekly dengue cases",
    linetype = NULL,
    title = "Observed and Fitted Dengue Cases: San Juan"
  ) +
  
  theme_minimal()



# ============================================================
# 16. PLOT ESTIMATED DISPERSION OVER TIME
# ============================================================

ggplot(
  dat_results,
  aes(
    x = week_start_date,
    y = fitted_dispersion
  )
) +
  
  geom_line() +
  
  labs(
    x = "Date",
    y = "Estimated dispersion",
    title = "Estimated Time-Varying Dispersion"
  ) +
  
  theme_minimal()



# ============================================================
# 17. DISPERSION BY MONTH
# ============================================================

dat_results %>%
  group_by(month) %>%
  summarise(
    mean_dispersion = mean(
      fitted_dispersion,
      na.rm = TRUE
    ),
    .groups = "drop"
  )



# ============================================================
# 18. PLOT DISPERSION BY MONTH
# ============================================================

dat_results %>%
  group_by(month) %>%
  summarise(
    mean_dispersion = mean(
      fitted_dispersion,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  ggplot(
    aes(
      x = month,
      y = mean_dispersion,
      group = 1
    )
  ) +
  
  geom_line() +
  geom_point() +
  
  labs(
    x = "Month",
    y = "Estimated dispersion",
    title = "Monthly Variation in Dengue Dispersion"
  ) +
  
  theme_minimal()
# ============================================================
# QUASI-POISSON VS DGLM COMPARISON
#
# Run this AFTER running the complete DGLM code.
#
# Existing objects:
#   dat_model = dengue analysis dataset
#   model.mean = mean model formula
#   fit = fitted DGLM
# ============================================================


# ------------------------------------------------------------
# 1. Rename the existing DGLM fit
# ------------------------------------------------------------

fit_dglm <- fit



# ============================================================
# 2. FIT STANDALONE QUASI-POISSON MODEL
# ============================================================

fit_qp <- glm(
  formula = model.mean,
  family = quasipoisson(link = "log"),
  data = dat_model
)



# ------------------------------------------------------------
# 3. Quasi-Poisson model results
# ------------------------------------------------------------

summary(fit_qp)



# ------------------------------------------------------------
# 4. DGLM mean-model results
# ------------------------------------------------------------

summary(fit_dglm$ED_glm)



# ------------------------------------------------------------
# 5. DGLM dispersion-model results
# ------------------------------------------------------------

summary(fit_dglm$d_glm)



# ============================================================
# 6. COMPARE DISPERSION
# ============================================================


# Quasi-Poisson:
# one constant dispersion value for all weeks

phi_qp <- summary(fit_qp)$dispersion

phi_qp



# DGLM:
# different estimated dispersion values over time

phi_dglm <- fitted(fit_dglm$d_glm)

summary(phi_dglm)

range(phi_dglm)

sd(phi_dglm)



# ============================================================
# 7. CREATE COMPARISON DATASET
# ============================================================

comparison_data <- dat_model %>%
  mutate(
    
    # Observed dengue cases
    observed = total_cases,
    
    # Fitted mean from Quasi-Poisson
    fitted_qp = fitted(fit_qp),
    
    # Fitted mean from DGLM
    fitted_dglm = fitted(fit_dglm$ED_glm),
    
    # Constant Quasi-Poisson dispersion
    phi_qp = phi_qp,
    
    # Time-varying DGLM dispersion
    phi_dglm = phi_dglm
  )



# Look at first observations

comparison_data %>%
  select(
    week_start_date,
    observed,
    fitted_qp,
    fitted_dglm,
    phi_qp,
    phi_dglm
  ) %>%
  head(20)



# ============================================================
# 8. PLOT QUASI-POISSON VS DGLM DISPERSION
# ============================================================

ggplot(
  comparison_data,
  aes(
    x = week_start_date,
    y = phi_dglm
  )
) +
  
  geom_line() +
  
  geom_hline(
    yintercept = phi_qp,
    linetype = "dashed"
  ) +
  
  labs(
    x = "Date",
    y = "Estimated dispersion",
    title = "Quasi-Poisson vs DGLM Dispersion",
    subtitle = "Dashed line = constant Quasi-Poisson dispersion"
  ) +
  
  theme_minimal()



# ============================================================
# 9. CALCULATE STANDARDIZED RESIDUALS
# ============================================================


# Observed dengue cases

y <- dat_model$total_cases



# ------------------------------------------------------------
# Quasi-Poisson residuals
#
# r = (Y - mu) / sqrt(phi * mu)
# ------------------------------------------------------------

mu_qp <- fitted(fit_qp)

resid_qp <- (
  y - mu_qp
) /
  sqrt(
    phi_qp * mu_qp
  )



# ------------------------------------------------------------
# DGLM residuals
#
# r = (Y - mu) / sqrt(phi_t * mu)
# ------------------------------------------------------------

mu_dglm <- fitted(fit_dglm$ED_glm)

resid_dglm <- (
  y - mu_dglm
) /
  sqrt(
    phi_dglm * mu_dglm
  )



# Add residuals to comparison dataset

comparison_data <- comparison_data %>%
  mutate(
    resid_qp = resid_qp,
    resid_dglm = resid_dglm
  )



# ============================================================
# 10. COMPARE OVERALL RESIDUAL VARIANCE
# ============================================================


residual_variance_comparison <- tibble(
  
  Model = c(
    "Quasi-Poisson",
    "DGLM"
  ),
  
  Standardized_Residual_Variance = c(
    var(resid_qp, na.rm = TRUE),
    var(resid_dglm, na.rm = TRUE)
  )
)


residual_variance_comparison