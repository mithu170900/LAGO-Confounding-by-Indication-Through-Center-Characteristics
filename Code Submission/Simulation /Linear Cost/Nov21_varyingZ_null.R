# ==============================================================================
# NULL HYPOTHESIS SIMULATION - VARYING Z_J VERSION
# Z_j values are generated NEW for EACH simulation iteration
# ==============================================================================
library(dplyr)
library(sandwich)
library(LAGO)
library(lmtest)
library(parallel)
library(foreach)
library(doParallel)

start_time <- Sys.time()
source("OptInterLinGeneral_extraCVD.R")
source("OptInterLinGeneral_extraCVD 2.R")

# ==============================================================================
# SETUP PARALLEL PROCESSING
# ==============================================================================

n_cores <- detectCores() - 1
cat("Setting up parallel processing with", n_cores, "cores\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(dplyr)
  library(sandwich)
  library(LAGO)
  library(lmtest)
})

# ==============================================================================
# SETUP AND PARAMETERS - NULL HYPOTHESIS
# ==============================================================================
set.seed(10) 

# ==============================================================================
# SIMULATION PARAMETERS - NULL HYPOTHESIS
# ==============================================================================

# NULL HYPOTHESIS: No intervention effects
beta_true_int1 <- 0.00  # NULL: β₁* = 0
beta_true_int2 <- 0.00  # NULL: β₂* = 0
beta_z <- 2.42          # Keep confounder effect (same as main simulation)

# Target correlations (same as main simulation)
rho_1 <- 0.1  # Correlation between A1 and Z
rho_2 <- 0.2  # Correlation between A2 and Z

# Calculate eta values to achieve target correlations
# Since Z ~ N(0,1), sigma_Z = 1
# Formula: eta = rho / sqrt(1 - rho^2)
eta_1 <- rho_1 / sqrt(1 - rho_1^2)  # = 0.0501
eta_2 <- rho_2 / sqrt(1 - rho_2^2)  # = 0.0702

x_min <- c(0, 0)
x_max <- c(4, 3)

x1_rec_stage2 <- median(c(x_min[1], x_max[1]))
x2_rec_stage2 <- median(c(x_min[2], x_max[2]))

scenario_grid <- expand.grid(n1j = c(50, 100),     # Per center Stage 1
                             n2j = c(100, 200),    # Per center Stage 2
                             J = c(6, 10, 20))     # Number of centers
n_sims <- 2000
alpha <- 0.05  # Significance level

# ==============================================================================
# CALCULATE TRUE OPTIMAL UNDER NULL (using expected value of z_j = 0)
# ==============================================================================

my_beta <- c(beta_true_int1, beta_true_int2)  

trueXopt <- OptInterLinGeneral(
  beta.vec = -my_beta,  
  gamma.vec = -beta_z,
  lin.cost.coef = c(1, 0.5),
  p.bar = 5,
  x.min = x_min,
  x.max = x_max,
  z = 0,  # Using expected value of z_j
  intercept = FALSE
)$est.x.opt

cat("True optimal interventions (under null):", trueXopt, "\n")

clusterExport(cl, "trueXopt")

#======================
# Function for Stage 1 - NULL WITH 1:1 RANDOMIZATION
#======================

simulate_stage_1 <- function(n1j, J, z_j_current) {
  # Initialize data frame
  stage1_data <- data.frame()
  
  # Loop through each center
  for (j in 1:J) {
    # Generate participant data for center j
    # Half intervention (1), half control (0) - 1:1 ratio
    center_data <- data.frame(
      center = j,
      z_j = z_j_current[j],
      participant_id = 1:n1j,
      studyarm = rep(c(1, 0), each = n1j/2)  # 1 = intervention, 0 = control
    )
    
    # Initialize intervention columns
    center_data$int1 <- 0
    center_data$int2 <- 0
    
    # Generate interventions only for intervention arm (studyarm == 1)
    intervention_idx <- center_data$studyarm == 1
    n_intervention <- sum(intervention_idx)
    
    if(n_intervention > 0) {
      # Int1: int1 = x1_rec_stage2 + eta_1 * z_j + xi_ij
      center_data$int1[intervention_idx] <- x1_rec_stage2 + 
        eta_1 * z_j_current[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Int2: int2 = x2_rec_stage2 + eta_2 * z_j + xi_ij
      center_data$int2[intervention_idx] <- x2_rec_stage2 + 
        eta_2 * z_j_current[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Constrain to bounds 
      center_data$int1[intervention_idx] <- pmax(x_min[1], pmin(x_max[1], center_data$int1[intervention_idx]))
      center_data$int2[intervention_idx] <- pmax(x_min[2], pmin(x_max[2], center_data$int2[intervention_idx]))
    }
    
    # Control group (studyarm == 0) keeps int1 = 0, int2 = 0
    
    # Simulate outcomes under NULL HYPOTHESIS:
    # Y = 0*int1 + 0*int2 + beta_z * z_j + epsilon
    center_data$y <- beta_true_int1 * center_data$int1 +  # = 0
      beta_true_int2 * center_data$int2 +                 # = 0
      beta_z * center_data$z_j + 
      rnorm(n1j, mean = 0, sd = 1)
    
    # Add to overall dataset
    stage1_data <- rbind(stage1_data, center_data)
  }
  
  return(stage1_data)
}

#======================
# Function for Stage 2 - NULL WITH 1:1 RANDOMIZATION
#======================

simulate_stage_2 <- function(n2j, J, z_j_current, x1_rec_stage2, x2_rec_stage2) {
  # Initialize data frame
  stage2_data <- data.frame()
  
  # Loop through each center
  for (j in 1:J) {
    # Generate participant data for center j
    # Half intervention (1), half control (0) - 1:1 ratio
    center_data <- data.frame(
      center = j,
      z_j = z_j_current[j],
      participant_id = 1:n2j,
      studyarm = rep(c(1, 0), each = n2j/2)  # 1 = intervention, 0 = control
    )
    
    # Initialize intervention columns
    center_data$int1 <- 0
    center_data$int2 <- 0
    
    # Generate interventions only for intervention arm (studyarm == 1)
    intervention_idx <- center_data$studyarm == 1
    n_intervention <- sum(intervention_idx)
    
    if(n_intervention > 0) {
      # Int1: stage 2 recommended intervention comp from LAGO optimization at the end of stage 1
      center_data$int1[intervention_idx] <- x1_rec_stage2 + 
        eta_1 * z_j_current[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Int2: stage 2 recommended intervention comp from LAGO optimization at the end of stage 1
      center_data$int2[intervention_idx] <- x2_rec_stage2 + 
        eta_2 * z_j_current[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Constrain to bounds
      center_data$int1[intervention_idx] <- pmax(x_min[1], pmin(x_max[1], center_data$int1[intervention_idx]))
      center_data$int2[intervention_idx] <- pmax(x_min[2], pmin(x_max[2], center_data$int2[intervention_idx]))
    }
    
    # Control group (studyarm == 0) keeps int1 = 0, int2 = 0
    
    # Simulate outcomes under NULL HYPOTHESIS:
    # Y = 0*int1 + 0*int2 + beta_z * z_j + epsilon
    center_data$y <- beta_true_int1 * center_data$int1 +  # = 0
      beta_true_int2 * center_data$int2 +                 # = 0
      beta_z * center_data$z_j + 
      rnorm(n2j, mean = 0, sd = 1)
    
    # Add to overall dataset
    stage2_data <- rbind(stage2_data, center_data)
  }
  
  return(stage2_data)
}

# Export simulation functions to cluster
clusterExport(cl, c("simulate_stage_1", "simulate_stage_2"))

# Export global variables to cluster 
clusterExport(cl, c("beta_true_int1", "beta_true_int2", "beta_z",
                    "eta_1", "eta_2",
                    "x_min", "x_max", "x1_rec_stage2", "x2_rec_stage2",
                    "alpha"))

#======================
# Single simulation function - NULL HYPOTHESIS (for parallelization)
#======================

run_single_simulation <- function(sim_id, n1j, n2j, J) {
  # Set seed for reproducibility within each simulation
  set.seed(10 + sim_id * 1000 + n1j * 10 + n2j + J)
  
  tryCatch({
    # Generate new Z_j values for this simulation
    z_j_current <- rnorm(J, 0, 1)
    
    # Stage 1 simulation
    sim_data_stage1 <- simulate_stage_1(n1j, J, z_j_current)
    
    opt_stage1 <- lago_optimization(
      data = sim_data_stage1,
      outcome_name = "y",
      outcome_type = "continuous",
      glm_family = "gaussian",
      link = "identity",
      intervention_components = c("int1", "int2"),
      intervention_lower_bounds = x_min,
      intervention_upper_bounds = x_max,
      include_center_effects = TRUE,
      cost_list_of_vectors = list(c(0, 1), c(0, 0.5)),
      center_weights_for_outcome_goal = rep(1/J, J),
      outcome_goal = -5,
      outcome_goal_intention = "minimize",
      optimization_method = "grid_search",
      optimization_grid_search_step_size = c(0.1, 0.1),
      include_confidence_set = FALSE
    )
    
    x1_rec_stage2 <- opt_stage1$rec_int[1]
    x2_rec_stage2 <- opt_stage1$rec_int[2] 
    
    # Stage 1 model
    model_stage1 <- lm(y ~ int1 + int2 + z_j - 1, data = sim_data_stage1)
    stage1_coef_results <- coeftest(model_stage1, vcov = sandwich)
    
    # Stage 2 simulation (using same z_j_current)
    sim_data_stage2 <- simulate_stage_2(n2j, J, z_j_current, x1_rec_stage2, x2_rec_stage2)
    
    # Combine data
    sim_data_stage1$phase <- "Stage 1"
    sim_data_stage2$phase <- "Stage 2"
    
    common_cols <- c("center", "int1", "int2", "y", "z_j", "studyarm", "phase")
    combined_data <- rbind(
      sim_data_stage1[, common_cols],
      sim_data_stage2[, common_cols]
    )
    
    # Stage 2 optimization
    opt_stage2 <- lago_optimization(
      data = combined_data,
      outcome_name = "y",
      outcome_type = "continuous",
      glm_family = "gaussian",
      link = "identity",
      intervention_components = c("int1", "int2"),
      intervention_lower_bounds = x_min,
      intervention_upper_bounds = x_max,
      include_center_effects = TRUE,
      cost_list_of_vectors = list(c(0, 1), c(0, 0.5)),
      outcome_goal = -5,
      center_weights_for_outcome_goal = rep(1/J, J),
      outcome_goal_intention = "minimize",
      optimization_method = "grid_search",
      optimization_grid_search_step_size = c(0.1, 0.1),
      include_confidence_set = FALSE
    )
    
    # Combined model for hypothesis testing
    model_combined <- lm(y ~ int1 + int2 + z_j - 1, data = combined_data)
    combined_coef_results <- coeftest(model_combined, vcov = sandwich)
    
    # Extract p-values for hypothesis tests
    p_value_int1 <- combined_coef_results["int1", "Pr(>|t|)"]
    p_value_int2 <- combined_coef_results["int2", "Pr(>|t|)"]
    
    # Joint test using Wald test
    vcov_sandwich <- sandwich(model_combined)
    coef_indices <- which(names(coef(model_combined)) %in% c("int1", "int2"))
    beta_hat <- coef(model_combined)[coef_indices]
    vcov_beta <- vcov_sandwich[coef_indices, coef_indices]
    
    # Wald statistic for joint test H₀: β₁ = 0 & β₂ = 0
    wald_stat <- t(beta_hat) %*% solve(vcov_beta) %*% beta_hat
    p_value_joint <- pchisq(wald_stat, df = 2, lower.tail = FALSE)
    
    # Type I error indicators
    reject_h0_int1 <- as.numeric(p_value_int1 < alpha)
    reject_h0_int2 <- as.numeric(p_value_int2 < alpha)
    reject_h0_joint <- as.numeric(p_value_joint < alpha)
    
    # Return structured results
    return(list(
      sim_id = sim_id,
      
      # Coefficient estimates
      stage1_coefs = c(stage1_coef_results["int1", "Estimate"], 
                       stage1_coef_results["int2", "Estimate"]),
      combined_coefs = c(combined_coef_results["int1", "Estimate"], 
                         combined_coef_results["int2", "Estimate"]),
      
      # Standard errors
      stage1_ses = c(stage1_coef_results["int1", "Std. Error"], 
                     stage1_coef_results["int2", "Std. Error"]),
      combined_ses = c(combined_coef_results["int1", "Std. Error"], 
                       combined_coef_results["int2", "Std. Error"]),
      
      # P-values
      p_values = c(p_value_int1, p_value_int2, p_value_joint),
      
      # Type I error indicators
      type1_errors = c(reject_h0_int1, reject_h0_int2, reject_h0_joint),
      
      # Store z_j values used
      z_j_values = z_j_current,
      
      success = TRUE
    ))
    
  }, error = function(e) {
    return(list(
      sim_id = sim_id,
      success = FALSE,
      error = e$message
    ))
  })
}

# Export the single simulation function
clusterExport(cl, "run_single_simulation")

#======================
# Main simulation procedure with parallel processing
#======================

all_results <- list()

for (scenario_idx in 1:nrow(scenario_grid)) {
  current_n1j <- scenario_grid$n1j[scenario_idx]
  current_n2j <- scenario_grid$n2j[scenario_idx] 
  current_J <- scenario_grid$J[scenario_idx]
  
  cat("\n=== Processing Scenario", scenario_idx, "of", nrow(scenario_grid), "===\n")
  cat("n1j =", current_n1j, ", n2j =", current_n2j, ", J =", current_J, "\n")
  cat("Running", n_sims, "NULL HYPOTHESIS simulations in parallel...\n")
  
  # Run simulations in parallel
  sim_start_time <- Sys.time()
  
  parallel_results <- foreach(i = 1:n_sims, 
                              .packages = c("dplyr", "sandwich", "LAGO", "lmtest"),
                              .errorhandling = 'pass') %dopar% {
                                run_single_simulation(i, current_n1j, current_n2j, current_J)
                              }
  
  sim_end_time <- Sys.time()
  sim_time <- difftime(sim_end_time, sim_start_time, units = "secs")
  cat("  Completed", n_sims, "simulations in", round(sim_time, 2), "seconds\n")
  
  # Filter successful results
  successful_results <- Filter(function(x) isTRUE(x$success), parallel_results)
  n_sims_actual <- length(successful_results)
  
  if(n_sims_actual < n_sims) {
    cat("  Warning:", n_sims - n_sims_actual, "simulations failed\n")
  }
  
  # Organize results
  stage1_coefs <- t(sapply(successful_results, function(x) x$stage1_coefs))
  combined_coefs <- t(sapply(successful_results, function(x) x$combined_coefs))
  stage1_ses <- t(sapply(successful_results, function(x) x$stage1_ses))
  combined_ses <- t(sapply(successful_results, function(x) x$combined_ses))
  p_values_matrix <- t(sapply(successful_results, function(x) x$p_values))
  type1_errors_matrix <- t(sapply(successful_results, function(x) x$type1_errors))
  
  # Extract z_j values used in each simulation
  z_j_values <- lapply(successful_results, function(x) x$z_j_values)
  
  all_results[[scenario_idx]] <- list(
    scenario_params = scenario_grid[scenario_idx, ],
    stage1_coefs = stage1_coefs,
    combined_coefs = combined_coefs,
    stage1_ses = stage1_ses,
    combined_ses = combined_ses,
    p_values_matrix = p_values_matrix,
    type1_errors_matrix = type1_errors_matrix,
    z_j_values = z_j_values,  # Store z_j values for each simulation
    n_sims_actual = n_sims_actual
  )
  
  cat("  Completed scenario", scenario_idx, "\n")
}

# Stop the cluster
stopCluster(cl)

cat("\n=== ALL NULL HYPOTHESIS SIMULATIONS COMPLETED ===\n")

# ==============================================================================
# NULL HYPOTHESIS ANALYSIS
# ==============================================================================

cat("=== CALCULATING NULL HYPOTHESIS METRICS ===\n")

null_results <- data.frame()

for (scenario_idx in 1:length(all_results)) {
  result <- all_results[[scenario_idx]]
  scenario_params <- result$scenario_params
  
  # Bias calculations (multiply by 1000 for reporting)
  bias_beta1_combined <- (mean(result$combined_coefs[, 1]) - beta_true_int1) * 1000
  bias_beta2_combined <- (mean(result$combined_coefs[, 2]) - beta_true_int2) * 1000
  
  # Type I error rates
  type1_error_alpha1 <- mean(result$type1_errors_matrix[, 1])
  type1_error_alpha2 <- mean(result$type1_errors_matrix[, 2])
  type1_error_combined <- mean(result$type1_errors_matrix[, 3])
  
  # Mean p-values
  mean_p_value_int1 <- mean(result$p_values_matrix[, 1])
  mean_p_value_int2 <- mean(result$p_values_matrix[, 2])
  mean_p_value_joint <- mean(result$p_values_matrix[, 3])
  
  # SE/EMP.SD ratios
  se_emp_sd_beta1 <- mean(result$combined_ses[, 1]) / sd(result$combined_coefs[, 1]) * 100
  se_emp_sd_beta2 <- mean(result$combined_ses[, 2]) / sd(result$combined_coefs[, 2]) * 100
  
  null_results <- rbind(null_results, data.frame(
    n1j = scenario_params$n1j,
    n2j = scenario_params$n2j,
    J = scenario_params$J,
    se_emp_sd_beta1 = se_emp_sd_beta1,
    se_emp_sd_beta2 = se_emp_sd_beta2,
    bias_beta1 = bias_beta1_combined,
    bias_beta2 = bias_beta2_combined,
    alpha_1 = type1_error_alpha1,
    alpha_2 = type1_error_alpha2,
    alpha_combined = type1_error_combined,
    mean_p_int1 = mean_p_value_int1,
    mean_p_int2 = mean_p_value_int2,
    mean_p_joint = mean_p_value_joint
  ))
}

null_results <- null_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\n=== NULL HYPOTHESIS SIMULATION RESULTS ===\n")
cat("True β₁* =", beta_true_int1, ", True β₂* =", beta_true_int2, "\n")
cat("Expected Type I error rate (α) =", alpha, "\n\n")

print(null_results)

# Type I error assessment
cat("\n=== TYPE I ERROR ASSESSMENT ===\n")
acceptable_lower <- alpha - 0.02  
acceptable_upper <- alpha + 0.02

within_bounds <- null_results %>%
  summarise(
    alpha1_ok = sum(alpha_1 >= acceptable_lower & alpha_1 <= acceptable_upper),
    alpha2_ok = sum(alpha_2 >= acceptable_lower & alpha_2 <= acceptable_upper),
    combined_ok = sum(alpha_combined >= acceptable_lower & alpha_combined <= acceptable_upper)
  )

cat("Scenarios with acceptable Type I error rates [", acceptable_lower, ",", acceptable_upper, "]:\n")
cat("  α₁ (int1):", within_bounds$alpha1_ok, "/", nrow(null_results), "\n")
cat("  α₂ (int2):", within_bounds$alpha2_ok, "/", nrow(null_results), "\n") 
cat("  α_combined:", within_bounds$combined_ok, "/", nrow(null_results), "\n")

# Overall summary
cat("\n=== OVERALL SUMMARY ===\n")
overall_summary <- null_results %>%
  summarise(
    mean_bias_beta1 = mean(bias_beta1),
    mean_bias_beta2 = mean(bias_beta2),
    mean_alpha1 = mean(alpha_1),
    mean_alpha2 = mean(alpha_2),
    mean_alpha_combined = mean(alpha_combined),
    mean_se_ratio_beta1 = mean(se_emp_sd_beta1),
    mean_se_ratio_beta2 = mean(se_emp_sd_beta2)
  ) %>%
  mutate(across(everything(), ~ round(.x, 4)))

print(overall_summary)

# Save results
write.csv(null_results, "Table_NULL_linear_varying_zj_highconf.csv", row.names = FALSE)

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")
cat("\nNull hypothesis simulation time:", round(total_time, 2), "minutes\n")

cat("\n=== SIMULATION COMPLETE ===\n")
cat("Results saved to: Table_NULL_linear_varying_zj.csv\n")

# ==============================================================================
# VERIFY CORRELATIONS ACROSS ALL SIMULATIONS
# ==============================================================================
cat("\n=== VERIFYING CORRELATIONS ACROSS ALL SIMULATIONS ===\n")

all_cor_int1_z <- numeric(length(all_results) * n_sims)
all_cor_int2_z <- numeric(length(all_results) * n_sims)
counter <- 1

# Note: We don't have combined_data stored, so we'll just report the expected correlations
cat("\nExpected Correlations (based on simulation setup):\n")
cat("Cor(int1, z_j): Target =", round(rho_1, 4), "\n")
cat("Cor(int2, z_j): Target =", round(rho_2, 4), "\n")
cat("\nNote: Z_j values vary per simulation iteration in this design.\n")