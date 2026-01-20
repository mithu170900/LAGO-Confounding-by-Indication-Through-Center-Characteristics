# ==============================================================================
# TWO-STAGE SIMULATION - FIXED Z_J VERSION WITH CUBIC COSTS
# Z_j values are generated ONCE and remain CONSTANT across all simulations
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
# SETUP AND PARAMETERS
# ==============================================================================
set.seed(10) 

center_A <- -2.63
center_B <- 0.58
center_C <- 2.11
avg_center_effect <- mean(c(center_A, center_B, center_C))

# Generate center-level parameters (FIXED across all simulations)
z_j <- rnorm(20, 0, 1)

# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

beta_true_int1 <- -1.70
beta_true_int2 <- -0.70
beta_z <- 2.42  

# Target correlations
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

# Cubic cost structure
cost_list <- list(
  c(0, 1.25, 0, -0.04, 0.0055),
  c(0, 0.63, 0, -0.09, 0.026)
)

# ==============================================================================
# CALCULATE TRUE OPTIMAL (using expected value of z_j = 0)
# ==============================================================================

# Function to calculate true optimum with CUBIC costs using grid search
calculate_cubic_true_optimum <- function(true_beta1, true_beta2, x.min, x.max) {
  # Create fine grid for optimization (matching LAGO's approach)
  step_size <- 0.01
  x1_grid <- seq(x.min[1], x.max[1], by = step_size)
  x2_grid <- seq(x.min[2], x.max[2], by = step_size)
  grid <- expand.grid(x1 = x1_grid, x2 = x2_grid)
  
  # Calculate outcomes using TRUE coefficients with effective intercept (E[z_j] = 0)
  grid$outcome <- true_beta1 * grid$x1 + true_beta2 * grid$x2 + beta_z * 0
  
  # Calculate CUBIC costs (matching LAGO exactly)
  grid$cost <- (1.25 * grid$x1 - 0.04 * grid$x1^3 + 0.0055 * grid$x1^4) +
    (0.63 * grid$x2 - 0.09 * grid$x2^3 + 0.026 * grid$x2^4)
  
  # Find feasible solutions (outcome <= -5) and get the minimum cost solution
  feasible <- grid[grid$outcome <= -5, ]
  
  if(nrow(feasible) == 0) {
    # If no feasible solution, return boundary values
    return(c(x.max[1], x.max[2]))
  }
  
  best <- feasible[which.min(feasible$cost), ]
  return(c(best$x1, best$x2))
}

# Calculate true optimum
trueXopt <- calculate_cubic_true_optimum(beta_true_int1, beta_true_int2, x_min, x_max)
cat("True optimal interventions:", trueXopt, "\n")

#======================
# Function for Stage 1
#======================

simulate_stage_1 <- function(n1j, J) {
  # Initialize data frame
  stage1_data <- data.frame()
  
  z_j_subset <- z_j[1:J]
  
  # Loop through each center
  for (j in 1:J) {
    # Generate participant data for center j
    # Half intervention (1), half control (0) - 1:1 ratio
    center_data <- data.frame(
      center = j,
      z_j = z_j_subset[j],
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
        eta_1 * z_j_subset[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Int2: int2 = x2_rec_stage2 + eta_2 * z_j + xi_ij
      center_data$int2[intervention_idx] <- x2_rec_stage2 + 
        eta_2 * z_j_subset[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Constrain to bounds 
      center_data$int1[intervention_idx] <- pmax(x_min[1], pmin(x_max[1], center_data$int1[intervention_idx]))
      center_data$int2[intervention_idx] <- pmax(x_min[2], pmin(x_max[2], center_data$int2[intervention_idx]))
    }
    
    # Control group (studyarm == 0) keeps int1 = 0, int2 = 0
    
    # Simulate outcomes for all participants:
    # Y = beta_true_int1 * int1 + beta_true_int2 * int2 + beta_z * z_j + epsilon
    center_data$y <- beta_true_int1 * center_data$int1 + 
      beta_true_int2 * center_data$int2 + 
      beta_z * center_data$z_j + 
      rnorm(n1j, mean = 0, sd = 1)
    
    # Add to overall dataset
    stage1_data <- rbind(stage1_data, center_data)
  }
  
  return(stage1_data)
}

#======================
# Function for Stage 2
#======================

simulate_stage_2 <- function(n2j, J, x1_rec_stage2, x2_rec_stage2) {
  # Select first J centers (same centers as Stage 1)
  z_j_subset <- z_j[1:J]
  
  # Initialize data frame
  stage2_data <- data.frame()
  
  # Loop through each center
  for (j in 1:J) {
    # Generate participant data for center j
    # Half intervention (1), half control (0) - 1:1 ratio
    center_data <- data.frame(
      center = j,
      z_j = z_j_subset[j],
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
        eta_1 * z_j_subset[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Int2: stage 2 recommended intervention comp from LAGO optimization at the end of stage 1
      center_data$int2[intervention_idx] <- x2_rec_stage2 + 
        eta_2 * z_j_subset[j] + 
        rnorm(n_intervention, mean = 0, sd = 1)
      
      # Constrain to bounds
      center_data$int1[intervention_idx] <- pmax(x_min[1], pmin(x_max[1], center_data$int1[intervention_idx]))
      center_data$int2[intervention_idx] <- pmax(x_min[2], pmin(x_max[2], center_data$int2[intervention_idx]))
    }
    
    # Control group (studyarm == 0) keeps int1 = 0, int2 = 0
    
    # Simulate outcomes for all participants:
    # Y = beta_true_int1 * int1 + beta_true_int2 * int2 + beta_z * z_j + epsilon
    center_data$y <- beta_true_int1 * center_data$int1 + 
      beta_true_int2 * center_data$int2 + 
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
                    "eta_1", "eta_2", "z_j",
                    "x_min", "x_max", "x1_rec_stage2", "x2_rec_stage2", 
                    "cost_list", "trueXopt"))

#======================
# Single simulation function (for parallelization)
#======================

run_single_simulation <- function(sim_id, n1j, n2j, J) {
  # Set seed for reproducibility within each simulation
  set.seed(10 + sim_id * 1000 + n1j * 10 + n2j + J)
  
  # Stage 1 simulation
  sim_data_stage1 <- simulate_stage_1(n1j, J)
  
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
    cost_list_of_vectors = cost_list,
    outcome_goal = -5,
    center_weights_for_outcome_goal = rep(1/J, J),
    outcome_goal_intention = "minimize",
    optimization_method = "grid_search",
    optimization_grid_search_step_size = c(0.1, 0.1),
    include_confidence_set = FALSE
  )
  
  x1_rec_stage2 <- opt_stage1$rec_int[1]
  x2_rec_stage2 <- opt_stage1$rec_int[2] 
  
  # Stage 1 model
  model_stage1 <- lm(y ~ int1 + int2 + z_j - 1, data = sim_data_stage1)
  stage1_results <- list(
    coef = coeftest(model_stage1, vcov = sandwich),
    ci = coefci(model_stage1, vcov. = sandwich, level = 0.95),
    x1_rec_stage2 = x1_rec_stage2,
    x2_rec_stage2 = x2_rec_stage2
  )
  
  # Stage 2 simulation
  sim_data_stage2 <- simulate_stage_2(n2j, J, x1_rec_stage2, x2_rec_stage2)
  
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
    cost_list_of_vectors = cost_list,
    outcome_goal = -5,
    center_weights_for_outcome_goal = rep(1/J, J),
    outcome_goal_intention = "minimize",
    optimization_method = "grid_search",
    optimization_grid_search_step_size = c(0.1, 0.1),
    include_confidence_set = FALSE
  )
  
  x1_final_est_opt <- opt_stage2$rec_int[1]
  x2_final_est_opt <- opt_stage2$rec_int[2]
  
  # Combined model
  model_combined <- lm(y ~ int1 + int2 + z_j - 1, data = combined_data)
  combined_results <- list(
    coef = coeftest(model_combined, vcov = sandwich),
    ci = coefci(model_combined, vcov. = sandwich, level = 0.95),
    data = combined_data
  )
  
  return(list(
    x_rec_stage2 = c(x1_rec_stage2 = x1_rec_stage2, x2_rec_stage2 = x2_rec_stage2),
    x_final_est_opt = c(x1_final_est_opt = x1_final_est_opt, x2_final_est_opt = x2_final_est_opt),
    stage1_results = stage1_results,
    combined_results = combined_results
  ))
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
  cat("Running", n_sims, "simulations in parallel...\n")
  
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
  
  # Extract results
  x_rec_stage2 <- data.frame(
    x1_rec_stage2 = sapply(parallel_results, function(x) x$x_rec_stage2["x1_rec_stage2"]),
    x2_rec_stage2 = sapply(parallel_results, function(x) x$x_rec_stage2["x2_rec_stage2"])
  )
  
  x_final_est_opt <- data.frame(
    x1_final_est_opt = sapply(parallel_results, function(x) x$x_final_est_opt["x1_final_est_opt"]),
    x2_final_est_opt = sapply(parallel_results, function(x) x$x_final_est_opt["x2_final_est_opt"])
  )
  
  stage1_sim_results <- lapply(parallel_results, function(x) x$stage1_results)
  combined_sim_results <- lapply(parallel_results, function(x) x$combined_results)
  
  all_results[[scenario_idx]] <- list(
    scenario_params = scenario_grid[scenario_idx, ],
    x_rec_stage2 = x_rec_stage2,
    x_final_est_opt = x_final_est_opt,
    stage1_sim_results = stage1_sim_results,
    combined_sim_results = combined_sim_results
  )
  
  cat("  Completed scenario", scenario_idx, "\n")
}

# Stop the cluster
stopCluster(cl)

cat("\n=== ALL SIMULATIONS COMPLETED ===\n")

# ==============================================================================
# METRICS CALCULATION AND TABLE CREATION
# ==============================================================================

cat("=== CALCULATING METRICS ACROSS ALL SCENARIOS ===\n")

# Function to calculate metrics for a single scenario
calculate_scenario_metrics <- function(scenario_idx, all_results, trueXopt, beta_true_int1, beta_true_int2, beta_z, z_j, x_min, x_max, n_sims) {
  result <- all_results[[scenario_idx]]
  scenario_params <- result$scenario_params
  n_centers <- scenario_params$J
  stage1_n <- scenario_params$n1j
  stage2_n <- scenario_params$n2j
  combined_n_total <- stage1_n + stage2_n
  
  # Table A metrics
  betahat_stage2 <- cbind(
    sapply(result$combined_sim_results, function(x) x$coef["int1", "Estimate"]),
    sapply(result$combined_sim_results, function(x) x$coef["int2", "Estimate"])
  )
  
  sehat_stage2 <- cbind(
    sapply(result$combined_sim_results, function(x) x$coef["int1", "Std. Error"]),
    sapply(result$combined_sim_results, function(x) x$coef["int2", "Std. Error"])
  )
  
  z_score <- 1.96
  beta1ci <- cbind(
    betahat_stage2[,1] - z_score * sehat_stage2[,1],
    betahat_stage2[,1] + z_score * sehat_stage2[,1]
  )
  beta2ci <- cbind(
    betahat_stage2[,2] - z_score * sehat_stage2[,2],
    betahat_stage2[,2] + z_score * sehat_stage2[,2]
  )
  
  final_relbias_int1 <- mean(100 * (betahat_stage2[,1] - beta_true_int1) / abs(beta_true_int1))
  final_se_ratio_int1 <- mean(sehat_stage2[,1]) / sd(betahat_stage2[,1]) * 100
  final_coverage_int1 <- sum((beta1ci[,1] <= beta_true_int1) & (beta1ci[,2] >= beta_true_int1)) / n_sims * 100
  
  final_relbias_int2 <- mean(100 * (betahat_stage2[,2] - beta_true_int2) / abs(beta_true_int2))
  final_se_ratio_int2 <- mean(sehat_stage2[,2]) / sd(betahat_stage2[,2]) * 100
  final_coverage_int2 <- sum((beta2ci[,1] <= beta_true_int2) & (beta2ci[,2] >= beta_true_int2)) / n_sims * 100
  
  # Table B metrics
  x_rec_stage2 <- result$x_rec_stage2
  x_final_est_opt <- result$x_final_est_opt
  
  stage2_bias_int1 <- mean(x_rec_stage2$x1_rec_stage2) - trueXopt[1]
  stage2_bias_int2 <- mean(x_rec_stage2$x2_rec_stage2) - trueXopt[2]
  stage2_euclidean_distances <- apply(x_rec_stage2, 1, function(row) sqrt(sum((trueXopt - row)^2)))
  stage2_overall_rmse <- mean(stage2_euclidean_distances)
  
  final_bias_int1 <- mean(x_final_est_opt$x1_final_est_opt) - trueXopt[1]
  final_bias_int2 <- mean(x_final_est_opt$x2_final_est_opt) - trueXopt[2]
  final_euclidean_distances <- apply(x_final_est_opt, 1, function(row) sqrt(sum((trueXopt - row)^2)))
  final_overall_rmse <- mean(final_euclidean_distances)
  
  # Table C metrics
  z_j_current <- z_j[1:n_centers]
  
  stage2_expected_outcomes <- (beta_true_int1 * x_rec_stage2$x1_rec_stage2 + 
                                 beta_true_int2 * x_rec_stage2$x2_rec_stage2 +
                                 beta_z * mean(z_j_current))
  
  final_expected_outcomes <- (beta_true_int1 * x_final_est_opt$x1_final_est_opt + 
                                beta_true_int2 * x_final_est_opt$x2_final_est_opt +
                                beta_z * mean(z_j_current))
  
  # Confidence set calculations
  setCP_vec <- numeric(n_sims)
  setperc_vec <- numeric(n_sims)
  cbCP_vec <- numeric(n_sims)
  
  grid_size <- 0.1
  x_grid <- expand.grid(
    int1 = seq(x_min[1], x_max[1], by = grid_size),
    int2 = seq(x_min[2], x_max[2], by = grid_size)
  )
  n_grid_points <- nrow(x_grid)
  trueXopt_rounded <- c(round(trueXopt[1], 1), round(trueXopt[2], 1))
  
  for(sim in 1:n_sims) {
    tryCatch({
      model_combined <- lm(y ~ int1 + int2 + z_j - 1, data = result$combined_sim_results[[sim]]$data)
      mean_z_j_data <- mean(result$combined_sim_results[[sim]]$data$z_j)
      X_grid <- cbind(x_grid$int1, x_grid$int2, mean_z_j_data)
      vcov_sandwich <- sandwich(model_combined)
      se_sandwich <- sqrt(diag(X_grid %*% vcov_sandwich %*% t(X_grid)))
      pred_avg <- X_grid %*% coef(model_combined)
      ci_lower <- pred_avg - qnorm(0.975) * se_sandwich
      ci_upper <- pred_avg + qnorm(0.975) * se_sandwich
      
      aim_score <- beta_true_int1 * round(trueXopt[1], 1) + 
        beta_true_int2 * round(trueXopt[2], 1) + 
        beta_z * mean_z_j_data
      
      in_interval <- (ci_lower <= aim_score) & (ci_upper >= aim_score)
      setperc_vec[sim] <- sum(in_interval) / n_grid_points
      
      cs_indices <- which(in_interval)
      if(length(cs_indices) > 0) {
        cs <- x_grid[cs_indices, ]
        cs_rounded <- data.frame(int1 = round(cs[,1], 1), int2 = round(cs[,2], 1))
        in_cs <- sum((cs_rounded$int1 == trueXopt_rounded[1]) & (cs_rounded$int2 == trueXopt_rounded[2])) > 0
        setCP_vec[sim] <- as.numeric(in_cs)
      } else {
        setCP_vec[sim] <- 0
      }
      
      chi_critical <- sqrt(qchisq(0.95, df = 3))
      cb_lower <- pred_avg - chi_critical * se_sandwich
      cb_upper <- pred_avg + chi_critical * se_sandwich
      
      true_mean_outcome <- beta_true_int1 * x_grid$int1 + beta_true_int2 * x_grid$int2 + beta_z * mean_z_j_data
      all_in_bands <- all((cb_lower <= true_mean_outcome) & (cb_upper >= true_mean_outcome))
      cbCP_vec[sim] <- as.numeric(all_in_bands)
      
    }, error = function(e) {
      setperc_vec[sim] <- NA
      setCP_vec[sim] <- NA
      cbCP_vec[sim] <- NA
    })
  }
  
  setCP95 <- (sum(setCP_vec, na.rm = TRUE) / sum(!is.na(setCP_vec))) * 100
  setperc <- mean(setperc_vec, na.rm = TRUE) * 100
  bandsCP95 <- (sum(cbCP_vec, na.rm = TRUE) / sum(!is.na(cbCP_vec))) * 100
  
  return(list(
    table_a = data.frame(
      n_centers = n_centers,
      stage1_n = stage1_n,
      stage2_n = stage2_n,
      combined_n = combined_n_total,
      final_relbias_int1 = final_relbias_int1,
      final_se_ratio_int1 = final_se_ratio_int1,
      final_coverage_int1 = final_coverage_int1,
      final_relbias_int2 = final_relbias_int2,
      final_se_ratio_int2 = final_se_ratio_int2,
      final_coverage_int2 = final_coverage_int2
    ),
    table_b = data.frame(
      n_centers = n_centers,
      stage1_n = stage1_n,
      stage2_n = stage2_n,
      combined_n = combined_n_total,
      true_opt_int1 = trueXopt[1],
      true_opt_int2 = trueXopt[2],
      stage2_bias_int1 = stage2_bias_int1,
      stage2_bias_int2 = stage2_bias_int2,
      stage2_overall_rmse = stage2_overall_rmse,
      final_bias_int1 = final_bias_int1,
      final_bias_int2 = final_bias_int2,
      final_overall_rmse = final_overall_rmse
    ),
    table_c = data.frame(
      n_centers = n_centers,
      stage1_n = stage1_n,
      stage2_n = stage2_n,
      combined_n = combined_n_total,
      stage2_mean_outcome = mean(stage2_expected_outcomes),
      stage2_q025 = quantile(stage2_expected_outcomes, 0.025),
      stage2_q975 = quantile(stage2_expected_outcomes, 0.975),
      final_mean_outcome = mean(final_expected_outcomes),
      final_q025 = quantile(final_expected_outcomes, 0.025),
      final_q975 = quantile(final_expected_outcomes, 0.975),
      SetCP95 = setCP95,
      SetPerc = setperc,
      BandsCP95 = bandsCP95
    )
  ))
}

# Sequential processing - no cluster needed
metrics_results <- list()

for(i in 1:length(all_results)) {
  cat("Calculating metrics for scenario", i, "of", length(all_results), "...\n")
  
  metrics_results[[i]] <- calculate_scenario_metrics(
    i, all_results, trueXopt, beta_true_int1, beta_true_int2, 
    beta_z, z_j, x_min, x_max, n_sims
  )
}

cat("Metrics calculation completed.\n")

# Combine results into tables
Table1A_intervention_effects <- do.call(rbind, lapply(metrics_results, function(x) x$table_a))
Table1B_optimization_performance <- do.call(rbind, lapply(metrics_results, function(x) x$table_b))
Table1C_confidence_and_outcomes <- do.call(rbind, lapply(metrics_results, function(x) x$table_c))

cat("\n=== TABLE 1A: INTERVENTION EFFECT ESTIMATION PERFORMANCE ===\n")
print(Table1A_intervention_effects, digits = 3)

cat("\n=== TABLE 1B: OPTIMIZATION PERFORMANCE ===\n")
print(Table1B_optimization_performance, digits = 3)

cat("\n=== TABLE 1C: CONFIDENCE METRICS AND EXPECTED OUTCOMES ===\n")
print(Table1C_confidence_and_outcomes, digits = 3)

write.csv(Table1A_intervention_effects, "Table1A_sim_cubic_fixed_zj_highconf.csv", row.names = FALSE)
write.csv(Table1B_optimization_performance, "Table1B_sim_cubic_fixed_zj_highconf.csv", row.names = FALSE)
write.csv(Table1C_confidence_and_outcomes, "Table1C_sim_cubic_fixed_zj_highconf.csv", row.names = FALSE)

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")
cat("\nTotal simulation time:", round(total_time, 2), "minutes\n")
cat("Parallel processing used", n_cores, "cores for simulations\n")

cat("\n=== SIMULATION COMPLETE ===\n")
cat("All results saved to CSV files with '_cubic_fixed_zj' suffix\n")

# ==============================================================================
# VERIFY CORRELATIONS ACROSS ALL SIMULATIONS
# ==============================================================================
cat("\n=== VERIFYING CORRELATIONS ACROSS ALL SIMULATIONS ===\n")

all_cor_int1_z <- numeric(length(all_results) * n_sims)
all_cor_int2_z <- numeric(length(all_results) * n_sims)
counter <- 1

for(scenario_idx in 1:length(all_results)) {
  for(sim_idx in 1:n_sims) {
    combined_data <- all_results[[scenario_idx]]$combined_sim_results[[sim_idx]]$data
    
    # Get intervention group only
    intervention_data <- combined_data[combined_data$studyarm == 1, ]
    
    # Calculate correlations
    all_cor_int1_z[counter] <- cor(intervention_data$int1, intervention_data$z_j)
    all_cor_int2_z[counter] <- cor(intervention_data$int2, intervention_data$z_j)
    
    counter <- counter + 1
  }
}

cat("\nCorrelations Summary (Intervention Group, All Simulations):\n")
cat("Cor(int1, z_j): Mean =", round(mean(all_cor_int1_z), 4), 
    ", SD =", round(sd(all_cor_int1_z), 4), "(Target: 0.05)\n")
cat("Cor(int2, z_j): Mean =", round(mean(all_cor_int2_z), 4), 
    ", SD =", round(sd(all_cor_int2_z), 4), "(Target: 0.07)\n")