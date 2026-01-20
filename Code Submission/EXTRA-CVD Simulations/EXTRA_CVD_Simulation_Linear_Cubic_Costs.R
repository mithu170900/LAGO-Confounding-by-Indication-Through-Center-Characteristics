# ==================================
# TWO-STAGE SIMULATION - PARALLEL VERSION
# ==================================
# This simulation compares LAGO optimization vs non-LAGO approach in a two-stage
# trial design after the EXTRA-CVD study.
# Stage 1: Use original intervention data, simulate outcomes and
#          estimate recommended interventions for Stage 2
# Stage 2:
#     For LAGO: Apply LAGO optimized recommended interventions after stage 1
#     For nonLAGO: Apply stage 1 recommended interventions.
# Notes:
# (x1_rec_stage1, x2_rec_stage1) are the stage 1 recommended interventions estimated using EXTRA-CVD data
# (x1_rec_stage2, x2_rec_stage2) are the stage 2 recommended interventions, obtained from the LAGO optimization at the end of stage 1 using simulated stage 1 data.
# (x1_final_est_opt, x2_final_est_opt) are the final estimated optimal interventions, obtained from the LAGO optimization by the end of stage 2 using simulated stage 1+2 data.
source("OptInterLinGeneral_extraCVD 2.R")

# ==============================================================================
# SETUP AND PARAMETERS
# ==============================================================================

start_time <- Sys.time()

# Load required libraries
library(dplyr)        # Data manipulation
library(sandwich)     # Robust standard errors
library(LAGO)         # Learning-Assisted Goal Optimization
library(lmtest)       # Coefficient testing
library(foreach)      # Parallel processing
library(doParallel)   # Parallel backend

# Setup parallel backend - use all available cores except one
n_cores <- detectCores() - 1  # Leave one core free for system processes
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat(sprintf("Using %d cores for parallel processing\n", n_cores))

# Load EXTRA-CVD month-12 data as Stage 1 baseline data.
extracvd_data <- read.csv('final_results_m1_m12.csv')
extracvd_data$center <- factor(extracvd_data$center)      # Centers A, B, C
extracvd_data$pid    <- factor(extracvd_data$pid)         # Patient ID
extracvd_data$studyarm <- as.integer(extracvd_data$studyarm)  # 0=control, 1=intervention

# Select only relevant columns for simulation
extracvd_data <- extracvd_data[, c('center', 'pid', 'enroll_date_numeric',
                                   'wd_date', 'studyarm',
                                   'avg_counseling_m1_m12',           # Intervention component 1
                                   'avg_homebp_measurements_m1_m12',  # Intervention component 2
                                   'delta_sbp_m1_m12')]               # Outcome: SBP change

# Fit outcome model to original EXTRA-CVD data to estimate parameters and calculate residuals
# Model: outcome = beta1*counseling + beta2*homebp + center_effects + error
model_extracvd <- lm(delta_sbp_m1_m12 ~ 0 + avg_counseling_m1_m12 + avg_homebp_measurements_m1_m12 + center, data = extracvd_data)
resid_extracvd <- residuals(model_extracvd) 

# True intervention effect parameters in simulation (estimated from real data - Table 3 in main text)
beta_true_int1 <- -1.59  # Effect of counseling sessions on month-12 SBP reduction (mmHg per session)
beta_true_int2 <- -0.59  # Effect of home BP measurements on month-12 SBP reduction (mmHg per measurement)
center_A <- -2.63        # Center A effect 
center_B <-  0.58        # Center B effect
center_C <-  2.11        # Center C effect
avg_center_effect <- mean(c(center_A, center_B, center_C))  # Equal weights across centers


# =================================================
# RESIDUALS and CENTER-SPECIFIC DEVIATIONS (ETA)
# =================================================
# Select intervention group data from EXTRA-CVD data.
extracvd_intervention_data <- extracvd_data[extracvd_data$studyarm == 1, ]

# Model intervention component 1 (counseling) based on center using EXTRA-CVD intervention data.
model_A1_extracvd <- lm(avg_counseling_m1_m12 ~ 0 + center, data = extracvd_intervention_data)

# Model intervention component 2 (home BP) based on center using EXTRA-CVD intervention data.
model_A2_extracvd <- lm(avg_homebp_measurements_m1_m12 ~ 0 + center, data = extracvd_intervention_data)

# Data for use later, including residuals from the intervention model for each participant.
# These residuals capture patient-level variability around their center's mean
paired_intervention_residuals <- data.frame(
  pid    = extracvd_intervention_data$pid,
  center = extracvd_intervention_data$center,
  xi1    = residuals(model_A1_extracvd) ,  # Residuals in counseling
  xi2    = residuals(model_A2_extracvd)  # Residuals in home BP
)

# Chose the stage 1 recommended interventions
x1_rec_stage1 <- mean(coef(model_A1_extracvd))  # Mean counseling across centers
x2_rec_stage1 <- mean(coef(model_A2_extracvd))  # Mean home BP across centers


# Center-specific deviations (eta) around recommended intervention for counseling
# Formula: actual_center_value = recommended_value + eta + xi
eta1_centerA <- coef(model_A1_extracvd)[1] - x1_rec_stage1
eta1_centerB <- coef(model_A1_extracvd)[2] - x1_rec_stage1
eta1_centerC <- coef(model_A1_extracvd)[3] - x1_rec_stage1

# Center-specific deviations (eta) around recommended intervention for home BP
eta2_centerA <- coef(model_A2_extracvd)[1] - x2_rec_stage1
eta2_centerB <- coef(model_A2_extracvd)[2] - x2_rec_stage1
eta2_centerC <- coef(model_A2_extracvd)[3] - x2_rec_stage1

# Expected outcome under the stage1 recommended interventions (x1_rec_stage1, x2_rec_stage1)
expected_out_rec_stage1_centerA <- beta_true_int1 * x1_rec_stage1 + beta_true_int2 * x2_rec_stage1 + center_A
expected_out_rec_stage1_centerB <- beta_true_int1 * x1_rec_stage1 + beta_true_int2 * x2_rec_stage1 + center_B
expected_out_rec_stage1_centerC <- beta_true_int1 * x1_rec_stage1 + beta_true_int2 * x2_rec_stage1 + center_C
expected_out_rec_stage1 <- (expected_out_rec_stage1_centerA + expected_out_rec_stage1_centerB + expected_out_rec_stage1_centerC) / 3

# Expected outcome under the stage1 actual interventions
# This can be calculated once before the simulation loop since Stage 1 uses fixed
# intervention assignments from EXTRA-CVD data

# Stage 1 intervention group - using original EXTRA-CVD data
extracvd_intervention_data_for_exp <- extracvd_data[extracvd_data$studyarm == 1, ]

# Calculate mean interventions per center for intervention group
stage1_intervention_means_by_center <- aggregate(
  cbind(avg_counseling_m1_m12, avg_homebp_measurements_m1_m12) ~ center,
  data = extracvd_intervention_data_for_exp,
  FUN = mean
)

# Calculate expected outcome under actual interventions for each center for stage 1
exp_act_centerA_stage1_constant <- beta_true_int1 * stage1_intervention_means_by_center[stage1_intervention_means_by_center$center == 'A', 'avg_counseling_m1_m12'] +
  beta_true_int2 * stage1_intervention_means_by_center[stage1_intervention_means_by_center$center == 'A', 'avg_homebp_measurements_m1_m12'] +
  center_A

exp_act_centerB_stage1_constant <- beta_true_int1 * stage1_intervention_means_by_center[stage1_intervention_means_by_center$center == 'B', 'avg_counseling_m1_m12'] +
  beta_true_int2 * stage1_intervention_means_by_center[stage1_intervention_means_by_center$center == 'B', 'avg_homebp_measurements_m1_m12'] +
  center_B

exp_act_centerC_stage1_constant <- beta_true_int1 * stage1_intervention_means_by_center[stage1_intervention_means_by_center$center == 'C', 'avg_counseling_m1_m12'] +
  beta_true_int2 * stage1_intervention_means_by_center[stage1_intervention_means_by_center$center == 'C', 'avg_homebp_measurements_m1_m12'] +
  center_C
expected_out_act_int_stage1 <- (exp_act_centerA_stage1_constant + exp_act_centerB_stage1_constant + exp_act_centerC_stage1_constant) / 3

# # ====================================================
# # CALCULATE TRUE OPTIMAL INTERVENTIONS FOR SIMULATIONS
# # ====================================================
# # Functions built by Daniel, used in both Daniel's and Ante's work
# # Flipping the signs since we are minimizing. Has been verified with the LAGO package.
# 
# my_beta <- c(beta_true_int1, beta_true_int2)  # (-1.59, -0.59)
# 
# trueXopt <- OptInterLinGeneral(
#   beta.vec      = -my_beta,    # Flipped the sign because we're minimizing 
#   gamma.vec     = 0,           
#   lin.cost.coef = c(1, 0.5),   # Linear cost function: cost = 1*counseling + 0.5*homebp
#   p.bar         = 7,           # Outcome goal
#   x.min         = c(0, 0),     # Minimum intervention values
#   x.max         = c(6.5, 3),   # Maximum intervention values
#   z             = avg_center_effect,  # Average center effect
#   intercept     = FALSE
# )$est.x.opt
# ====================================================
# CALCULATE TRUE OPTIMAL INTERVENTIONS FOR SIMULATIONS
# ====================================================
# Functions built by Daniel's, used in both Daniel's and Ante's
# Flipping the signs since we are minimizing. Has been verified with 
# the LAGO package.
cost_list <- list(
  c(0, 1.25, 0, -0.043, 0.0055),
  c(0, 0.63, 0, -0.09, 0.026)
)
calculate_cubic_true_optimum <- function(beta1, beta2, avg_center, x.min, x.max) {
  step_size <- 0.1
  x1_grid <- seq(x.min[1], x.max[1], by = step_size)
  x2_grid <- seq(x.min[2], x.max[2], by = step_size)
  grid <- expand.grid(x1 = x1_grid, x2 = x2_grid)
  
  # Calculate outcomes using the true model
  grid$outcome <- beta1 * grid$x1 + beta2 * grid$x2 + avg_center
  
  # Calculate CUBIC costs
  grid$cost <- (1.25 * grid$x1 - 0.043 * grid$x1^3 + 0.0055* grid$x1^4) +
    (0.63 * grid$x2 - 0.09 * grid$x2^3 + 0.026 * grid$x2^4)
  
  # Apply constraint: outcome <= -7, otherwise set cost to very high value
  grid$cost[grid$outcome > -9] <- 1e10
  
  # Find minimum cost solution
  best_idx <- which.min(grid$cost)
  
  # If no feasible solution found, return upper bounds
  if(grid$cost[best_idx] >= 1e10) {
    return(c(x.max[1], x.max[2]))
  }
  
  return(c(grid$x1[best_idx], grid$x2[best_idx]))
}

trueXopt <- calculate_cubic_true_optimum(
  beta_true_int1,
  beta_true_int2,
  mean(c(center_A, center_B, center_C)),           
  c(0, 0),
  c(6.5, 3)
)
# ======================================================
# SIMULATION FUNCTIONS
# ======================================================

# STAGE 1 SIMULATION FUNCTION.
# Generate new outcomes using original intervention package but new random errors
simulate_stage1 <- function(extracvd_baseline_data, seed_offset = 0) {
  set.seed(123 + seed_offset)
  
  sim_stage1_data <- extracvd_baseline_data  # Keep original intervention assignments from EXTRA-CVD
  
  # Re-sample outcome residuals (errors) from EXTRA-CVD data
  # This introduces variability in outcomes while preserving the error distribution
  errors <- sample(resid_extracvd, nrow(sim_stage1_data), replace = TRUE)
  
  # Apply center-level effects - each center has different baseline outcomes
  center_effects <- ifelse(sim_stage1_data$center == 'A', center_A,
                           ifelse(sim_stage1_data$center == 'B', center_B,
                                  ifelse(sim_stage1_data$center == 'C', center_C, NA)))
  
  # Generate new outcomes using true simulation model parameters
  # Outcome = beta1*counseling + beta2*homebp + center_effect + error
  sim_stage1_data$delta_sbp_m1_m12 <-
    beta_true_int1 * sim_stage1_data$avg_counseling_m1_m12 +
    beta_true_int2 * sim_stage1_data$avg_homebp_measurements_m1_m12 +
    center_effects +
    errors
  
  return(sim_stage1_data)
}

# STAGE 2 SIMULATION FUNCTION
simulate_stage2 <- function(recommended_x1, recommended_x2, sim_stage1_data, method = 'LAGO', seed_offset = 0) {
  set.seed(456 + seed_offset)  # Ensure reproducibility, different seed from stage 1.
  
  # Create Stage 2 data (3 centers) with 1:1 randomization within center
  sim_stage2_data <- data.frame(
    pid = 1:300,  # 300 participants (larger than Stage 1)
    center = factor(rep(c('A', 'B', 'C'), each = 100)),  # 100 per center
    studyarm = rep(c(rep(c(0, 1), each = 50)), 3)  # 50 intervention, 50 control in each center
  )
  
  # Initialize components at 0 (control group receives no intervention)
  sim_stage2_data$avg_counseling_m1_m12 <- 0
  sim_stage2_data$avg_homebp_measurements_m1_m12 <- 0
  
  # Identify intervention group - only they receive the intervention package
  intervention_indices <- sim_stage2_data$studyarm == 1
  
  n_intervention <- sum(intervention_indices)
  
  # Sample patient-level intervention residuals from EXTRA-CVD intervention patients
  sampled_patient_indices <- sample(nrow(paired_intervention_residuals), n_intervention, replace = TRUE)
  sampled_intervention_residuals <- paired_intervention_residuals[sampled_patient_indices, ]
  
  # Get the center assignment for each Stage 2 intervention patient
  intervention_centers <- sim_stage2_data$center[intervention_indices]
  
  # Apply center-specific effects for counseling
  eta1_for_stage2 <- ifelse(intervention_centers == 'A', eta1_centerA,
                            ifelse(intervention_centers == 'B', eta1_centerB,
                                   ifelse(intervention_centers == 'C', eta1_centerC, NA)))
  
  # Apply center-specific effects for home BP
  eta2_for_stage2 <- ifelse(intervention_centers == 'A', eta2_centerA,
                            ifelse(intervention_centers == 'B', eta2_centerB,
                                   ifelse(intervention_centers == 'C', eta2_centerC, NA)))
  
  # Calculate final intervention components
  sim_stage2_data$avg_counseling_m1_m12[intervention_indices] <-
    recommended_x1 + eta1_for_stage2 + sampled_intervention_residuals$xi1
  
  sim_stage2_data$avg_homebp_measurements_m1_m12[intervention_indices] <-
    recommended_x2 + eta2_for_stage2 + sampled_intervention_residuals$xi2
  
  # Apply bounds to keep interventions within feasible range
  sim_stage2_data$avg_counseling_m1_m12[intervention_indices] <-
    pmax(0, sim_stage2_data$avg_counseling_m1_m12[intervention_indices])
  
  sim_stage2_data$avg_homebp_measurements_m1_m12[intervention_indices] <-
    pmax(0, sim_stage2_data$avg_homebp_measurements_m1_m12[intervention_indices])
  
  # Generate Stage 2 outcomes using the same model as Stage 1
  center_effects <- ifelse(sim_stage2_data$center == 'A', center_A,
                           ifelse(sim_stage2_data$center == 'B', center_B,
                                  ifelse(sim_stage2_data$center == 'C', center_C, NA)))
  
  errors <- sample(resid_extracvd, nrow(sim_stage2_data), replace = TRUE)
  
  sim_stage2_data$delta_sbp_m1_m12 <-
    beta_true_int1 * sim_stage2_data$avg_counseling_m1_m12 +
    beta_true_int2 * sim_stage2_data$avg_homebp_measurements_m1_m12 +
    center_effects + errors
  
  return(sim_stage2_data)
}

# =================
# RUN SIMULATIONS
# =================

STAGE2_METHOD <- 'LAGO'  # 'LAGO' or 'nonLAGO'

n_sims  <- 2000  # Number of simulation iterations 

outcome_goal <- -9  # Target outcome: achieve at least -7 mmHg SBP reduction

cat(sprintf("Starting %d simulations with method: %s\n", n_sims, STAGE2_METHOD))

# Run parallel simulations
sim_results <- foreach(
  sim_id = 1:n_sims,
  .packages = c('dplyr', 'sandwich', 'LAGO', 'lmtest'),
  .export = c('simulate_stage1', 'simulate_stage2', 'extracvd_data',
              'resid_extracvd', 'beta_true_int1', 'beta_true_int2',
              'center_A', 'center_B', 'center_C', 'avg_center_effect',
              'paired_intervention_residuals', 'x1_rec_stage1', 'x2_rec_stage1',
              'eta1_centerA', 'eta1_centerB', 'eta1_centerC',
              'eta2_centerA', 'eta2_centerB', 'eta2_centerC',
              'STAGE2_METHOD', 'outcome_goal','cost_list')
) %dopar% {
  
  # -----------------------------------------------------------------------------------------
  # 1) Simulate Stage 1, then optimize on Stage 1 data to get recommendations FOR Stage 2
  # -----------------------------------------------------------------------------------------
  sim_stage1_data <- simulate_stage1(extracvd_data, seed_offset = sim_id)
  
  # Calculate average observed outcome in Stage 1 (intervention arm only)
  stage1_intervention_data <- sim_stage1_data[sim_stage1_data$studyarm == 1, ]
  avg_obs_stage1 <- mean(stage1_intervention_data$delta_sbp_m1_m12)
  
  # Optimize based on Stage 1 data to get recommendations for Stage 2
  opt_stage1 <- lago_optimization(
    data = sim_stage1_data,
    outcome_name = 'delta_sbp_m1_m12',
    outcome_type = 'continuous',
    glm_family = 'gaussian',
    link = 'identity',
    intervention_components = c('avg_counseling_m1_m12', 'avg_homebp_measurements_m1_m12'),
    prev_recommended_interventions = c(x1_rec_stage1, x2_rec_stage1), 
    #intervention_lower_bounds = c(0,0),  
    intervention_lower_bounds = c(x1_rec_stage1, x2_rec_stage1), 
    intervention_upper_bounds = c(6.5, 3.0),  
    include_center_effects = TRUE,  
    #cost_list_of_vectors = list(c(0, 1), c(0, 0.5)),
    cost_list_of_vectors = cost_list,
    outcome_goal = outcome_goal,
    center_weights_for_outcome_goal = c(1/3, 1/3, 1/3),  
    optimization_method = 'grid_search',
    optimization_grid_search_step_size = c(0.1, 0.1),  
    outcome_goal_intention = 'minimize',  
    include_confidence_set = FALSE
  )
  
  # Extract recommended intervention components and estimated outcome goal
  x1_rec_stage2 <- opt_stage1$rec_int[1]
  x2_rec_stage2 <- opt_stage1$rec_int[2]
  est_outcome_goal_stage1 <- opt_stage1$est_outcome_goal
  
  # ========== NEW: EXTRACT COST FROM STAGE 1 OPTIMIZATION ==========
  cost_rec_stage1 <- opt_stage1$rec_int_cost
  # =================================================================
  
  # Fit Stage 1 model to extract coefficients and 95% CI
  model_stage1 <- lm(delta_sbp_m1_m12 ~ 0 + avg_counseling_m1_m12
                     + avg_homebp_measurements_m1_m12 + center,
                     data = sim_stage1_data)
  
  stage1_coef <- coeftest(model_stage1, vcov = sandwich)  
  stage1_ci   <- coefci(model_stage1, vcov. = sandwich, level = 0.95)
  
  # --------------------------------------------------------------------------------
  # 2) Simulate Stage 2: Choose between LAGO or nonLAGO.
  # --------------------------------------------------------------------------------
  if (STAGE2_METHOD == 'LAGO') {
    sim_stage2_data <- simulate_stage2(x1_rec_stage2, x2_rec_stage2, sim_stage1_data, method = STAGE2_METHOD, seed_offset = sim_id)
  } else if (STAGE2_METHOD == 'nonLAGO') {
    sim_stage2_data <- simulate_stage2(x1_rec_stage1, x2_rec_stage1, sim_stage1_data, method = STAGE2_METHOD, seed_offset = sim_id)
  }
  
  # Calculate average observed outcome in Stage 2 (intervention arm only)
  stage2_intervention_data <- sim_stage2_data[sim_stage2_data$studyarm == 1, ]
  avg_obs_stage2 <- mean(stage2_intervention_data$delta_sbp_m1_m12)
  
  # ========== NEW: CALCULATE ACTUAL DELIVERED COST IN STAGE 2 ==========
  total_cost_stage2_actual <- mean(1.0 * stage2_intervention_data$avg_counseling_m1_m12 + 
                                     0.5 * stage2_intervention_data$avg_homebp_measurements_m1_m12)
  # =====================================================================
  
  # ----------------------------------------------------------------------------------------
  # 3) Combine data (Stage 1 + Stage 2) and optimize for final estimated optimal intervention package.
  # ----------------------------------------------------------------------------------------
  
  sim_stage1_data$stage <- 1
  sim_stage2_data$stage <- 2
  common_cols <- c('center', 'studyarm',
                   'avg_counseling_m1_m12',
                   'avg_homebp_measurements_m1_m12',
                   'delta_sbp_m1_m12', 'stage')
  
  combined_data <- rbind(sim_stage1_data[, common_cols],
                         sim_stage2_data[, common_cols])
  
  # Final optimization using complete Stage 1 + Stage 2 data
  opt_stage2 <- lago_optimization(
    data = combined_data,
    outcome_name = 'delta_sbp_m1_m12',
    outcome_type = 'continuous',
    glm_family = 'gaussian',
    link = 'identity',
    intervention_components = c('avg_counseling_m1_m12', 'avg_homebp_measurements_m1_m12'),
    prev_recommended_interventions = c(x1_rec_stage1, x2_rec_stage1), 
    #intervention_lower_bounds = c(0,0),  
    intervention_lower_bounds = c(x1_rec_stage1, x2_rec_stage1), 
    intervention_upper_bounds = c(6.5, 3.0),
    include_center_effects = TRUE,
    #cost_list_of_vectors = list(c(0, 1), c(0, 0.5)),
    cost_list_of_vectors = cost_list,
    outcome_goal = outcome_goal,
    center_weights_for_outcome_goal = c(1/3, 1/3, 1/3),
    optimization_method = 'grid_search',
    optimization_grid_search_step_size = c(0.1, 0.1),
    outcome_goal_intention = 'minimize',
    include_confidence_set = FALSE
  )
  
  # Extract final estimated optimal interventions
  x1_final_est_opt <- opt_stage2$rec_int[1]
  x2_final_est_opt <- opt_stage2$rec_int[2]
  est_outcome_goal_stage2 <- opt_stage2$est_outcome_goal
  
  # ========== NEW: EXTRACT COST FROM STAGE 2 OPTIMIZATION ==========
  cost_rec_stage2 <- opt_stage2$rec_int_cost
  # =================================================================
  
  # Calculate expected outcomes for Stage 2 intervention group
  sim_stage2_interventions_data <- sim_stage2_data[sim_stage2_data$studyarm == 1, ]
  
  sim_stage2_intervention_means_by_center <- aggregate(
    cbind(avg_counseling_m1_m12, avg_homebp_measurements_m1_m12) ~ center,
    data = sim_stage2_interventions_data,
    FUN = mean
  )
  
  # Calculate expected outcome for each center based on actual interventions delivered
  exp_act_centerA_stage2 <- beta_true_int1 * sim_stage2_intervention_means_by_center[sim_stage2_intervention_means_by_center$center == 'A', 'avg_counseling_m1_m12'] +
    beta_true_int2 * sim_stage2_intervention_means_by_center[sim_stage2_intervention_means_by_center$center == 'A', 'avg_homebp_measurements_m1_m12'] +
    center_A
  
  exp_act_centerB_stage2 <- beta_true_int1 * sim_stage2_intervention_means_by_center[sim_stage2_intervention_means_by_center$center == 'B', 'avg_counseling_m1_m12'] +
    beta_true_int2 * sim_stage2_intervention_means_by_center[sim_stage2_intervention_means_by_center$center == 'B', 'avg_homebp_measurements_m1_m12'] +
    center_B
  
  exp_act_centerC_stage2 <- beta_true_int1 * sim_stage2_intervention_means_by_center[sim_stage2_intervention_means_by_center$center == 'C', 'avg_counseling_m1_m12'] +
    beta_true_int2 * sim_stage2_intervention_means_by_center[sim_stage2_intervention_means_by_center$center == 'C', 'avg_homebp_measurements_m1_m12'] +
    center_C
  
  expected_out_act_int_stage2 <- (exp_act_centerA_stage2 + exp_act_centerB_stage2 + exp_act_centerC_stage2) / 3
  
  # Expected outcome under the recommendations used in Stage 2
  if (STAGE2_METHOD == 'LAGO') {
    expected_out_rec_stage2 <- beta_true_int1 * x1_rec_stage2 + beta_true_int2 * x2_rec_stage2 + avg_center_effect
  } else if (STAGE2_METHOD == 'nonLAGO') {
    expected_out_rec_stage2 <- beta_true_int1 * x1_rec_stage1 + beta_true_int2 * x2_rec_stage1 + avg_center_effect
  }
  
  # Expected outcome under final estimated optimal interventions
  expected_out_est_opt_int_stage2 <- beta_true_int1 * x1_final_est_opt + beta_true_int2 * x2_final_est_opt + avg_center_effect
  
  # Model fit on combined data for coefficient estimation metrics
  model_combined <- lm(delta_sbp_m1_m12 ~ 0 + avg_counseling_m1_m12
                       + avg_homebp_measurements_m1_m12 + center,
                       data = combined_data)
  
  combined_coef <- coeftest(model_combined, vcov = sandwich)
  combined_ci   <- coefci(model_combined, vcov. = sandwich, level = 0.95)
  
  # Return results for this simulation iteration
  list(
    x1_rec_stage2 = x1_rec_stage2,  
    x2_rec_stage2 = x2_rec_stage2,
    x1_final_est_opt = x1_final_est_opt,  
    x2_final_est_opt = x2_final_est_opt,
    expected_out_act_int_stage2 = expected_out_act_int_stage2,  
    expected_out_rec_stage2 = expected_out_rec_stage2,  
    expected_out_est_opt_int_stage2 = expected_out_est_opt_int_stage2,  
    est_outcome_goal_stage1 = est_outcome_goal_stage1,  
    est_outcome_goal_stage2 = est_outcome_goal_stage2, 
    stage1_coef = stage1_coef,  
    stage1_ci = stage1_ci,
    combined_coef = combined_coef,
    combined_ci = combined_ci,
    combined_data = combined_data,
    avg_obs_stage1 = avg_obs_stage1,
    avg_obs_stage2 = avg_obs_stage2,
    # ========== NEW: COST METRICS ==========
    cost_rec_stage1 = cost_rec_stage1,
    cost_rec_stage2 = cost_rec_stage2,
    total_cost_stage2_actual = total_cost_stage2_actual
    # =======================================
  )
}

# Stop the parallel cluster to free up resources
stopCluster(cl)

cat('\n\nSimulations completed. Processing results...\n')

# ==============================
# EXTRACT AND ANALYZE RESULTS
# ==============================

# Extract optimization results across all simulation iterations
x1_rec_stage2_vec <- sapply(sim_results, function(x) x$x1_rec_stage2)
x2_rec_stage2_vec <- sapply(sim_results, function(x) x$x2_rec_stage2)
x1_final_est_opt_vec <- sapply(sim_results, function(x) x$x1_final_est_opt)
x2_final_est_opt_vec <- sapply(sim_results, function(x) x$x2_final_est_opt)

# Extract estimated outcome goal values for feasibility assessment
est_outcome_goal_stage1_vec <- sapply(sim_results, function(x) x$est_outcome_goal_stage1)
est_outcome_goal_stage2_vec <- sapply(sim_results, function(x) x$est_outcome_goal_stage2)

# Extract average observed outcomes
avg_obs_stage1_vec <- sapply(sim_results, function(x) x$avg_obs_stage1)
avg_obs_stage2_vec <- sapply(sim_results, function(x) x$avg_obs_stage2)

# ========== NEW: EXTRACT COST VECTORS ==========
cost_rec_stage1_vec <- sapply(sim_results, function(x) x$cost_rec_stage1)
cost_rec_stage2_vec <- sapply(sim_results, function(x) x$cost_rec_stage2)
total_cost_stage2_actual_vec <- sapply(sim_results, function(x) x$total_cost_stage2_actual)
# ==============================================

# Calculate unfeasibility percentages
perc_unfeasible_stage1 <- mean(est_outcome_goal_stage1_vec > outcome_goal) * 100
perc_unfeasible_stage2 <- mean(est_outcome_goal_stage2_vec > outcome_goal) * 100

# Extract model coefficients and confidence intervals for all simulations
combined_sim_results <- lapply(sim_results, function(x) list(
  coef = x$combined_coef,
  ci   = x$combined_ci,
  data = x$combined_data
))

# ==========================================
# TABLE 1A: COEFFICIENT ESTIMATION METRICS
# ==========================================

# Extract Stage 1 coefficient estimates across all simulations
betahats_stage1 <- cbind(
  sapply(sim_results, function(x) x$stage1_coef['avg_counseling_m1_m12', 'Estimate']),
  sapply(sim_results, function(x) x$stage1_coef['avg_homebp_measurements_m1_m12', 'Estimate'])
)

# Extract Stage 1 standard errors
sehats_stage1 <- cbind(
  sapply(sim_results, function(x) x$stage1_coef['avg_counseling_m1_m12', 'Std. Error']),
  sapply(sim_results, function(x) x$stage1_coef['avg_homebp_measurements_m1_m12', 'Std. Error'])
)

# Extract Stage 2 (combined) coefficient estimates across all simulations
betahats_stage2 <- cbind(
  sapply(combined_sim_results, function(x) x$coef['avg_counseling_m1_m12', 'Estimate']),
  sapply(combined_sim_results, function(x) x$coef['avg_homebp_measurements_m1_m12', 'Estimate'])
)

# Extract Stage 2 standard errors
sehats_stage2 <- cbind(
  sapply(combined_sim_results, function(x) x$coef['avg_counseling_m1_m12', 'Std. Error']),
  sapply(combined_sim_results, function(x) x$coef['avg_homebp_measurements_m1_m12', 'Std. Error'])
)

# Extract confidence intervals for coverage calculations
beta1ci <- cbind(
  sapply(combined_sim_results, function(x) x$ci['avg_counseling_m1_m12', 1]),
  sapply(combined_sim_results, function(x) x$ci['avg_counseling_m1_m12', 2])
)
beta2ci <- cbind(
  sapply(combined_sim_results, function(x) x$ci['avg_homebp_measurements_m1_m12', 1]),
  sapply(combined_sim_results, function(x) x$ci['avg_homebp_measurements_m1_m12', 2])
)

# Create Table 1A with key estimation metrics
table1A <- data.frame(
  Method = STAGE2_METHOD,
  Parameter = c('Beta1 (Counseling)', 'Beta2 (HomeBP)'),
  Relative_Bias = c(
    mean(100 * (betahats_stage2[,1] - beta_true_int1) / abs(beta_true_int1)),
    mean(100 * (betahats_stage2[,2] - beta_true_int2) / abs(beta_true_int2))
  ),
  SE_EmpSD_Ratio = c(
    mean(sehats_stage2[,1]) / sd(betahats_stage2[,1]) * 100,
    mean(sehats_stage2[,2]) / sd(betahats_stage2[,2]) * 100
  ),
  Coverage = c(
    mean((beta1ci[,1] <= beta_true_int1) & (beta1ci[,2] >= beta_true_int1)) * 100,
    mean((beta2ci[,1] <= beta_true_int2) & (beta2ci[,2] >= beta_true_int2)) * 100
  )
)

# ==========================================
# TABLE 1B: OPTIMIZATION PERFORMANCE
# ==========================================

table1B <- data.frame(
  Method = STAGE2_METHOD,
  Stage = c(1, 2),
  Bias_Int1 = c(
    mean(x1_rec_stage2_vec) - trueXopt[1],
    mean(x1_final_est_opt_vec) - trueXopt[1]
  ),
  Bias_Int2 = c(
    mean(x2_rec_stage2_vec) - trueXopt[2],
    mean(x2_final_est_opt_vec) - trueXopt[2]
  ),
  RMSE = c(
    sqrt(mean((x1_rec_stage2_vec - trueXopt[1])^2 + (x2_rec_stage2_vec - trueXopt[2])^2)),
    sqrt(mean((x1_final_est_opt_vec - trueXopt[1])^2 + (x2_final_est_opt_vec - trueXopt[2])^2))
  )
)

# =====================================================
# TABLE 1C: CONFIDENCE METRICS AND EXPECTED OUTCOMES
# =====================================================

# Extract expected outcomes under recommended and final estimated optimal interventions
expected_out_rec_int_stage2_vec <- sapply(sim_results, function(x) x$expected_out_rec_stage2)
expected_out_est_opt_int_stage2_vec <- sapply(sim_results, function(x) x$expected_out_est_opt_int_stage2)

# Calculate percentage of simulations meeting the clinical goal
perc_meet_goal_stage1 <- mean(expected_out_rec_int_stage2_vec <= outcome_goal) * 100
perc_meet_goal_stage2 <- mean(expected_out_est_opt_int_stage2_vec <= outcome_goal) * 100

# Confidence set and simultaneous confidence bands calculation
setCP_vec  <- numeric(n_sims)  
setperc_vec<- numeric(n_sims)  
cbCP_vec   <- numeric(n_sims)  

# Create grid over intervention space for confidence set calculations
grid_size <- 0.1
x_grid <- expand.grid(
  int1 = seq(0, 6.5, by = grid_size),
  int2 = seq(0, 3,   by = grid_size)
)

n_grid_points <- nrow(x_grid)
trueXopt_rounded <- c(round(trueXopt[1], 1), round(trueXopt[2], 1))

# Loop over simulations to calculate confidence metrics
for (sim in 1:n_sims) {
  tryCatch({
    model_combined <- lm(delta_sbp_m1_m12 ~ 0 + avg_counseling_m1_m12
                         + avg_homebp_measurements_m1_m12 + center,
                         data = combined_sim_results[[sim]]$data)
    
    all_coefs <- coef(model_combined)
    
    all_param_names <- c('avg_counseling_m1_m12',
                         'avg_homebp_measurements_m1_m12',
                         'centerA',
                         'centerB',
                         'centerC')
    
    intervention_coefs <- all_coefs[all_param_names]
    
    x_grid_full <- cbind(
      x_grid$int1,
      x_grid$int2,
      1/3,
      1/3,
      1/3
    )
    
    predicted_average_outcome <- x_grid_full %*% intervention_coefs
    
    vcov_full <- sandwich(model_combined)
    vcov_full <- vcov_full[all_param_names, all_param_names]
    
    se_full <- sqrt(diag(x_grid_full %*% vcov_full %*% t(x_grid_full)))
    
    ci_lower <- predicted_average_outcome - qnorm(0.975) * se_full
    ci_upper <- predicted_average_outcome + qnorm(0.975) * se_full
    
    in_interval <- (ci_lower <= outcome_goal) & (ci_upper >= outcome_goal)
    setperc_vec[sim] <- sum(in_interval) / n_grid_points
    
    cs_indices <- which(in_interval)
    if (length(cs_indices) > 0) {
      cs <- x_grid[cs_indices, ]
      cs_rounded <- data.frame(int1 = round(cs[,1], 1), int2 = round(cs[,2], 1))
      
      in_cs <- any((cs_rounded$int1 == trueXopt_rounded[1]) &
                     (cs_rounded$int2 == trueXopt_rounded[2]))
      
      setCP_vec[sim] <- as.numeric(in_cs)
    } else {
      setCP_vec[sim] <- 0
    }
    
    chi_critical <- sqrt(qchisq(0.95, df = 5))
    cb_lower <- predicted_average_outcome - chi_critical * se_full
    cb_upper <- predicted_average_outcome + chi_critical * se_full
    
    true_mean_outcome <- beta_true_int1 * x_grid$int1 +
      beta_true_int2 * x_grid$int2 +
      avg_center_effect
    
    cbCP_vec[sim] <- as.numeric(all((cb_lower <= true_mean_outcome) &
                                      (cb_upper >= true_mean_outcome)))
  }, error = function(e) {
    setperc_vec[sim] <- NA
    setCP_vec[sim]   <- NA
    cbCP_vec[sim]    <- NA
  })
}

# Calculate summary statistics for confidence metrics
setCP95 <- (sum(setCP_vec) / sum(!is.na(setCP_vec))) * 100
setperc <- mean(setperc_vec, na.rm = TRUE) * 100
bandsCP95 <- (sum(cbCP_vec) / sum(!is.na(cbCP_vec))) * 100

# Create Table 1C with clinical and confidence metrics
table1C <- data.frame(
  Method = STAGE2_METHOD,
  Stage = c(1, 2),
  ExpectedOutActInt = c(
    expected_out_act_int_stage1,
    mean(sapply(sim_results, function(x) x$expected_out_act_int_stage2))
  ),
  MedianExpectedOutActInt = c(
    expected_out_act_int_stage1,
    median(sapply(sim_results, function(x) x$expected_out_act_int_stage2))
  ),
  ExpectedOutRecInt = c(
    expected_out_rec_stage1,
    mean(expected_out_rec_int_stage2_vec)
  ),
  MedianExpectedOutRecInt = c(
    expected_out_rec_stage1,
    median(expected_out_rec_int_stage2_vec)
  ),
  ExpectedOutEstOptInt = c(NA, mean(expected_out_est_opt_int_stage2_vec)),
  MedianExpectedOutEstOptInt = c(NA, median(expected_out_est_opt_int_stage2_vec)),
  AvgObsOut = c(
    mean(avg_obs_stage1_vec),
    mean(avg_obs_stage2_vec)
  ),
  MedianObsOut = c(
    median(avg_obs_stage1_vec),
    median(avg_obs_stage2_vec)
  ),
  # ========== NEW: COST METRICS ==========
  MeanCostRec = c(
    mean(cost_rec_stage1_vec),
    mean(cost_rec_stage2_vec)
  ),
  MedianCostRec = c(
    median(cost_rec_stage1_vec),
    median(cost_rec_stage2_vec)
  ),
  MeanCostActual = c(
    NA,
    mean(total_cost_stage2_actual_vec)
  ),
  MedianCostActual = c(
    NA,
    median(total_cost_stage2_actual_vec)
  ),
  # =======================================
  PercMeetGoal = c(perc_meet_goal_stage1, perc_meet_goal_stage2),
  PercUnfeasible = c(perc_unfeasible_stage1, perc_unfeasible_stage2),
  Q025 = c(quantile(expected_out_rec_int_stage2_vec, 0.025),
           quantile(expected_out_est_opt_int_stage2_vec, 0.025)),
  Q975 = c(quantile(expected_out_rec_int_stage2_vec, 0.975),
           quantile(expected_out_est_opt_int_stage2_vec, 0.975))
)

# Add confidence metrics (only available for Stage 2)
table1C$SetCP95   <- c(NA, setCP95)
table1C$SetPerc   <- c(NA, setperc)
table1C$BandsCP95 <- c(NA, bandsCP95)

rownames(table1C) <- NULL

# ==============================================================================
# PRINT RESULTS
# ==============================================================================

cat(sprintf('\n=== TABLE 1A: COEFFICIENT METRICS (GOAL %d) — %s ===\n', outcome_goal, STAGE2_METHOD))
print(table1A)

cat(sprintf('\n=== TABLE 1B: OPTIMIZATION PERFORMANCE (GOAL %d) — %s ===\n', outcome_goal, STAGE2_METHOD))
print(table1B)

cat(sprintf('\n=== TABLE 1C: EXPECTED OUTCOMES AND CONFIDENCE METRICS (GOAL %d) — %s ===\n', outcome_goal, STAGE2_METHOD))
print(table1C)

cat("\n=== OBSERVED VS EXPECTED OUTCOMES ===\n")
cat(sprintf("Stage 1 - Mean Observed: %.2f, Median Observed: %.2f, Expected (Actual Int): %.2f\n", 
            mean(avg_obs_stage1_vec), median(avg_obs_stage1_vec), expected_out_act_int_stage1))
cat(sprintf("Stage 2 - Mean Observed: %.2f, Median Observed: %.2f, Expected (Actual Int): %.2f, Median Expected: %.2f\n", 
            mean(avg_obs_stage2_vec), median(avg_obs_stage2_vec), 
            mean(sapply(sim_results, function(x) x$expected_out_act_int_stage2)),
            median(sapply(sim_results, function(x) x$expected_out_act_int_stage2))))

# ========== NEW: COST SUMMARY ==========
cat("\n=== INTERVENTION COSTS ===\n")
cat(sprintf("Stage 1 Recommended Cost - Mean: %.2f, Median: %.2f\n", 
            mean(cost_rec_stage1_vec), median(cost_rec_stage1_vec)))
cat(sprintf("Stage 2 Recommended Cost - Mean: %.2f, Median: %.2f\n", 
            mean(cost_rec_stage2_vec), median(cost_rec_stage2_vec)))
cat(sprintf("Stage 2 Actual Delivered Cost - Mean: %.2f, Median: %.2f\n", 
            mean(total_cost_stage2_actual_vec), median(total_cost_stage2_actual_vec)))
# =======================================

# ==============================================================================
# SAVE RESULTS
# ============================================================================== 

write.csv(table1A, paste0('Table1A_', STAGE2_METHOD, '_goal', abs(outcome_goal), '_cubic_lowerbound.csv'), row.names = FALSE)
write.csv(table1B, paste0('Table1B_', STAGE2_METHOD, '_goal', abs(outcome_goal), '_cubic_lowerbound.csv'), row.names = FALSE)
write.csv(table1C, paste0('Table1C_', STAGE2_METHOD, '_goal', abs(outcome_goal), '_cubic_lowerbound.csv'), row.names = FALSE)

# Report total computation time
end_time <- Sys.time()
cat('\nTotal time:', round(difftime(end_time, start_time, units = 'mins'), 2), 'minutes\n')