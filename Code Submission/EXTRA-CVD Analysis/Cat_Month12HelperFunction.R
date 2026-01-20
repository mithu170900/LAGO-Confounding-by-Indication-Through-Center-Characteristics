#---------------------------------
# Avg Counseling Minutes M9-M12
#---------------------------------
calculate_counseling_metrics_m9_m12 <- function(data) {
  results <- data.frame(
    center = character(),
    pid = character(),
    enroll_date_numeric = numeric(),
    total_minutes_m9_m12 = numeric(),
    total_active_days_m9_m12 = numeric(),
    avg_counseling_m9_m12 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get all unique PIDs in the dataset
  all_m12_pids <- data %>% 
    filter(in_window_m9_m12 == 1) %>%
    pull(pid) %>%
    unique()
  
  for (pid in all_m12_pids) {
    # Get all data for this patient
    patient_data <- data %>% filter(pid == !!pid)
    center <- unique(patient_data$center)
    
    # Get patient's window dates
    window_start_patient_i <- unique(patient_data$window_start_m9_m12)
    window_end_patient_i <- unique(patient_data$window_end_m9_m12)
    
    facility_data <- data[data$center == center,]
    
    # Get all data from patient's center within their window
    subset_patients_data <- facility_data %>% 
      filter(
        records_taken_date_counseling >= window_start_patient_i &
          records_taken_date_counseling <= window_end_patient_i &
          in_window_m9_m12 == 1 &
          pc_date >= 240 &
          pc_date <= 359)
    
    total_minutes <- sum(subset_patients_data$timeDuration_pc_counsel_bp, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_meds, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_sideeffects, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_risk, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_diet, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_phys, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_smoke, na.rm = TRUE)
    
    # Calculate length of stay for each patient
    patient_stays <- subset_patients_data %>%
      group_by(pid) %>%
      slice(1) %>%
      mutate(length_of_stay = pmax(0, pmin(window_end_patient_i, window_end_m9_m12) - 
                                     pmax(window_start_patient_i, window_start_m9_m12) + 1)) %>%
      pull(length_of_stay)
    
    # Sum them up
    total_length_of_stay <- sum(patient_stays)
    
    avg_minutes <- total_minutes / total_length_of_stay * 10
    
    results <- rbind(results, data.frame(
      center = center,
      pid = pid,
      enroll_date_numeric = unique(patient_data$enroll_date_numeric),
      total_minutes_m9_m12 = total_minutes,
      total_active_days_m9_m12 = total_length_of_stay,
      avg_counseling_m9_m12 = avg_minutes,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

#---------------------------------
# Avg Home BP Measurements M9-M12
#---------------------------------
calculate_homebp_m9_m12 <- function(data) {
  results <- data.frame(
    center = character(),
    pid = character(),
    enroll_date_numeric = numeric(),
    total_homebp_measurements_m9_m12 = numeric(),
    total_length_of_stay_m9_m12 = numeric(),
    avg_homebp_measurements_m9_m12 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get all unique PIDs in the dataset with valid M9-M12 window
  all_m12_pids <- data %>% 
    filter(in_window_m9_m12 == 1) %>%
    pull(pid) %>%
    unique()
  
  for (pid in all_m12_pids) {
    # Get all data for this patient
    patient_data <- data %>% filter(pid == !!pid)
    center <- unique(patient_data$center)
    
    # Get patient's window dates
    window_start_patient_i <- unique(patient_data$window_start_m9_m12)
    window_end_patient_i <- unique(patient_data$window_end_m9_m12)
    
    facility_data <- data[data$center == center,]
    # Filter for homebp records in the M9-M12 window using exact matching
    subset_patients_data <- facility_data %>%
      filter(
        records_taken_date_homebp >= window_start_patient_i &
          records_taken_date_homebp <= window_end_patient_i &
          in_window_m9_m12 == 1 &
          homebp_dateentry >= 240 &
          homebp_dateentry <= 359
      )
    
    
    # Count homebp measurements
    total_homebp_measurements <- nrow(subset_patients_data[!is.na(subset_patients_data$homebp_sysbp), ])
    
    # Calculate length of stay for each patient
    patient_stays <- subset_patients_data %>%
      group_by(pid) %>%
      slice(1) %>%
      mutate(length_of_stay = pmax(0, pmin(window_end_patient_i, window_end_m9_m12) - 
                                     pmax(window_start_patient_i, window_start_m9_m12) + 1)) %>%
      pull(length_of_stay)
    
    # Sum them up
    total_length_of_stay <- sum(patient_stays)
    
    avg_homebp <- total_homebp_measurements / total_length_of_stay * 10
    
    # Add to results
    results <- rbind(results, data.frame(
      center = center,
      pid = pid,
      enroll_date_numeric = unique(patient_data$enroll_date_numeric),
      total_homebp_measurements_m9_m12 = total_homebp_measurements,
      total_length_of_stay_m9_m12 = total_length_of_stay,
      avg_homebp_measurements_m9_m12 = avg_homebp,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

