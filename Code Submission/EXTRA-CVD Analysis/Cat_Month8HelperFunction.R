#---------------------------------
# Avg Counseling Minutes M5-M8
#---------------------------------
calculate_counseling_metrics_m5_m8 <- function(data) {
  results <- data.frame(
    center = character(),
    pid = character(),
    enroll_date_numeric = numeric(),
    total_minutes_m5_m8 = numeric(),
    total_active_days_m5_m8 = numeric(),
    avg_counseling_m5_m8 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get all unique PIDs in the dataset
  all_m8_pids <- data %>% 
    filter(in_window_m5_m8 == 1) %>%
    pull(pid) %>%
    unique()
  
  for (pid in all_m8_pids) {
    # Get all data for this patient
    patient_data <- data %>% filter(pid == !!pid)
    center <- unique(patient_data$center)
    
    # Get patient's window dates
    window_start_patient_i <- unique(patient_data$window_start_m5_m8)
    window_end_patient_i <- unique(patient_data$window_end_m5_m8)
    
    facility_data <- data[data$center == center,]
    
    # Get all data from patient's center within their window
    subset_patients_data <- facility_data %>% 
      filter(
        records_taken_date_counseling >= window_start_patient_i &
          records_taken_date_counseling <= window_end_patient_i &
          in_window_m5_m8 == 1 &
          pc_date >= 120 &
          pc_date <= 239)
    
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
      mutate(length_of_stay = pmax(0, pmin(window_end_patient_i, window_end_m5_m8) - 
                                     pmax(window_start_patient_i, window_start_m5_m8) + 1)) %>%
      pull(length_of_stay)
    
    # Sum them up
    total_length_of_stay <- sum(patient_stays)
    
    avg_minutes <- total_minutes / total_length_of_stay * 10
    
    results <- rbind(results, data.frame(
      center = center,
      pid = pid,
      enroll_date_numeric = unique(patient_data$enroll_date_numeric),
      total_minutes_m5_m8 = total_minutes,
      total_active_days_m5_m8 = total_length_of_stay,
      avg_counseling_m5_m8 = avg_minutes,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}


#---------------------------------
# Avg Home BP Measurements M5-M8
#---------------------------------
calculate_homebp_m5_m8 <- function(data) {
  results <- data.frame(
    center = character(),
    pid = character(),
    enroll_date_numeric = numeric(),
    total_homebp_measurements_m5_m8 = numeric(),
    total_length_of_stay_m5_m8 = numeric(),
    avg_homebp_measurements_m5_m8 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get all unique PIDs in the dataset with valid Month 5-8 window
  # Only include patients who were still in the study during this period (in_window_m5_m8 == 1)
  all_m8_pids <- data %>% 
    filter(in_window_m5_m8 == 1) %>%
    pull(pid) %>%
    unique()
  
  # Process each patient individually
  for (pid in all_m8_pids) {
    # Get all data for this patient
    patient_data <- data %>% filter(pid == !!pid)
    center <- unique(patient_data$center)
    
    # Extract patient's specific window dates for months 5-8
    window_start_patient_i <- unique(patient_data$window_start_m5_m8)
    window_end_patient_i <- unique(patient_data$window_end_m5_m8)
    
    facility_data <- data[data$center == center,]
    
    # Filter for homebp records in the M5-M8 window using exact matching
    subset_patients_data <- facility_data %>%
      filter(
          records_taken_date_homebp >= window_start_patient_i &
          records_taken_date_homebp <= window_end_patient_i &
          in_window_m5_m8 == 1 &
          homebp_dateentry >= 120 &
          homebp_dateentry <= 239
      )
    
    # Count homebp measurements
    total_homebp_measurements <- nrow(subset_patients_data[!is.na(subset_patients_data$homebp_sysbp), ])
    
    # Calculate length of stay for each patient
    patient_stays <- subset_patients_data %>%
      group_by(pid) %>%
      slice(1) %>% # one record per patient
      mutate(length_of_stay = pmax(0, pmin(window_end_patient_i, window_end_m5_m8) - 
                                     pmax(window_start_patient_i, window_start_m5_m8) + 1)) %>%
      pull(length_of_stay)
    
    # Sum them up
    total_length_of_stay <- sum(patient_stays)
    
    # Calculate the average homebp measurements per month
    avg_homebp <- total_homebp_measurements / total_length_of_stay * 10
    
    # Add this patient's results to the output dataframe
    results <- rbind(results, data.frame(
      center = center,
      pid = pid,
      enroll_date_numeric = unique(patient_data$enroll_date_numeric),
      total_homebp_measurements_m5_m8 = total_homebp_measurements,
      total_length_of_stay_m5_m8 = total_length_of_stay,
      avg_homebp_measurements_m5_m8 = avg_homebp,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

