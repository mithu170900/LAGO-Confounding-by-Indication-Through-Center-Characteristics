#---------------------------------
# Avg Counseling Minutes M1-M4
#---------------------------------
calculate_counseling_metrics_m1_m4 <- function(data) {
  results <- data.frame(
    center = character(),
    pid = character(),
    enroll_date_numeric = numeric(),
    total_minutes_m1_m4 = numeric(),
    total_active_days_m1_m4 = numeric(),
    avg_counseling_m1_m4 = numeric(),
    stringsAsFactors = FALSE
  )
  all_m4_pids <- unique(data$pid)
  
  for (pid in all_m4_pids) {
    
    patient_data <- data %>% filter(pid == !!pid)
    center <- unique(patient_data$center)
    window_start_patient_i <- unique(patient_data$window_start_m1_m4)
    window_end_patient_i <- unique(patient_data$window_end_m1_m4)
    
    facility_data <- data[data$center == center,]
    
    subset_patients_data <- facility_data %>% 
      filter(
        records_taken_date_counseling >= window_start_patient_i &
          records_taken_date_counseling <= window_end_patient_i &
          pc_date >= 0 &
          pc_date <= 119)
    
    total_minutes <- sum(subset_patients_data$timeDuration_pc_counsel_bp, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_meds, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_sideeffects, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_risk, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_diet, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_phys, na.rm = TRUE) +
      sum(subset_patients_data$timeDuration_pc_counsel_smoke, na.rm = TRUE)
    
    patient_stays <- subset_patients_data %>%
      group_by(pid) %>%
      slice(1) %>%
      mutate(length_of_stay = pmax(0, pmin(window_end_patient_i, window_end_m1_m4) - 
                                     pmax(window_start_patient_i, window_start_m1_m4) + 1)) %>%
      pull(length_of_stay)
    
    total_length_of_stay <- sum(patient_stays)
    
    avg_minutes <- total_minutes / total_length_of_stay * 10
    
    results <- rbind(results, data.frame(
      center = center,
      pid = pid,
      enroll_date_numeric = unique(patient_data$enroll_date_numeric),
      total_minutes_m1_m4 = total_minutes,
      total_active_days_m1_m4 = total_length_of_stay,
      avg_counseling_m1_m4 = avg_minutes,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}


#---------------------------------
# Avg Home BP Measurements M1-M4
#---------------------------------

calculate_homebp_m1_m4 <- function(data) {
  results <- data.frame(
    center = character(),
    pid = character(),
    enroll_date_numeric = numeric(),
    total_homebp_measurements_m1_m4 = numeric(),
    total_length_of_stay_m1_m4 = numeric(),
    avg_homebp_measurements_m1_m4 = numeric(),
    stringsAsFactors = FALSE
  )
  all_m4_pids <- unique(data$pid)
  
  for (pid in all_m4_pids) {
    
    # Get data from patient i
    patient_data <- data %>% filter(pid == !!pid)
    center <- unique(patient_data$center)
    window_start_patient_i <- unique(patient_data$window_start_m1_m4)
    window_end_patient_i <- unique(patient_data$window_end_m1_m4)
    
    facility_data <- data[data$center == center,]
    
    # Filter for homebp records in the M1-M4 window using exact matching
    subset_patients_data <- facility_data %>%
      filter(
          records_taken_date_homebp >= window_start_patient_i &
          records_taken_date_homebp <= window_end_patient_i &
          homebp_dateentry >= 0 & 
          homebp_dateentry <= 119
      )
    
    # Count homebp measurements
    total_homebp_measurements <- nrow(subset_patients_data[!is.na(subset_patients_data$homebp_sysbp), ])
    
    # Calculate length of stay for each patient
    patient_stays <- subset_patients_data %>%
      group_by(pid) %>%
      slice(1) %>%
      mutate(length_of_stay = pmax(0, pmin(window_end_patient_i, window_end_m1_m4) - 
                                     pmax(window_start_patient_i, window_start_m1_m4) + 1)) %>%
      pull(length_of_stay)
    
    # Sum them up
    total_length_of_stay <- sum(patient_stays)
    
    # Calculate the average homebp measurements per 10 days
    avg_homebp <- total_homebp_measurements / total_length_of_stay * 10
    
    results <- rbind(results, data.frame(
      center = center,
      pid = pid,
      enroll_date_numeric = unique(patient_data$enroll_date_numeric),
      total_homebp_measurements_m1_m4 = total_homebp_measurements,
      total_length_of_stay_m1_m4 = total_length_of_stay,
      avg_homebp_measurements_m1_m4 = avg_homebp,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}


