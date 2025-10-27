# Define perturbation function
add_perturbation <- function(column) {
  # Copy original column data
  new_col <- column
  max_attempts <- 100  # Maximum number of attempts
  
  # Get the value of the first row
  first_row <- column[1]
  
  if (first_row == 0) {
    # Add perturbation A to the second row
    perturbation <- rnorm(1, 0, 0.01)
    perturbation <- sign(perturbation) * min(abs(perturbation), 0.05)
    new_col[2] <- new_col[2] + perturbation
    
    if (new_col[2] < 0.2) {new_col[2] <- 0.2}
    
  } else if (first_row == 1) {
    # Add perturbation B to the second row, perturbation A to the third row
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_attempts && !success) {
      temp_col <- new_col
      
      # Perturbation B (second row)
      perturbation_b <- rnorm(1, 0, 3)
      perturbation_b <- round(sign(perturbation_b) * min(abs(perturbation_b), 2))
      temp_col[2] <- temp_col[2] + perturbation_b
      
      # Perturbation A (rows 3,4)
      for (i in 3:4) {
        perturbation_a <- rnorm(1, 0, 0.01)
        perturbation_a <- sign(perturbation_a) * min(abs(perturbation_a), 0.05)
        temp_col[i] <- temp_col[i] + perturbation_a
        
        if (temp_col[i] < 0.2) {temp_col[i] <- 0.2}
      }
      
      # Validation: second row between 1-401, and maintaining ascending order
      if (temp_col[2] >= 1 && temp_col[2] <= 401 ) {
        new_col <- temp_col
        success <- TRUE
      }
      attempt <- attempt + 1
    }
    
    if (!success) {
      warning(paste("Failed to find valid perturbation for column after", max_attempts, "attempts"))
    }
    
  } else if (first_row == 2) {
    # Add perturbation B to the second row, perturbation A to rows 4,5,6
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_attempts && !success) {
      temp_col <- new_col
      
      # Perturbation B (second row)
      perturbation_b <- rnorm(1, 0, 3)
      perturbation_b <- round(sign(perturbation_b) * min(abs(perturbation_b), 2))
      temp_col[2] <- temp_col[2] + perturbation_b
      
      # Perturbation A (rows 4,5,6)
      for (i in 4:6) {
        perturbation_a <- rnorm(1, 0, 0.01)
        perturbation_a <- sign(perturbation_a) * min(abs(perturbation_a), 0.05)
        temp_col[i] <- temp_col[i] + perturbation_a
        
        if (temp_col[i] < 0.2) {temp_col[i] <- 0.2}
      }
      
      # Validation: second row between 1-401, and maintaining ascending order
      if (temp_col[2] >= 1 && temp_col[3] <= 401 && 
          temp_col[3] > temp_col[2]) {  # Check if first 6 rows are strictly increasing
        new_col <- temp_col
        success <- TRUE
      }
      attempt <- attempt + 1
    }
    
    if (!success) {
      warning(paste("Failed to find valid perturbation for column after", max_attempts, "attempts"))
    }
    
  } else if (first_row == 3) {
    # Add perturbation B to rows 2,3,4, perturbation A to rows 5,6,7,8
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_attempts && !success) {
      temp_col <- new_col
      
      # Perturbation B (rows 2,3,4)
      for (i in 2:4) {
        perturbation_b <- rnorm(1, 0, 3)
        perturbation_b <- round(sign(perturbation_b) * min(abs(perturbation_b), 2))
        temp_col[i] <- temp_col[i] + perturbation_b
      }
      
      # Perturbation A (rows 5,6,7,8)
      for (i in 5:8) {
        perturbation_a <- rnorm(1, 0, 0.01)
        perturbation_a <- sign(perturbation_a) * min(abs(perturbation_a), 0.05)
        temp_col[i] <- temp_col[i] + perturbation_a
        
        if (temp_col[i] < 0.2) {temp_col[i] <- 0.2}
      }
      
      # Validation: all B-perturbed rows (2,3,4) between 1-401, and maintaining ascending order
      b_rows_valid <- all(temp_col[2:4] >= 1 & temp_col[2:4] <= 402)
      sorted_valid <- temp_col[3] > temp_col[2] & temp_col[4] > temp_col[3]  # Check if entire column is strictly increasing
      
      if (b_rows_valid && sorted_valid) {
        new_col <- temp_col
        success <- TRUE
      }
      attempt <- attempt + 1
    }
    
    if (!success) {
      warning(paste("Failed to find valid perturbation for column after", max_attempts, "attempts"))
    }
  }
  
  return(new_col)
}
# Define perturbation function
add_perturbation <- function(column) {
  # Copy original column data
  new_col <- column
  max_attempts <- 100  # Maximum number of attempts
  
  # Get the value of the first row
  first_row <- column[1]
  
  if (first_row == 0) {
    # Add perturbation A to the second row
    perturbation <- rnorm(1, 0, 0.01)
    perturbation <- sign(perturbation) * min(abs(perturbation), 0.05)
    new_col[2] <- new_col[2] + perturbation
    
    if (new_col[2] < 0.2) {new_col[2] <- 0.2}
    
  } else if (first_row == 1) {
    # Add perturbation B to the second row, perturbation A to the third row
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_attempts && !success) {
      temp_col <- new_col
      
      # Perturbation B (second row)
      perturbation_b <- rnorm(1, 0, 3)
      perturbation_b <- round(sign(perturbation_b) * min(abs(perturbation_b), 2))
      temp_col[2] <- temp_col[2] + perturbation_b
      
      # Perturbation A (rows 3,4)
      for (i in 3:4) {
        perturbation_a <- rnorm(1, 0, 0.01)
        perturbation_a <- sign(perturbation_a) * min(abs(perturbation_a), 0.05)
        temp_col[i] <- temp_col[i] + perturbation_a
        
        if (temp_col[i] < 0.2) {temp_col[i] <- 0.2}
      }
      
      # Validation: second row between 1-401, and maintaining ascending order
      if (temp_col[2] >= 1 && temp_col[2] <= 401 ) {
        new_col <- temp_col
        success <- TRUE
      }
      attempt <- attempt + 1
    }
    
    if (!success) {
      warning(paste("Failed to find valid perturbation for column after", max_attempts, "attempts"))
    }
    
  } else if (first_row == 2) {
    # Add perturbation B to the second row, perturbation A to rows 4,5,6
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_attempts && !success) {
      temp_col <- new_col
      
      # Perturbation B (second row)
      perturbation_b <- rnorm(1, 0, 3)
      perturbation_b <- round(sign(perturbation_b) * min(abs(perturbation_b), 2))
      temp_col[2] <- temp_col[2] + perturbation_b
      
      # Perturbation A (rows 4,5,6)
      for (i in 4:6) {
        perturbation_a <- rnorm(1, 0, 0.01)
        perturbation_a <- sign(perturbation_a) * min(abs(perturbation_a), 0.05)
        temp_col[i] <- temp_col[i] + perturbation_a
        
        if (temp_col[i] < 0.2) {temp_col[i] <- 0.2}
      }
      
      # Validation: second row between 1-401, and maintaining ascending order
      if (temp_col[2] >= 1 && temp_col[3] <= 401 && 
          temp_col[3] > temp_col[2]) {  # Check if first 6 rows are strictly increasing
        new_col <- temp_col
        success <- TRUE
      }
      attempt <- attempt + 1
    }
    
    if (!success) {
      warning(paste("Failed to find valid perturbation for column after", max_attempts, "attempts"))
    }
    
  } else if (first_row == 3) {
    # Add perturbation B to rows 2,3,4, perturbation A to rows 5,6,7,8
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_attempts && !success) {
      temp_col <- new_col
      
      # Perturbation B (rows 2,3,4)
      for (i in 2:4) {
        perturbation_b <- rnorm(1, 0, 3)
        perturbation_b <- round(sign(perturbation_b) * min(abs(perturbation_b), 2))
        temp_col[i] <- temp_col[i] + perturbation_b
      }
      
      # Perturbation A (rows 5,6,7,8)
      for (i in 5:8) {
        perturbation_a <- rnorm(1, 0, 0.01)
        perturbation_a <- sign(perturbation_a) * min(abs(perturbation_a), 0.05)
        temp_col[i] <- temp_col[i] + perturbation_a
        
        if (temp_col[i] < 0.2) {temp_col[i] <- 0.2}
      }
      
      # Validation: all B-perturbed rows (2,3,4) between 1-401, and maintaining ascending order
      b_rows_valid <- all(temp_col[2:4] >= 1 & temp_col[2:4] <= 402)
      sorted_valid <- temp_col[3] > temp_col[2] & temp_col[4] > temp_col[3]  # Check if entire column is strictly increasing
      
      if (b_rows_valid && sorted_valid) {
        new_col <- temp_col
        success <- TRUE
      }
      attempt <- attempt + 1
    }
    
    if (!success) {
      warning(paste("Failed to find valid perturbation for column after", max_attempts, "attempts"))
    }
  }
  
  # Post-processing for model 3 (first_row == 3)
  if (first_row == 3) {
    # Step 1: Check for NaN values in position parameters (rows 2,3,4)
    position_params <- new_col[2:4]
    nan_indices <- which(is.nan(position_params))
    
    if (length(nan_indices) > 0) {
      # Calculate mean of non-NaN values and round to integer
      non_nan_values <- position_params[!is.nan(position_params)]
      if (length(non_nan_values) > 0) {
        replacement_value <- round(mean(non_nan_values))
        # Replace NaN values with the calculated mean
        new_col[1 + nan_indices] <- replacement_value  # +1 because position_params starts from row 2
      } else {
        # If all are NaN, use default value 200
        new_col[2:4] <- 200
      }
    }
    
    # Step 2: Check for values greater than 401 in position parameters
    position_params <- new_col[2:4]
    gt_401_indices <- which(position_params > 401)
    
    if (length(gt_401_indices) > 0) {
      # Subtract 20 from values greater than 401
      new_col[1 + gt_401_indices] <- position_params[gt_401_indices] - 20
    }
    
    # Step 3: Ensure position parameters (rows 2,3,4) are sorted in ascending order
    position_params <- new_col[2:4]
    if (!all(diff(position_params) > 0)) {
      # Sort the position parameters
      sorted_positions <- sort(position_params)
      new_col[2:4] <- sorted_positions
    }
    
    # Final validation check
    position_params <- new_col[2:4]
    if (any(is.nan(position_params))) {
      # If there are still NaN values after replacement, set to default
      nan_indices <- which(is.nan(position_params))
      new_col[1 + nan_indices] <- 200
    }
    
    if (any(position_params > 401)) {
      # If there are still values > 401, set to 401
      gt_401_indices <- which(position_params > 401)
      new_col[1 + gt_401_indices] <- 401
    }
    
    # Final sort to ensure ascending order
    position_params <- new_col[2:4]
    if (!all(diff(position_params) > 0)) {
      new_col[2:4] <- sort(position_params)
    }
  }
  
  return(new_col)
}

# # Apply perturbation to each column of EnPara matrix
# set.seed(123)  # Set random seed for reproducible results
# EnPara_perturbed <- EnPara
# 
# for (i in 1:ncol(EnPara)) {
#   EnPara_perturbed[, i] <- add_perturbation(EnPara[, i])
# }
# 
# # Display first few columns for comparison
# cat("Original EnPara (first 5 columns):\n")
# print(EnPara[, 1:5])
# 
# cat("\nPerturbed EnPara (first 5 columns):\n")
# print(EnPara_perturbed[, 1:5])
# 
# # Validation results
# cat("\nValidation results:\n")
# for (i in 1:ncol(EnPara_perturbed)) {
#   col_data <- EnPara_perturbed[, i]
#   first_row <- col_data[1]
#   
#   # Check if B-perturbed rows are between 1-401
#   if (first_row == 1) {
#     b_rows <- 2
#   } else if (first_row == 2) {
#     b_rows <- 2
#   } else if (first_row == 3) {
#     b_rows <- 2:4
#   } else {
#     b_rows <- NULL
#   }
#   
#   if (!is.null(b_rows)) {
#     b_valid <- all(col_data[b_rows] >= 1 & col_data[b_rows] <= 401)
#     sorted_valid <- all(diff(col_data) > 0)
#     cat(sprintf("Column %d: B-rows in [1,401] = %s, Sorted = %s\n", 
#                 i, b_valid, sorted_valid))
#   }
# }
# 
# # Perturbation statistics
# cat("\nPerturbation summary:\n")
# perturbation_diff <- EnPara_perturbed - EnPara
# perturbation_nonzero <- perturbation_diff[perturbation_diff != 0]
# 
# if (length(perturbation_nonzero) > 0) {
#   cat("Number of perturbed elements:", length(perturbation_nonzero), "\n")
#   cat("Range of perturbations:", range(perturbation_nonzero), "\n")
#   cat("Mean absolute perturbation:", mean(abs(perturbation_nonzero)), "\n")
# } else {
#   cat("No perturbations applied\n")
# }