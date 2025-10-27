#' Estimate probability density function for mixture models
#'
#' @param EnPara Input data matrix
#' @param bw_adjust Bandwidth adjustment parameter, default is 1
#' @return Returns a list containing all model information, and a function to calculate probability density values
estimate_mixture_density <- function(EnPara, bw_adjust = 1) {
  bw_adjust <- bw_adjust*10
  # Check and load required packages
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("MASS package is required. Please run: install.packages('MASS')")
  }
  
  # Extract model indices (first row)
  model_indices <- EnPara[1, ]
  
  # Initialize list to store data for each model
  model_data <- list(
    model0 = list(),
    model1 = list(),
    model2 = list(),
    model3 = list()
  )
  
  # Store previous KDE estimates (for handling cases with zero particles)
  previous_kdes <- list()
  
  # Separate data for different models
  for (i in 1:ncol(EnPara)) {
    model_idx <- model_indices[i]
    
    if (model_idx == 0) {
      # Model 0: Only velocity parameter (row 2)
      speed <- as.numeric(EnPara[2, i])
      if (!is.na(speed)) {
        model_data$model0$speed <- c(model_data$model0$speed, speed)
      }
    }
    else if (model_idx == 1) {
      # Model 1: Position parameter (row 2), velocity parameters (rows 3-4)
      position <- as.numeric(EnPara[2, i])
      speed1 <- as.numeric(EnPara[3, i])
      speed2 <- as.numeric(EnPara[4, i])
      
      if (!is.na(position) && !is.na(speed1) && !is.na(speed2)) {
        model_data$model1$position <- c(model_data$model1$position, position)
        model_data$model1$speed <- cbind(
          model_data$model1$speed, 
          c(speed1, speed2)
        )
      }
    }
    else if (model_idx == 2) {
      # Model 2: Position parameters (rows 2-3), velocity parameters (rows 4-6)
      position1 <- as.numeric(EnPara[2, i])
      position2 <- as.numeric(EnPara[3, i])
      speed1 <- as.numeric(EnPara[4, i])
      speed2 <- as.numeric(EnPara[5, i])
      speed3 <- as.numeric(EnPara[6, i])
      
      if (!is.na(position1) && !is.na(position2) && 
          !is.na(speed1) && !is.na(speed2) && !is.na(speed3)) {
        model_data$model2$position <- cbind(
          model_data$model2$position, 
          c(position1, position2)
        )
        model_data$model2$speed <- cbind(
          model_data$model2$speed, 
          c(speed1, speed2, speed3)
        )
      }
    }
    else if (model_idx == 3) {
      # Model 3: Position parameters (rows 2-4), velocity parameters (rows 5-8)
      position1 <- as.numeric(EnPara[2, i])
      position2 <- as.numeric(EnPara[3, i])
      position3 <- as.numeric(EnPara[4, i])
      speed1 <- as.numeric(EnPara[5, i])
      speed2 <- as.numeric(EnPara[6, i])
      speed3 <- as.numeric(EnPara[7, i])
      speed4 <- as.numeric(EnPara[8, i])
      
      if (!is.na(position1) && !is.na(position2) && !is.na(position3) &&
          !is.na(speed1) && !is.na(speed2) && !is.na(speed3) && !is.na(speed4)) {
        model_data$model3$position <- cbind(
          model_data$model3$position, 
          c(position1, position2, position3)
        )
        model_data$model3$speed <- cbind(
          model_data$model3$speed, 
          c(speed1, speed2, speed3, speed4)
        )
      }
    }
  }
  
  # Calculate particle count for each model
  n0 <- ifelse(is.null(model_data$model0$speed), 0, length(model_data$model0$speed))
  n1 <- ifelse(is.null(model_data$model1$position), 0, length(model_data$model1$position))
  n2 <- ifelse(is.null(model_data$model2$position), 0, ncol(model_data$model2$position))
  n3 <- ifelse(is.null(model_data$model3$position), 0, ncol(model_data$model3$position))
  
  # Apply special requirement 1: if particle count is less than 1, set to 1
  n0_adj <- max(n0, 1)
  n1_adj <- max(n1, 1)
  n2_adj <- max(n2, 1)
  n3_adj <- max(n3, 1)
  
  total_particles <- n0_adj + n1_adj + n2_adj + n3_adj
  
  # Calculate model probabilities (normalized)
  model_probs <- c(
    n0_adj / total_particles,
    n1_adj / total_particles,
    n2_adj / total_particles,
    n3_adj / total_particles
  )
  
  # Create KDE estimates for each model
  kde_models <- list()
  
  # Helper function: safe bandwidth calculation
  safe_bw <- function(x) {
    if (length(unique(x)) <= 1) {
      return(0.1)  # If data has little variation, return default bandwidth
    }
    return(bw.nrd0(x))
  }
  
  # Helper function: safe 2D KDE
  safe_kde2d <- function(x, y, h_adjust = 1, n = 50, lims = c(0, 1, 0, 1)) {
    valid_data <- complete.cases(x, y)
    x_clean <- x[valid_data]
    y_clean <- y[valid_data]
    
    if (length(x_clean) < 2) {
      # If insufficient data, return uniform distribution
      x_seq <- seq(lims[1], lims[2], length.out = n)
      y_seq <- seq(lims[3], lims[4], length.out = n)
      z <- matrix(1/((lims[2]-lims[1])*(lims[4]-lims[3])), nrow = n, ncol = n)
      return(list(x = x_seq, y = y_seq, z = z))
    }
    
    # Calculate bandwidth, ensure not zero
    h1 <- max(h_adjust * safe_bw(x_clean), 0.01)
    h2 <- max(h_adjust * safe_bw(y_clean), 0.01)
    
    tryCatch({
      return(MASS::kde2d(x_clean, y_clean, h = c(h1, h2), n = n, lims = lims))
    }, error = function(e) {
      # If kde2d fails, use larger bandwidth
      warning("kde2d failed, using larger bandwidth: ", e$message)
      h1 <- max(h1 * 2, 0.1)
      h2 <- max(h2 * 2, 0.1)
      return(MASS::kde2d(x_clean, y_clean, h = c(h1, h2), n = n, lims = lims))
    })
  }
  
  # Helper function: safe multivariate KDE (using independent 1D KDE as alternative)
  safe_multivariate_kde <- function(data, dim_names, bounds, bw_adjust = 1) {
    if (nrow(data) < 2) {
      # Insufficient data, return uniform distribution
      result <- list()
      for (i in 1:length(dim_names)) {
        result[[dim_names[i]]] <- list(
          x = seq(bounds[1], bounds[2], length.out = 100),
          y = rep(1/(bounds[2]-bounds[1]), 100)
        )
      }
      return(result)
    }
    
    # Use independent 1D KDE
    result <- list()
    for (i in 1:ncol(data)) {
      x <- data[, i]
      if (length(unique(x)) <= 1) {
        # If data has no variation, use uniform distribution
        result[[dim_names[i]]] <- list(
          x = seq(bounds[1], bounds[2], length.out = 100),
          y = rep(1/(bounds[2]-bounds[1]), 100)
        )
      } else {
        result[[dim_names[i]]] <- density(
          x, 
          bw = bw_adjust * safe_bw(x),
          from = bounds[1], to = bounds[2],
          na.rm = TRUE
        )
      }
    }
    return(result)
  }
  
  # Model 0 KDE
  if (n0 > 0) {
    kde_models$model0 <- density(
      model_data$model0$speed, 
      bw = bw_adjust * safe_bw(model_data$model0$speed),
      from = 0, to = 1,
      na.rm = TRUE
    )
    previous_kdes$model0 <- kde_models$model0
  } else if (!is.null(previous_kdes$model0)) {
    # Special requirement 2: if particle count is 0 but previous KDE exists, use previous one
    kde_models$model0 <- previous_kdes$model0
  }
  
  # Model 1 KDE (2D: position and velocity)
  if (n1 > 0) {
    # Position parameter KDE
    kde_models$model1_position <- density(
      model_data$model1$position,
      bw = bw_adjust * safe_bw(model_data$model1$position),
      from = 0, to = 400,
      na.rm = TRUE
    )
    
    # Velocity parameter KDE (2D) - using safe KDE function
    speed_data <- t(model_data$model1$speed)
    kde_models$model1_speed <- safe_kde2d(
      speed_data[, 1], speed_data[, 2],
      h_adjust = bw_adjust,
      lims = c(0, 1, 0, 1)
    )
    
    previous_kdes$model1_position <- kde_models$model1_position
    previous_kdes$model1_speed <- kde_models$model1_speed
  } else if (!is.null(previous_kdes$model1_position)) {
    kde_models$model1_position <- previous_kdes$model1_position
    kde_models$model1_speed <- previous_kdes$model1_speed
  }
  
  # Model 2 KDE
  if (n2 > 0) {
    # Position parameter KDE (2D) - using safe KDE function
    position_data <- t(model_data$model2$position)
    kde_models$model2_position <- safe_kde2d(
      position_data[, 1], position_data[, 2],
      h_adjust = bw_adjust,
      lims = c(0, 400, 0, 400)
    )
    
    # Velocity parameter KDE (3D) - using independent 1D KDE
    speed_data <- t(model_data$model2$speed)
    kde_models$model2_speed <- safe_multivariate_kde(
      speed_data,
      dim_names = paste0("speed", 1:3),
      bounds = c(0, 1),
      bw_adjust = bw_adjust
    )
    
    previous_kdes$model2_position <- kde_models$model2_position
    previous_kdes$model2_speed <- kde_models$model2_speed
  } else if (!is.null(previous_kdes$model2_position)) {
    kde_models$model2_position <- previous_kdes$model2_position
    kde_models$model2_speed <- previous_kdes$model2_speed
  }
  
  # Model 3 KDE
  if (n3 > 0) {
    # Position parameter KDE (3D) - using independent 1D KDE
    position_data <- t(model_data$model3$position)
    kde_models$model3_position <- safe_multivariate_kde(
      position_data,
      dim_names = paste0("pos", 1:3),
      bounds = c(0, 400),
      bw_adjust = bw_adjust
    )
    
    # Velocity parameter KDE (4D) - using independent 1D KDE
    speed_data <- t(model_data$model3$speed)
    kde_models$model3_speed <- safe_multivariate_kde(
      speed_data,
      dim_names = paste0("speed", 1:4),
      bounds = c(0, 1),
      bw_adjust = bw_adjust
    )
    
    previous_kdes$model3_position <- kde_models$model3_position
    previous_kdes$model3_speed <- kde_models$model3_speed
  } else if (!is.null(previous_kdes$model3_position)) {
    kde_models$model3_position <- previous_kdes$model3_position
    kde_models$model3_speed <- previous_kdes$model3_speed
  }
  
  # Print metrics
  cat("=== Mixture Model KDE Estimation Results ===\n")
  cat("Model probabilities:\n")
  cat(sprintf("Model 0: %.4f (Original particle count: %d)\n", model_probs[1], n0))
  cat(sprintf("Model 1: %.4f (Original particle count: %d)\n", model_probs[2], n1))
  cat(sprintf("Model 2: %.4f (Original particle count: %d)\n", model_probs[3], n2))
  cat(sprintf("Model 3: %.4f (Original particle count: %d)\n", model_probs[4], n3))
  cat(sprintf("Bandwidth adjustment parameter: %.2f\n", bw_adjust))
  cat(sprintf("Total particles (adjusted): %d\n", total_particles))
  
  # Create function to calculate probability density values
  calculate_density <- function(para_vector) {
    # Parameter check
    if (is.null(para_vector) || length(para_vector) < 2) {
      stop("Parameter vector must contain at least model index and one parameter")
    }
    
    # Extract model index (first element)
    model_idx <- para_vector[1]
    
    # Helper function: interpolation in KDE
    interp_1d <- function(kde, x) {
      if (is.null(kde) || is.null(kde$x) || is.null(kde$y)) return(1)  # Return 1 instead of 0 to avoid product being 0
      idx <- which.min(abs(kde$x - x))
      return(max(kde$y[idx], 1e-10))  # Avoid returning 0
    }
    
    # Helper function: 2D interpolation
    interp_2d <- function(kde2d, x, y) {
      if (is.null(kde2d) || is.null(kde2d$x) || is.null(kde2d$y) || is.null(kde2d$z)) return(1)
      x_idx <- which.min(abs(kde2d$x - x))
      y_idx <- which.min(abs(kde2d$y - y))
      return(max(kde2d$z[x_idx, y_idx], 1e-10))
    }
    
    if (model_idx == 0) {
      # Model 0: Only needs velocity parameter from row 2
      if (length(para_vector) < 2) {
        stop("Model 0 requires at least 2 parameters: [model index, velocity]")
      }
      speed <- para_vector[2]
      if (is.null(kde_models$model0)) {
        return(1e-10)  # Return very small value instead of 0
      }
      speed_dens <- interp_1d(kde_models$model0, speed)
      return(model_probs[1] * speed_dens)
    }
    else if (model_idx == 1) {
      # Model 1: Row 2 position parameter, rows 3-4 velocity parameters
      if (length(para_vector) < 4) {
        stop("Model 1 requires at least 4 parameters: [model index, position, velocity1, velocity2]")
      }
      position <- para_vector[2]
      speed1 <- para_vector[3]
      speed2 <- para_vector[4]
      
      if (is.null(kde_models$model1_position) || is.null(kde_models$model1_speed)) {
        return(1e-10)
      }
      # Position density
      pos_dens <- interp_1d(kde_models$model1_position, position)
      
      # Velocity density (2D interpolation)
      speed_dens <- interp_2d(kde_models$model1_speed, speed1, speed2)
      
      return(model_probs[2] * pos_dens * speed_dens)
    }
    else if (model_idx == 2) {
      # Model 2: Rows 2-3 position parameters, rows 4-6 velocity parameters
      if (length(para_vector) < 6) {
        stop("Model 2 requires at least 6 parameters: [model index, position1, position2, velocity1, velocity2, velocity3]")
      }
      position1 <- para_vector[2]
      position2 <- para_vector[3]
      speed1 <- para_vector[4]
      speed2 <- para_vector[5]
      speed3 <- para_vector[6]
      
      if (is.null(kde_models$model2_position) || is.null(kde_models$model2_speed)) {
        return(1e-10)
      }
      # Position density (2D interpolation)
      pos_dens <- interp_2d(kde_models$model2_position, position1, position2)
      
      # Velocity density - using independent 1D KDE
      speed_dens <- 1
      for (i in 1:3) {
        comp_name <- paste0("speed", i)
        speed_val <- para_vector[3 + i]  # speed1 at position 4, speed2 at position 5, speed3 at position 6
        speed_dens <- speed_dens * interp_1d(kde_models$model2_speed[[comp_name]], speed_val)
      }
      
      return(model_probs[3] * pos_dens * speed_dens)
    }
    else if (model_idx == 3) {
      # Model 3: Rows 2-4 position parameters, rows 5-8 velocity parameters
      if (length(para_vector) < 8) {
        stop("Model 3 requires at least 8 parameters: [model index, position1, position2, position3, velocity1, velocity2, velocity3, velocity4]")
      }
      position1 <- para_vector[2]
      position2 <- para_vector[3]
      position3 <- para_vector[4]
      speed1 <- para_vector[5]
      speed2 <- para_vector[6]
      speed3 <- para_vector[7]
      speed4 <- para_vector[8]
      
      if (is.null(kde_models$model3_position) || is.null(kde_models$model3_speed)) {
        return(1e-10)
      }
      # Position density - using independent 1D KDE
      pos_dens <- 1
      for (i in 1:3) {
        comp_name <- paste0("pos", i)
        pos_val <- para_vector[1 + i]  # position1 at position 2, position2 at position 3, position3 at position 4
        pos_dens <- pos_dens * interp_1d(kde_models$model3_position[[comp_name]], pos_val)
      }
      
      # Velocity density - using independent 1D KDE
      speed_dens <- 1
      for (i in 1:4) {
        comp_name <- paste0("speed", i)
        speed_val <- para_vector[4 + i]  # speed1 at position 5, speed2 at position 6, speed3 at position 7, speed4 at position 8
        speed_dens <- speed_dens * interp_1d(kde_models$model3_speed[[comp_name]], speed_val)
      }
      
      return(model_probs[4] * pos_dens * speed_dens)
    }
    else {
      stop("Invalid model index, must be one of 0,1,2,3")
    }
  }
  
  # Return results
  result <- list(
    model_probs = model_probs,
    kde_models = kde_models,
    model_data = model_data,
    calculate_density = calculate_density,
    bw_adjust = bw_adjust
  )
  
  class(result) <- "mixture_kde"
  return(result)
}

# Usage example
# result <- estimate_mixture_density(EnPara, bw_adjust = 1.0)
# 
# # Calculate probability density at a point
# # Model 0 example
# para_vector0 <- c(0, 0.5)
# dens0 <- result$calculate_density(para_vector0)
# cat(sprintf("Model 0 density at velocity=0.5: %.6f\n", dens0))
# 
# # Model 1 example  
# para_vector1 <- c(1, 100, 0.3, 0.4)
# dens1 <- result$calculate_density(para_vector1)
# cat(sprintf("Model 1 density at position=100, velocity=(0.3,0.4): %.6f\n", dens1))