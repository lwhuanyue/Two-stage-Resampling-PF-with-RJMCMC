#'
#' Compute Acceptance Statistics from Transition Matrix
#'
#' This function analyzes a 2xN transition matrix `Mtx_Accep` and computes various counts 
#' of acceptance types and their combinations. Specifically designed for RJMCMC algorithms 
#' tracking multiple proposal types and their mutual acceptance.
#'
#' @param Mtx_Accep A 2 x N numeric matrix containing transition states. 
#'        Row 1: Previous state (values 1-4). Row 2: New state (values 1-4).
#'        Special value 0 indicates rejection, values 1-4 represent different acceptance types.
#'
#' @return A numeric vector of length 14 with the following elements:
#' \itemize{
#'   \item [1-4] Count of proposal types (a)-(d) 
#'   \item [5-8] Count of accepted proposal types (a)-(d)
#'   \item [9] Sum of accepted proposal
#'   \item [10-13] Acceptance rates for types (a)-(d)
#'   \item [14] Placeholder for future use
#' }
#'
#' @examples
#' \dontrun{
#' Mtx_Accep <- matrix(c(1,1, 2,1, 3,0, 4,1), nrow=2)
#' result <- compute_accept_stats(Mtx_Accep)
#' }
#'
#' @export
compute_accept_stats <- function(Mtx_Accep) {
  # Input validation
  if (!is.matrix(Mtx_Accep) || nrow(Mtx_Accep) != 2) {
    stop("Mtx_Accep must be a matrix with 2 rows")
  }
  
  # Initialize output vector (now length 14)
  temp_vector <- numeric(14)
  
  # Count individual acceptance types (elements 1-4)
  temp_vector[1] <- sum(Mtx_Accep[1,] == 1)
  temp_vector[2] <- sum(Mtx_Accep == 2) 
  temp_vector[3] <- sum(Mtx_Accep == 3)
  temp_vector[4] <- sum(Mtx_Accep == 4)
  
  # Count conditional transitions (elements 5-8)
  temp_vector[5] <- sum(Mtx_Accep[1,] == 1 & Mtx_Accep[2,] == 1)
  temp_vector[6] <- sum(Mtx_Accep[1,] == 2 & Mtx_Accep[2,] == 1)
  temp_vector[7] <- sum(Mtx_Accep[1,] == 3 & Mtx_Accep[2,] == 1)
  temp_vector[8] <- sum(Mtx_Accep[1,] == 4 & Mtx_Accep[2,] == 1)
  
  # Sum of conditional transitions (element 9)
  temp_vector[9] <- sum(temp_vector[5:8])
  
  # Calculate ratios: elements 10-13 = (5-8) / (1-4)
  for (i in 1:4) {
    if (temp_vector[i] != 0) {
      temp_vector[9 + i] <- temp_vector[4 + i] / temp_vector[i]
    } else {
      temp_vector[9 + i] <- 0  # Avoid division by zero
    }
  }
  
  # Element 14 remains as 0 (placeholder)
  
  return(temp_vector)
}
