# Utility functions for EMMREML package

#' Check if a matrix is (effectively) an identity matrix
#'
#' @param K Matrix to check
#' @param tol Tolerance for comparison
#'
#' @return Logical indicating whether K is an identity matrix
#'
#' @keywords internal
is_identity_matrix <- function(K, tol = 1e-8) {
  # If Matrix package is available, use its efficient methods
  if (requireNamespace("Matrix", quietly = TRUE) && 
      inherits(K, "Matrix") && Matrix::isSymmetric(K)) {
    if (Matrix::isDiagonal(K)) {
      # Check if diagonal elements are 1
      return(all(abs(Matrix::diag(K) - 1) < tol))
    }
    return(FALSE)
  }
  
  # For large matrices, check a few key properties
  n <- nrow(K)
  
  # Check diagonal elements are 1
  diag_check <- all(abs(diag(K) - 1) < tol)
  
  # For small matrices, check all elements
  if (n < 500) {
    I <- diag(n)
    return(diag_check && all(abs(K - I) < tol))
  } 
  
  # For large matrices, sample some elements
  if (!diag_check) {
    return(FALSE)
  }
  
  # Sample some rows/columns
  sample_size <- min(100, n)
  idx <- sample(n, sample_size)
  I_sample <- diag(sample_size)
  K_sample <- K[idx, idx]
  
  return(all(abs(K_sample - I_sample) < tol))
}

#' Check if a matrix is (effectively) diagonal
#'
#' @param K Matrix to check
#' @param tol Tolerance for comparison
#'
#' @return Logical indicating whether K is a diagonal matrix
#'
#' @keywords internal
is_diagonal_matrix <- function(K, tol = 1e-8) {
  # Quick check using Matrix package if available
  if (requireNamespace("Matrix", quietly = TRUE)) {
    if (Matrix::isDiagonal(K)) {
      return(TRUE)
    }
  }
  
  n <- nrow(K)
  
  # For small matrices, check all elements
  if (n < 500) {
    is_diag <- TRUE
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j && abs(K[i, j]) > tol) {
          is_diag <- FALSE
          break
        }
      }
      if (!is_diag) break
    }
    return(is_diag)
  }
  
  # For large matrices, sample some elements
  sample_size <- min(100, n)
  idx <- sample(n, sample_size)
  K_sample <- K[idx, idx]
  
  is_diag <- TRUE
  for (i in 1:sample_size) {
    for (j in 1:sample_size) {
      if (i != j && abs(K_sample[i, j]) > tol) {
        is_diag <- FALSE
        break
      }
    }
    if (!is_diag) break
  }
  
  return(is_diag)
}
