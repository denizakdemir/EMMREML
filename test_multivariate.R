# EMMREML Multivariate Implementation Test 
cat("\n==== Testing EMMREML Multivariate Implementation ====\n")

library(EMMREML)
options(warn = 1)  # Show warnings as they occur

cat("\nCreating test data...\n")
set.seed(12345)

# Generate test data
n <- 20  # Number of subjects
p <- 2   # Number of traits
m <- 2   # Number of fixed effects

# Fixed effects design matrix with intercept
X <- cbind(1, rnorm(n))
colnames(X) <- c("intercept", "x")

# Simple kinship matrix
K <- matrix(0.5, n, n)  # Simple kinship with 0.5 covariance
diag(K) <- 1  # Diagonal is 1 (self-relationships)

# Create covariance matrices
Vg <- matrix(c(1.5, 0.5, 
               0.5, 1.0), p, p)  # Genetic covariance
Ve <- matrix(c(0.8, 0.2, 
               0.2, 0.6), p, p)  # Residual covariance

# Fixed effects coefficients
B <- matrix(c(10, 2,   # Effect 1 (intercept) for traits 1 & 2
              5, 1), 2, p)  # Effect 2 (x) for traits 1 & 2

# Create the phenotype matrix
# 1. Generate genetic effects using a Cholesky decomposition
Vg_chol <- chol(Vg)
G <- matrix(rnorm(n*p), n, p) %*% Vg_chol

# 2. Generate residual effects
Ve_chol <- chol(Ve)
E <- matrix(rnorm(n*p), n, p) %*% Ve_chol

# 3. Combine to get the phenotype matrix
Y <- X %*% B + G + E

# Create different scales for the two traits to test standardization
# We'll multiply trait 1 by 10 to create a scale difference
Y_scaled <- Y
Y_scaled[, 1] <- Y_scaled[, 1] * 10  
# Now trait 1 has much larger values and variance than trait 2

# Convert to data frame for formula interface
df <- data.frame(
  y1 = Y[,1],
  y2 = Y[,2],
  x = X[,2]
)

df_scaled <- data.frame(
  y1 = Y_scaled[,1],
  y2 = Y_scaled[,2],
  x = X[,2]
)

test_models <- function(test_name, data_df, should_standardize) {
  cat(paste0("\n===== ", test_name, " =====\n"))
  
  tryCatch({
    # Run the model
    result <- emmremlMultivariate(
      formula = cbind(y1, y2) ~ 1 + x,
      data = data_df,
      K = K,
      maxIter = 50,
      verbose = TRUE,
      standardize = should_standardize
    )
    
    # Show estimated parameters
    cat("\nFixed effects estimates:\n")
    print(result$Bhat)
    cat("\nTrue fixed effects (for comparison):\n")
    print(B)
    
    # Calculate R-squared as a simple goodness-of-fit measure
    Y_data <- as.matrix(data_df[, c("y1", "y2")])
    Y_hat <- t(result$Bhat) %*% t(X)
    SST <- sum((Y_data - mean(Y_data))^2)
    SSE <- sum((Y_data - t(Y_hat))^2)
    r_squared <- 1 - SSE/SST
    
    cat("\nR-squared:", r_squared, "\n")
    
    # Show convergence info
    cat("Iterations completed:", result$iterations, "\n")
    cat("Convergence achieved:", result$convergence, "\n")
    
    cat("\nVariance components estimates\n")
    cat("Estimated genetic covariance matrix (Vg):\n")
    print(result$Vg)
    cat("\nTrue genetic covariance matrix (original scale):\n")
    print(Vg)
    
    cat("\nEstimated residual covariance matrix (Ve):\n")
    print(result$Ve)
    cat("\nTrue residual covariance matrix (original scale):\n")
    print(Ve)
    
    # Display scaling info if available
    if (!is.null(result$scaling_info)) {
      cat("\nScaling information:\n")
      cat("Means:", result$scaling_info$means, "\n")
      cat("Standard deviations:", result$scaling_info$sds, "\n")
    }
    
    # Calculate relative errors for variance components
    # For scaled data, we need to adjust the true variances
    if (identical(data_df, df_scaled) && !should_standardize) {
      # If using scaled data without standardization, adjust true matrices for comparison
      scale_factors <- c(10, 1)
      scale_matrix <- diag(scale_factors)
      adjusted_Vg <- scale_matrix %*% Vg %*% scale_matrix
      adjusted_Ve <- scale_matrix %*% Ve %*% scale_matrix
      
      rel_error_Vg <- sum(abs(result$Vg - adjusted_Vg)) / sum(abs(adjusted_Vg))
      rel_error_Ve <- sum(abs(result$Ve - adjusted_Ve)) / sum(abs(adjusted_Ve))
      
      cat("\nAdjusted true genetic covariance matrix (scaled):\n")
      print(adjusted_Vg)
      cat("\nAdjusted true residual covariance matrix (scaled):\n")
      print(adjusted_Ve)
    } else {
      rel_error_Vg <- sum(abs(result$Vg - Vg)) / sum(abs(Vg))
      rel_error_Ve <- sum(abs(result$Ve - Ve)) / sum(abs(Ve))
    }
    
    cat("\nRelative error in Vg:", rel_error_Vg, "\n")
    cat("Relative error in Ve:", rel_error_Ve, "\n")
    
    # Determine overall test success
    if (result$convergence && rel_error_Vg < 0.5 && rel_error_Ve < 0.5) {
      cat("\nTest PASSED: Model successfully converged with reasonable variance estimates.\n")
    } else if (result$method == "OLS fallback") {
      cat("\nTest PARTIALLY PASSED: Used OLS fallback with basic variance estimates.\n")
    } else {
      cat("\nTest FAILED: Model did not converge or variance estimates are poor.\n")
    }
    
  }, error = function(e) {
    cat("Error in test:", conditionMessage(e), "\n")
    cat("Test FAILED due to error.\n")
  })
}

# Test 1: Original data without standardization
test_models("Original data without standardization", df, FALSE)

# Test 2: Original data with standardization
test_models("Original data with standardization", df, TRUE)

# Test 3: Scaled data without standardization
test_models("Scaled data without standardization", df_scaled, FALSE)

# Test 4: Scaled data with standardization
test_models("Scaled data with standardization", df_scaled, TRUE)

cat("\n==== Test Complete ====\n")