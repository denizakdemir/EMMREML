# Simple test to debug the multivariate model

library(EMMREML)
options(warn = 1)  # Print warnings immediately
set.seed(123)  # For reproducibility

# Generate simple example data
cat("Generating multivariate example data...\n")
n <- 20  # sample size
p <- 2   # number of traits

# Create fixed effects
X <- cbind(1, rnorm(n))
colnames(X) <- c("intercept", "x")

# Create kinship matrix
K_raw <- matrix(rnorm(n*20), n, 20)
K <- tcrossprod(K_raw)
diag(K) <- diag(K) + 1e-6  # ensure positive definiteness
K <- K / mean(diag(K))  # Scale K

# Create Z matrix (identity)
Z <- diag(n)

# Create covariance matrices
Vg <- matrix(c(1.5, 0.5, 
               0.5, 1.0), p, p)
Ve <- matrix(c(0.8, 0.2, 
               0.2, 0.6), p, p)

# True fixed effects
B <- matrix(c(10, 2, 
              5, 1), 2, p)  
cat("Fixed effects matrix (B):\n")
print(B)

# Simulate data
# Generate the Cholesky decomposition for simulation
Vg_chol <- chol(Vg)
Ve_chol <- chol(Ve)

# Simulate multivariate phenotypes
G <- matrix(rnorm(n*p), n, p) %*% Vg_chol  # random genetic effects
E <- matrix(rnorm(n*p), n, p) %*% Ve_chol  # residual errors
Y <- X %*% B + G + E

# Create a version with different scales for traits
Y_scaled <- Y
Y_scaled[, 1] <- Y_scaled[, 1] * 10  # Trait 1 is now on a different scale

cat("Matrix dimensions check:\n")
cat("Y dimensions:", dim(Y)[1], "x", dim(Y)[2], "\n")
cat("X dimensions:", dim(X)[1], "x", dim(X)[2], "\n")
cat("Z dimensions:", dim(Z)[1], "x", dim(Z)[2], "\n")
cat("K dimensions:", dim(K)[1], "x", dim(K)[2], "\n")

# Convert to data frame for formula interface
df <- data.frame(Y)
colnames(df) <- c("y1", "y2")
df$x <- X[,2]

df_scaled <- data.frame(Y_scaled)
colnames(df_scaled) <- c("y1", "y2")
df_scaled$x <- X[,2]

# Create formula
formula_obj <- as.formula("cbind(y1, y2) ~ 1 + x")

run_test <- function(test_name, Y_data, should_standardize) {
  cat(paste0("\n===== ", test_name, " =====\n"))
  
  # Test with the core function directly
  tryCatch({
    result_core <- emmremlMultivariate_core(
      Y = Y_data, 
      X = X, 
      Z = Z, 
      K = K,
      maxIter = 20,
      standardize = should_standardize
    )
    
    cat("Direct call succeeded!\n")
    cat("Estimated Vg:\n")
    print(result_core$Vg)
    
    # For scaled data without standardization, adjust true matrices for comparison
    if (identical(Y_data, Y_scaled) && !should_standardize) {
      scale_factors <- c(10, 1)
      scale_matrix <- diag(scale_factors)
      adjusted_Vg <- scale_matrix %*% Vg %*% scale_matrix
      adjusted_Ve <- scale_matrix %*% Ve %*% scale_matrix
      
      cat("Adjusted True Vg (for scaled data):\n")
      print(adjusted_Vg)
      
      cat("Estimated Ve:\n")
      print(result_core$Ve)
      cat("Adjusted True Ve (for scaled data):\n")
      print(adjusted_Ve)
    } else {
      cat("True Vg:\n")
      print(Vg)
      
      cat("Estimated Ve:\n")
      print(result_core$Ve)
      cat("True Ve:\n")
      print(Ve)
    }
    
    cat("Estimated B:\n")
    print(result_core$Bhat)
    cat("True B (transposed for comparison):\n")
    print(t(B))
    
    # Display scaling info if available
    if (!is.null(result_core$scaling_info)) {
      cat("\nScaling information:\n")
      cat("Means:", result_core$scaling_info$means, "\n")
      cat("Standard deviations:", result_core$scaling_info$sds, "\n")
    }
    
  }, error = function(e) {
    cat("Error in direct core function call:", conditionMessage(e), "\n")
    cat("Traceback:\n")
    print(sys.calls())
  })
}

# Test with different configurations
run_test("Original data without standardization", Y, FALSE)
run_test("Original data with standardization", Y, TRUE)
run_test("Scaled data without standardization", Y_scaled, FALSE)
run_test("Scaled data with standardization", Y_scaled, TRUE)

# Test with formula interface
cat("\nTesting emmremlMultivariate formula interface...\n")
tryCatch({
  # With original data
  result_formula <- emmremlMultivariate(
    formula = formula_obj,
    data = df,
    K = K,
    maxIter = 20,
    verbose = TRUE,
    standardize = TRUE
  )
  
  cat("Formula interface call succeeded with original data!\n")
  cat("Estimated Vg:\n")
  print(result_formula$Vg)
  cat("Estimated Ve:\n")
  print(result_formula$Ve)
  cat("Estimated B:\n")
  print(result_formula$Bhat)
  
  # With scaled data
  result_formula_scaled <- emmremlMultivariate(
    formula = formula_obj,
    data = df_scaled,
    K = K,
    maxIter = 20,
    verbose = TRUE,
    standardize = TRUE
  )
  
  cat("\nFormula interface call succeeded with scaled data!\n")
  cat("Estimated Vg:\n")
  print(result_formula_scaled$Vg)
  cat("Estimated Ve:\n")
  print(result_formula_scaled$Ve)
  cat("Estimated B:\n")
  print(result_formula_scaled$Bhat)
  
  # Compare results
  cat("\nDifference in Vg between original and scaled data (should be small):\n")
  print(result_formula$Vg - result_formula_scaled$Vg)
  
}, error = function(e) {
  cat("Error in formula interface call:", conditionMessage(e), "\n")
  cat("Traceback:\n")
  print(sys.calls())
})

cat("\nTest complete.\n")