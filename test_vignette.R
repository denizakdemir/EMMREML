# Test script for EMMREML vignette code

library(EMMREML)

# Example Data
set.seed(123)
n <- 50  # Using a smaller sample size for quicker testing
message("Generating example data...")

# Create fixed effects
X <- cbind(1, rnorm(n))
colnames(X) <- c("intercept", "x")

# Create a realistic kinship matrix
Z <- diag(n)
K_raw <- matrix(rnorm(n*20), n, 20)
K <- tcrossprod(K_raw)
diag(K) <- diag(K) + 1e-6  # ensure positive definiteness
K <- K / mean(diag(K))  # Scale K

# Simulate data from the model
u <- crossprod(chol(K), rnorm(n)) * sqrt(2)  # random effect with Var=2*K
e <- rnorm(n)  # residual with Var=1
y <- 10 + 2 * X[, 2] + u + e

# Create data frame
df <- data.frame(y = y, x = X[, 2], id = 1:n)

# Fit model
message("\nFitting basic mixed model...")
tryCatch({
  # Using core implementation directly
  result <- emmreml_core(y = df$y, 
                       X = cbind(1, df$x),  # Intercept + x
                       Z = diag(nrow(df)),  # Identity for observations
                       K = K)
  
  message("Variance components:")
  message("Genetic variance (Vu):", result$Vu)
  message("Residual variance (Ve):", result$Ve)
  message("\nFixed effects:")
  names(result$betahat) <- c("Intercept", "x")
  print(result$betahat)
}, error = function(e) {
  message("Error in basic model:", e$message)
})

# Multi-Kernel Models
message("\nTesting multi-kernel model...")
tryCatch({
  # Create a second kinship matrix
  set.seed(456)
  K2_raw <- matrix(rnorm(n*15), n, 15)
  K2 <- tcrossprod(K2_raw)
  diag(K2) <- diag(K2) + 1e-6
  K2 <- K2 / mean(diag(K2))
  
  # Simulate new data with two random components
  u1 <- crossprod(chol(K), rnorm(n)) * sqrt(1.5)  # first random effect
  u2 <- crossprod(chol(K2), rnorm(n)) * sqrt(0.8)  # second random effect
  e <- rnorm(n) * sqrt(0.7)  # residual
  y_mk <- 10 + 2 * X[, 2] + u1 + u2 + e
  
  # Create data frame
  df_mk <- data.frame(y = y_mk, x = X[, 2], id = 1:n)
  
  # List of kernel matrices
  K_list <- list(K1 = K, K2 = K2)
  
  # Debugging information
  message("K matrices dimensions:")
  message("K1: ", nrow(K), "x", ncol(K))
  message("K2: ", nrow(K2), "x", ncol(K2))
  
  # Create Zlist for multi-kernel model
  Zlist <- list(Z1 = diag(nrow(df_mk)), Z2 = diag(nrow(df_mk)))
  
  # Fit multi-kernel model using main implementation
  X_mk <- cbind(1, df_mk$x)
  
  message("Testing regular emmremlMultiKernel function...")
  tryCatch({
    mk_result <- emmremlMultiKernel(y = df_mk$y,
                                  X = X_mk,
                                  Zlist = Zlist,
                                  Klist = K_list)
    
    message("Variance components:")
    print(mk_result$weights)
    message("\nFixed effects:")
    names(mk_result$betahat) <- c("Intercept", "x")
    print(mk_result$betahat)
  }, error = function(e) {
    message("Error in emmremlMultiKernel: ", e$message)
  })
}, error = function(e) {
  message("Error in multi-kernel model:", e$message)
})

message("\nTest completed!")
