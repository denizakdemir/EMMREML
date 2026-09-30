# Test that standardization works correctly in the multivariate REML model

library(EMMREML)
options(warn = 1)  # Show warnings as they occur
set.seed(12345)

# Generate test data
n <- 50  # Number of subjects
p <- 2   # Number of traits
m <- 2   # Number of fixed effects

# Fixed effects design matrix with intercept
X <- cbind(1, rnorm(n))
colnames(X) <- c("intercept", "x")

# Create a kinship matrix
K_raw <- matrix(rnorm(n*30), n, 30)
K <- tcrossprod(K_raw)
K <- K / mean(diag(K))  # Scale K so diagonal has mean 1
diag(K) <- diag(K) + 1e-6  # Ensure positive definiteness

# Create covariance matrices
Vg <- matrix(c(1.5, 0.5, 
               0.5, 1.0), p, p)  # Genetic covariance
Ve <- matrix(c(0.8, 0.2, 
               0.2, 0.6), p, p)  # Residual covariance

# Fixed effects coefficients
B <- matrix(c(10, 2,   # Effect 1 (intercept) for traits 1 & 2
              5, 1), 2, p)  # Effect 2 (x) for traits 1 & 2

# Create the phenotype matrix
Vg_chol <- chol(Vg)
G <- matrix(rnorm(n*p), n, p) %*% Vg_chol
Ve_chol <- chol(Ve)
E <- matrix(rnorm(n*p), n, p) %*% Ve_chol
Y <- X %*% B + G + E

# Create a scaled version with trait 1 having much larger values (multiplied by 50)
Y_scaled <- Y
scale_factor <- 50
Y_scaled[, 1] <- Y_scaled[, 1] * scale_factor

cat("Summary of original data:\n")
cat("Trait 1 mean:", mean(Y[,1]), "sd:", sd(Y[,1]), "\n")
cat("Trait 2 mean:", mean(Y[,2]), "sd:", sd(Y[,2]), "\n")

cat("Summary of scaled data:\n")
cat("Trait 1 mean:", mean(Y_scaled[,1]), "sd:", sd(Y_scaled[,1]), "\n")
cat("Trait 2 mean:", mean(Y_scaled[,2]), "sd:", sd(Y_scaled[,2]), "\n")

# Test with regular (non-standardized) approach on original data
cat("\nTesting original data without standardization:\n")
result1 <- emmremlMultivariate_core(Y = Y, X = X, Z = diag(n), K = K, standardize = FALSE)
cat("Converged:", result1$convergence, "in", result1$iterations, "iterations\n")
cat("Estimated Vg:\n")
print(result1$Vg)
cat("True Vg:\n")
print(Vg)
rel_error_Vg1 <- sum(abs(result1$Vg - Vg)) / sum(abs(Vg))
cat("Relative error in Vg:", rel_error_Vg1, "\n")

# Test with standardization on original data
cat("\nTesting original data with standardization:\n")
result2 <- emmremlMultivariate_core(Y = Y, X = X, Z = diag(n), K = K, standardize = TRUE)
cat("Converged:", result2$convergence, "in", result2$iterations, "iterations\n")
cat("Estimated Vg:\n")
print(result2$Vg)
cat("True Vg:\n")
print(Vg)
rel_error_Vg2 <- sum(abs(result2$Vg - Vg)) / sum(abs(Vg))
cat("Relative error in Vg:", rel_error_Vg2, "\n")

# Test with non-standardized approach on scaled data
cat("\nTesting scaled data without standardization:\n")
result3 <- emmremlMultivariate_core(Y = Y_scaled, X = X, Z = diag(n), K = K, standardize = FALSE)
cat("Converged:", result3$convergence, "in", result3$iterations, "iterations\n")
cat("Estimated Vg:\n")
print(result3$Vg)

# Adjust true Vg for comparison with scaled data
scale_matrix <- diag(c(scale_factor, 1))
adjusted_Vg <- scale_matrix %*% Vg %*% scale_matrix
cat("Adjusted true Vg for scaled data:\n")
print(adjusted_Vg)

rel_error_Vg3 <- sum(abs(result3$Vg - adjusted_Vg)) / sum(abs(adjusted_Vg))
cat("Relative error in Vg compared to adjusted truth:", rel_error_Vg3, "\n")

# Test with standardization on scaled data
cat("\nTesting scaled data with standardization:\n")
result4 <- emmremlMultivariate_core(Y = Y_scaled, X = X, Z = diag(n), K = K, standardize = TRUE)
cat("Converged:", result4$convergence, "in", result4$iterations, "iterations\n")
cat("Estimated Vg:\n")
print(result4$Vg)
cat("True Vg:\n")
print(Vg)
rel_error_Vg4 <- sum(abs(result4$Vg - Vg)) / sum(abs(Vg))
cat("Relative error in Vg:", rel_error_Vg4, "\n")

# Compare results between original and scaled data when using standardization
cat("\nDifference between original and scaled data results with standardization:\n")
cat("Difference in Vg:\n")
print(result2$Vg - result4$Vg)
cat("Relative difference:", sum(abs(result2$Vg - result4$Vg)) / sum(abs(result2$Vg)), "\n")

cat("\nComparison of standardization vs. no standardization on ability to recover true parameters:\n")
cat("Original data - No standardization - Rel error in Vg:", rel_error_Vg1, "\n")
cat("Original data - With standardization - Rel error in Vg:", rel_error_Vg2, "\n")
cat("Scaled data - No standardization - Rel error in Vg (vs adjusted truth):", rel_error_Vg3, "\n") 
cat("Scaled data - With standardization - Rel error in Vg:", rel_error_Vg4, "\n")

cat("\nSummary:\n")
if (rel_error_Vg4 < rel_error_Vg3) {
  cat("PASSED: Standardization improves parameter recovery for scaled data\n")
} else {
  cat("FAILED: Standardization does not improve parameter recovery for scaled data\n")
}

if (sum(abs(result2$Vg - result4$Vg)) / sum(abs(result2$Vg)) < 0.1) {
  cat("PASSED: Standardization yields similar results for original and scaled data\n")
} else {
  cat("FAILED: Standardization does not yield similar results for original and scaled data\n")
}
