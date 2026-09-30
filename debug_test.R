# Simplified test case to debug the multivariate REML implementation

library(EMMREML)
options(warn = 1)  # Print warnings immediately

# Create very simple test data
n <- 10  # small sample size for easier debugging
p <- 2   # number of traits
m <- 2   # number of fixed effects

# Fixed effects
X <- cbind(1, 1:n/n)
colnames(X) <- c("intercept", "x")

# Create simple kinship matrix - identity for simplicity
K <- diag(n)

# Create Z matrix
Z <- diag(n)

# Create simple phenotype matrix
Y <- matrix(rnorm(n*p), n, p)

# Manually walk through the emmremlMultivariate_core function
# Use Trace-based approach to find where the code is failing

cat("Input dimensions:\n")
cat("Y:", dim(Y)[1], "x", dim(Y)[2], "\n")
cat("X:", dim(X)[1], "x", dim(X)[2], "\n")
cat("Z:", dim(Z)[1], "x", dim(Z)[2], "\n")
cat("K:", dim(K)[1], "x", dim(K)[2], "\n")

# Step 1: Compute core matrices
ZKZt <- Z %*% K %*% t(Z)
cat("ZKZt dimensions:", dim(ZKZt)[1], "x", dim(ZKZt)[2], "\n")

# Step 2: Initial covariance matrices
raw_cov <- cov(t(Y))
Vg <- raw_cov / 2
Ve <- raw_cov / 2
cat("Vg dimensions:", dim(Vg)[1], "x", dim(Vg)[2], "\n")
cat("Ve dimensions:", dim(Ve)[1], "x", dim(Ve)[2], "\n")

# Step 3: Create the phenotypic covariance matrix
V <- kronecker(ZKZt, Vg) + kronecker(diag(n), Ve)
cat("V dimensions:", dim(V)[1], "x", dim(V)[2], "\n")

# Step 4: Build the multivariate design matrix
X_mv <- matrix(0, n*p, m*p)
for (i in 1:p) {
  rows_i <- ((i-1)*n+1):(i*n)
  cols_i <- ((i-1)*m+1):(i*m)
  X_mv[rows_i, cols_i] <- X
}
cat("X_mv dimensions:", dim(X_mv)[1], "x", dim(X_mv)[2], "\n")

# Step 5: Compute inverse and check matrix compatibility
Vinv <- solve(V + diag(1e-6, nrow(V)))
cat("Vinv dimensions:", dim(Vinv)[1], "x", dim(Vinv)[2], "\n")

# Check that matrix dimensions are compatible for multiplication
XtVinv <- t(X_mv) %*% Vinv
cat("XtVinv dimensions:", dim(XtVinv)[1], "x", dim(XtVinv)[2], "\n")

cat("\nModified Debug Test:\n")

# Now try the fixed multivariate function with simple data
tryCatch({
  result <- emmremlMultivariate_core(
    Y = Y, 
    X = X, 
    Z = Z, 
    K = K,
    maxIter = 5,
    standardize = TRUE
  )
  
  cat("Test successful\!\n")
  cat("Iterations:", result$iterations, "\n")
  cat("Convergence:", result$convergence, "\n")
  
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")
  cat("Test failed\n")
})
