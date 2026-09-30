# Test script for GWAS functionality

library(EMMREML)

# Set seed for reproducibility
set.seed(123)

cat("\n==== Testing GWAS functionality ====\n")

# Generate a simple example dataset for GWAS testing
cat("Generating example data...\n")
n <- 50  # small sample size for testing
p <- 200  # small number of markers for testing

# Create fixed effects
X <- cbind(1, rnorm(n))
colnames(X) <- c("intercept", "x")

# Create a kinship matrix
Z <- diag(n)
K_raw <- matrix(rnorm(n*20), n, 20)
K <- tcrossprod(K_raw)
diag(K) <- diag(K) + 1e-6  # ensure positive definiteness
K <- K / mean(diag(K))  # Scale K

# Create markers
markers <- matrix(sample(c(0, 1, 2), n*p, replace=TRUE), nrow=n, ncol=p)
colnames(markers) <- paste0("SNP", 1:p)

# Select some causal variants
causal_idx <- c(20, 60, 140)
betas <- rep(0, p)
betas[causal_idx] <- c(0.5, 0.4, 0.3)

# Create genetic values from markers
g <- markers %*% betas

# Create phenotype
u <- crossprod(chol(K), rnorm(n)) * sqrt(2)  # random effect with Var=2*K
e <- rnorm(n)  # residual with Var=1
y_gwas <- 10 + g + u + e

# Create data frame
df_gwas <- data.frame(y = y_gwas, x = X[, 2], id = 1:n)

# Print dimensions for debugging
cat("Data dimensions check:\n")
cat("Data frame rows:", nrow(df_gwas), "\n")
cat("Markers matrix:", nrow(markers), "x", ncol(markers), "\n")
cat("Kinship matrix:", nrow(K), "x", ncol(K), "\n")

# Test GWAS function with explicit arguments to avoid formula parsing issues
cat("\nRunning GWAS with accelerated mode...\n")
tryCatch({
  # Explicitly use model matrices
  gwas_result <- gwas(
    formula = y ~ x,  # Simple formula
    random = ~ id,    # Random effect
    data = df_gwas,   # Data frame
    markers = markers,  # Marker matrix
    K = K,              # Kinship matrix
    accelerated = TRUE,  # Use accelerated mode
    verbose = TRUE,      # Show progress
    plots = FALSE        # Skip plotting for now
  )
  
  # Display results
  cat("\nGWAS results summary:\n")
  print(summary(gwas_result))
  
  # Check if causal variants were detected
  cat("\nRankings of causal variants:\n")
  for (idx in causal_idx) {
    variant_name <- paste0("SNP", idx)
    rank <- which(order(gwas_result$results$p_value) == which(colnames(markers) == variant_name))
    cat(variant_name, "rank:", rank, "\n")
  }
  
}, error = function(e) {
  cat("Error in GWAS:", conditionMessage(e), "\n")
  cat("Traceback:\n")
  print(sys.calls())
})

cat("\n==== GWAS test completed ====\n")
