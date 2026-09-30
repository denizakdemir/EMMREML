# Examples of using EMMREML package

#' Generate example data for EMMREML demonstrations
#'
#' @param n Number of samples
#' @param p Number of markers (for GWAS examples)
#' @param with_identity_K Whether to create an identity K matrix example 
#' @param with_gwas Whether to generate GWAS example data
#'
#' @return A list with example data
#'
#' @examples
#' \dontrun{
#' # Basic mixed model
#' example_data <- generate_example_data(n = 100)
#' 
#' # Fit model
#' mod <- emmreml(y ~ x, random = ~ id, data = example_data$data, K = example_data$K)
#' summary(mod)
#' }
#'
#' @export
generate_example_data <- function(n = 100, p = 1000, with_identity_K = TRUE, with_gwas = TRUE) {
  # Basic example data
  set.seed(123)
  
  # Create fixed effects
  X <- cbind(1, rnorm(n))
  colnames(X) <- c("intercept", "x")
  
  # Create a realistic kinship matrix
  Z <- diag(n)
  # Generate a matrix that's guaranteed to be positive definite
  K_raw <- matrix(rnorm(n*20), n, 20)
  K <- tcrossprod(K_raw)
  # Add a small value to diagonal to ensure positive definiteness
  diag(K) <- diag(K) + 1e-6
  K <- K / mean(diag(K))  # Scale K
  
  # Create scaled identity K for optimization examples
  if (with_identity_K) {
    K_identity <- diag(n)
  }
  
  # Simulate data from the model
  u <- crossprod(chol(K), rnorm(n)) * sqrt(2)  # random effect with Var=2*K
  e <- rnorm(n)  # residual with Var=1
  y <- 10 + 2 * X[, 2] + u + e
  
  # Create data frame
  df <- data.frame(y = y, x = X[, 2], id = 1:n)
  
  # Generate GWAS data if requested
  if (with_gwas) {
    # Create markers
    markers <- matrix(sample(c(0, 1, 2), n*p, replace=TRUE), nrow=n, ncol=p)
    colnames(markers) <- paste0("SNP", 1:p)
    
    # Select some causal variants
    causal_idx <- c(100, 300, 700)
    betas <- rep(0, p)
    betas[causal_idx] <- c(0.5, 0.4, 0.3)
    
    # Create genetic values from markers
    g <- markers %*% betas
    
    # Create phenotype for GWAS
    y_gwas <- 10 + g + u + e
    df_gwas <- data.frame(y = y_gwas, x = X[, 2], id = 1:n)
    
    # Return all data
    return(list(
      data = df,
      data_gwas = df_gwas,
      X = X,
      Z = Z,
      K = K,
      K_identity = if(with_identity_K) K_identity else NULL,
      markers = markers,
      causal_idx = causal_idx
    ))
  } else {
    # Return only basic data
    return(list(
      data = df,
      X = X,
      Z = Z,
      K = K,
      K_identity = if(with_identity_K) K_identity else NULL
    ))
  }
}

#' Example of using EMMREML with K optimization
#'
#' This function demonstrates how EMMREML can automatically detect when K is an 
#' identity matrix and optimize calculations accordingly.
#'
#' @return Nothing, runs examples and prints timing results
#'
#' @examples
#' \dontrun{
#' demo_K_optimization()
#' }
#'
#' @export
demo_K_optimization <- function() {
  # Generate example data
  cat("Generating example data...\n")
  example_data <- generate_example_data(n = 200)
  
  # Run with regular K
  cat("\nRunning with regular K matrix...\n")
  time_regular <- system.time({
    # Extract needed components manually for core function
    y <- example_data$data$y
    X <- cbind(1, example_data$data$x)  # Intercept + x
    Z <- diag(nrow(example_data$data))  # Identity for observations
    
    mod_regular <- emmreml_core(y = y, X = X, Z = Z, K = example_data$K)
  })
  
  # Run with identity K
  cat("\nRunning with identity K matrix (should be faster)...\n")
  time_identity <- system.time({
    # Extract needed components manually for core function
    y <- example_data$data$y
    X <- cbind(1, example_data$data$x)  # Intercept + x
    Z <- diag(nrow(example_data$data))  # Identity for observations
    
    mod_identity <- emmreml_core(y = y, X = X, Z = Z, K = example_data$K_identity)
  })
  
  # Compare timing
  cat("\nTiming comparison:\n")
  cat("Regular K time:", time_regular["elapsed"], "seconds\n")
  cat("Identity K time:", time_identity["elapsed"], "seconds\n")
  cat("Speedup factor:", time_regular["elapsed"] / time_identity["elapsed"], "\n")
  
  # Compare results
  cat("\nResults comparison:\n")
  cat("Regular K: Ve =", mod_regular$Ve, "Vu =", mod_regular$Vu, "\n")
  cat("Identity K: Ve =", mod_identity$Ve, "Vu =", mod_identity$Vu, "\n")
  
  invisible(NULL)
}

#' Example of using GWAS functionality
#'
#' This function demonstrates how to use the GWAS functionality in EMMREML.
#'
#' @return GWAS result object
#'
#' @examples
#' \dontrun{
#' gwas_result <- demo_gwas()
#' plot(gwas_result)
#' }
#'
#' @export
demo_gwas <- function() {
  # Generate example data
  cat("Generating GWAS example data...\n")
  example_data <- generate_example_data(n = 100, p = 500)
  
  # Run GWAS with accelerated mode
  cat("\nRunning GWAS with accelerated mode...\n")
  gwas_result <- gwas(
    formula = y ~ x, 
    random = ~ id, 
    data = example_data$data_gwas,
    markers = example_data$markers,
    K = example_data$K,
    accelerated = TRUE,
    verbose = TRUE,
    plots = TRUE
  )
  
  # Print summary
  summary(gwas_result)
  
  # Highlight causal variants in results
  cat("\nTrue causal variants:\n")
  causal_idx <- example_data$causal_idx
  causal_results <- gwas_result$results[causal_idx, ]
  print(causal_results)
  
  # Get rankings of true causal variants
  rankings <- sapply(causal_idx, function(idx) {
    which(order(gwas_result$results$p_value) == idx)
  })
  cat("\nRankings of causal variants:", rankings, "\n")
  
  # Return the GWAS result
  return(gwas_result)
}

#' Example of using the multivariate EMMREML
#'
#' This function demonstrates how to use the multivariate EMMREML functionality.
#'
#' @return Multivariate EMMREML result object
#'
#' @examples
#' \dontrun{
#' mv_result <- demo_multivariate()
#' summary(mv_result)
#' }
#'
#' @export
demo_multivariate <- function() {
  # Generate example data
  cat("Generating multivariate example data...\n")
  n <- 100
  p <- 3  # number of traits
  
  # Create fixed effects
  X <- cbind(1, rnorm(n))
  colnames(X) <- c("intercept", "x")
  
  # Create random effects
  Z <- diag(n)
  # Generate a matrix that's guaranteed to be positive definite
  K_raw <- matrix(rnorm(n*20), n, 20)
  K <- tcrossprod(K_raw)
  # Add a small value to diagonal to ensure positive definiteness
  diag(K) <- diag(K) + 1e-6
  K <- K/mean(diag(K))
  
  # Create covariance matrices
  Vg <- matrix(c(2, 1, 0.5, 1, 1.5, 0.75, 0.5, 0.75, 1), p, p)
  Ve <- matrix(c(1, 0.2, 0.1, 0.2, 0.8, 0.15, 0.1, 0.15, 0.5), p, p)
  
  # Simulate data
  B <- matrix(c(10, 2, 5, 1, 3, 0.5), 2, p)  # fixed effects for 3 traits
  G <- matrix(rnorm(n*p), n, p) %*% chol(Vg)  # random effects with Vg covariance
  E <- matrix(rnorm(n*p), n, p) %*% chol(Ve)  # errors with Ve covariance
  Y <- X %*% B + G + E
  
  # Prepare data frame
  df <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], x = X[,2], id = 1:n)
  
  # Fit multivariate model
  cat("\nFitting multivariate model...\n")
  mv_result <- emmremlMultivariate(
    formula = cbind(y1, y2, y3) ~ x, 
    random = ~ id, 
    data = df, 
    K = K,
    verbose = TRUE
  )
  
  # Print results
  cat("\nResults:\n")
  print(mv_result)
  cat("\nGenetic covariance matrix (Vg) - True:\n")
  print(Vg)
  cat("\nGenetic covariance matrix (Vg) - Estimated:\n")
  print(mv_result$Vg)
  cat("\nError covariance matrix (Ve) - True:\n")
  print(Ve)
  cat("\nError covariance matrix (Ve) - Estimated:\n")
  print(mv_result$Ve)
  
  return(mv_result)
}

#' Demonstrate the kinship matrix calculation functionality
#'
#' This function demonstrates how to use the kinship module to process marker data,
#' calculate different types of kinship matrices, and visualize the results.
#'
#' @return A list containing the example marker data, processed data, and kinship matrices
#'
#' @examples
#' \dontrun{
#' demo_results <- demo_kinship()
#' str(demo_results)
#' }
#'
#' @export
demo_kinship <- function() {
  # Check required packages
  has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
  if (!has_ggplot2) {
    message("The 'ggplot2' package is recommended for visualizations. ",
            "Install it with install.packages('ggplot2')")
  }
  
  # Generate example marker data
  cat("Generating example marker data...\n")
  set.seed(123)
  n_samples <- 100
  n_markers <- 1000
  
  # Create marker data with 0,1,2 coding (SNP data)
  cat("Creating marker data with 0,1,2 coding...\n")
  M_012 <- matrix(sample(c(0,1,2, NA), n_samples*n_markers, replace=TRUE, 
                        prob=c(0.4, 0.3, 0.2, 0.1)), 
                 nrow=n_samples, ncol=n_markers)
  colnames(M_012) <- paste0("SNP", 1:n_markers)
  rownames(M_012) <- paste0("Sample", 1:n_samples)
  
  # Identify marker coding
  cat("\nIdentifying marker coding...\n")
  coding_info <- identify_marker_coding(M_012)
  print(coding_info)
  
  # Calculate marker statistics
  cat("\nCalculating marker statistics...\n")
  stats <- marker_statistics(M_012, coding="012")
  cat("First few rows of marker statistics:\n")
  print(head(stats))
  
  # Create a summary table
  cat("\nCreating summary table...\n")
  summary_table <- summarize_marker_data(M_012, coding="012")
  print(summary_table)
  
  # Plot marker statistics if ggplot2 is available
  if (has_ggplot2) {
    cat("\nCreating marker statistics plots...\n")
    plots <- plot_marker_stats(M_012, coding="012", plot_type="both")
    cat("Created MAF and missing data distribution plots (not shown in console)\n")
  }
  
  # Filter markers
  cat("\nFiltering markers...\n")
  M_filtered <- filter_markers(M_012, min_maf=0.05, max_missing=0.05)
  
  # Impute missing values
  cat("\nImputing missing values...\n")
  M_imputed <- impute_markers(M_filtered, method="mean", coding="012")
  
  # Create different versions of the marker data
  cat("\nCreating different versions of the marker data...\n")
  M_n101 <- convert_marker_coding(M_imputed, from="012", to="-101")
  cat("Created M_n101 with -1,0,1 coding\n")
  
  # Calculate kinship matrices using different methods
  cat("\nCalculating kinship matrices using different methods...\n")
  
  # VanRaden method 1
  cat("- VanRaden method 1...\n")
  K_vr1 <- calculate_kinship(M_imputed, method="vanraden1", coding="012")
  
  # VanRaden method 2
  cat("- VanRaden method 2...\n")
  K_vr2 <- calculate_kinship(M_imputed, method="vanraden2", coding="012")
  
  # IBS method
  cat("- IBS method...\n")
  # K_ibs <- calculate_kinship(M_imputed, method="ibs", coding="012")
  
  # Correlation-based
  cat("- Correlation-based method...\n")
  K_cor <- calculate_kinship(M_imputed, method="correlation", coding="012")
  
  # Compare kinship matrices
  cat("\nComparing kinship matrices...\n")
  if (has_ggplot2) {
    comparison <- compare_kinship(
      list(K_vr1, K_vr2, K_ibs, K_cor),
      names=c("VanRaden1", "VanRaden2", "IBS", "Correlation")
    )
    cat("Created kinship comparison heatmap (not shown in console)\n")
  } else {
    cat("Correlation between kinship matrices:\n")
    comparison <- compare_kinship(
      list(K_vr1, K_vr2, K_ibs, K_cor),
      names=c("VanRaden1", "VanRaden2", "IBS", "Correlation"),
      plot=FALSE
    )
    print(comparison)
  }
  
  # Visualize one of the kinship matrices
  if (has_ggplot2) {
    cat("\nVisualizing VanRaden1 kinship matrix...\n")
    kinship_viz <- visualize_kinship(
      K_vr1, 
      clustering=TRUE,
      title="VanRaden Method 1 Kinship Matrix"
    )
    cat("Created kinship heatmap (not shown in console)\n")
  }
  
  # Demonstrate usage with EMMREML
  cat("\nDemonstrating usage with EMMREML...\n")
  
  # Make sure the kinship matrix is positive definite
  cat("Ensuring positive definiteness of the kinship matrix...\n")
  # Add a small value to the diagonal to ensure positive-definiteness
  diag(K_vr1) <- diag(K_vr1) + 1e-4
  
  # Simulate a phenotype based on the kinship
  set.seed(456)
  # Using eigen decomposition instead of Cholesky for robustness
  eig <- eigen(K_vr1)
  sqrt_K <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
  u <- sqrt_K %*% rnorm(n_samples) * sqrt(2)  # random effect with Var=2*K
  x <- rnorm(n_samples)  # covariate
  e <- rnorm(n_samples)  # residual with Var=1
  y <- 10 + 2 * x + u + e
  
  # Create data frame
  df <- data.frame(y = y, x = x, id = 1:n_samples)
  
  # Fit mixed model - Using core implementation directly
  cat("Fitting mixed model with VanRaden kinship...\n")
  # Prepare model matrices 
  model_X <- cbind(1, df$x)  # Intercept and x
  model_Z <- diag(nrow(df))  # Identity matrix for each observation
  
  # Call emmreml_core directly
  mod <- emmreml_core(y = df$y, X = model_X, Z = model_Z, K = K_vr1)
  
  # Print results
  cat("\nModel results:\n")
  print(summary(mod))
  
  # Return results
  return(list(
    M_012 = M_012,
    M_filtered = M_filtered,
    M_imputed = M_imputed,
    M_n101 = M_n101,
    kinship = list(
      vanraden1 = K_vr1,
      vanraden2 = K_vr2,
      ibs = K_ibs,
      correlation = K_cor
    ),
    model = mod
  ))
}
