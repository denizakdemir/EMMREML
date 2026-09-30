# GWAS Functionality for EMMREML Package

#' Genome-Wide Association Study Analysis
#'
#' @description
#' Performs genome-wide association analysis using mixed models with variance 
#' components estimated by REML. This function can handle both accelerated mode 
#' (where variance components are estimated once) and regular mode (where the 
#' full model is fitted for each marker).
#'
#' @param formula Formula for fixed effects excluding markers, e.g., y ~ x1 + x2
#' @param random Formula for random effects, e.g., ~ id (for kinship)
#' @param data A data frame containing phenotypes and fixed effect covariates
#' @param markers Matrix of marker genotypes (samples in rows, markers in columns)
#' @param K Kinship matrix or NULL (if no kinship correction is needed)
#' @param accelerated Logical, whether to use accelerated mode (default: TRUE)
#' @param pc.count Number of principal components to include as covariates
#' @param method Method for analysis: "mixed" (default), "lm", or "score"
#' @param marker.select Optional index vector for testing only a subset of markers
#' @param verbose Logical, whether to print progress information
#' @param plots Logical, whether to generate plots
#' @param plot.type Types of plots to generate: "manhattan", "qq", "pca", "effect"
#' @param tolParConv Tolerance for parameter convergence
#' @param ... Additional arguments
#'
#' @return An object of class 'emmreml_gwas' containing:
#' \item{results}{Data frame with marker results (effects, p-values, etc.)}
#' \item{manhattan}{Manhattan plot (if plots=TRUE)}
#' \item{qq}{QQ plot (if plots=TRUE)}
#' \item{base.model}{Base model (if accelerated=TRUE)}
#' \item{method}{Method used}
#' \item{call}{The matched call}
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' n <- 200  # number of samples
#' m <- 1000  # number of markers
#' 
#' # Create some markers
#' markers <- matrix(sample(c(0, 1, 2), n*m, replace=TRUE), nrow=n, ncol=m)
#' colnames(markers) <- paste0("SNP", 1:m)
#' 
#' # Create a kinship matrix
#' K <- tcrossprod(scale(markers)) / m
#' 
#' # Simulate phenotypes with causal variants
#' causal_idx <- c(100, 300, 700)
#' betas <- rep(0, m)
#' betas[causal_idx] <- c(0.5, 0.4, 0.3)
#' 
#' # Create genetic values
#' g <- markers %*% betas
#' 
#' # Add random effect from kinship
#' u <- MASS::mvrnorm(1, rep(0, n), 2 * K)
#' 
#' # Add residual error
#' e <- rnorm(n, 0, 1)
#' 
#' # Create phenotype
#' y <- 10 + g + u + e
#' 
#' # Prepare data
#' df <- data.frame(y = y, id = 1:n)
#' 
#' # Run GWAS
#' gwas_result <- gwas(y ~ 1, random = ~ id, data = df, 
#'                    markers = markers, K = K)
#'                    
#' # Plot results
#' plot(gwas_result)
#' }
#'
#' @export
gwas <- function(formula, random = NULL, data, markers, K = NULL, 
                accelerated = TRUE, pc.count = 0, 
                method = c("mixed", "lm", "score"),
                marker.select = NULL, verbose = TRUE, 
                plots = TRUE, plot.type = c("manhattan", "qq", "pca", "effect"),
                tolParConv = 1e-6, ...) {
  
  # Store the call
  call <- match.call()
  
  # Match method and plot type arguments
  method <- match.arg(method)
  plot.type <- match.arg(plot.type, several.ok = TRUE)
  
  # Input validation
  if (is.null(formula) || !inherits(formula, "formula")) stop("Fixed effects formula is required")
  if (missing(data)) stop("Data frame is required")
  if (missing(markers)) stop("Marker matrix is required")
  
  # Check if data is a data frame
  if (!is.data.frame(data)) stop("Data must be a data frame")
  
  # Check markers format
  if (!is.matrix(markers)) {
    if (is.data.frame(markers)) {
      if(verbose) cat("Converting markers from data frame to matrix...\n")
      markers <- as.matrix(markers)
    } else {
      stop("Markers must be a matrix or data frame")
    }
  }
  
  # Get marker dimensions
  n_samples <- nrow(markers)
  n_markers <- ncol(markers)
  if(verbose) cat("Number of samples:", n_samples, "Number of markers:", n_markers, "\n")
  
  # Check if sample sizes match
  if (nrow(data) != n_samples) {
    stop("Number of samples in data and markers must match")
  }
  
  # Check K matrix if provided
  if (!is.null(K)) {
    if (!is.matrix(K)) stop("K must be a matrix")
    if (nrow(K) != n_samples || ncol(K) != n_samples) {
      stop("Dimensions of K must match number of samples")
    }
  }
  
  # Process formulas and extract model matrices
  if (verbose) cat("Processing formulas and preparing model matrices...\n")
  
  # Process fixed effect formula
  mf <- model.frame(formula, data, na.action = na.pass)
  X <- model.matrix(formula, mf)
  y <- model.response(mf)
  
  # Check for missing values
  na_idx <- is.na(y)
  if (any(na_idx)) {
    if (verbose) cat("Removing", sum(na_idx), "rows with missing values in response variable\n")
    y <- y[!na_idx]
    X <- X[!na_idx, , drop = FALSE]
    markers <- markers[!na_idx, , drop = FALSE]
    data <- data[!na_idx, , drop = FALSE]
    if (!is.null(K)) {
      K <- K[!na_idx, !na_idx]
    }
  }
  
  # Update sample count after removing missing values
  n_samples <- nrow(markers)
  
  # Process random effects
  Z <- diag(n_samples)  # Default Z matrix for kinship
  
  # Handle the random effects if provided
  if (!is.null(random)) {
    if (verbose) cat("Processing random effects...\n")
    if (inherits(random, "formula")) {
      Z <- model.matrix(random, data)
    } else {
      stop("'random' must be a formula")
    }
  }
  
  # Apply PCA correction if requested
  if (pc.count > 0) {
    if (verbose) cat("Computing", pc.count, "principal components...\n")
    markers_scaled <- scale(markers)
    pca_result <- prcomp(markers_scaled, center = TRUE, scale. = FALSE)
    pcs <- pca_result$x[, 1:min(pc.count, ncol(pca_result$x))]
    
    # Add PCs as covariates to the design matrix
    if (verbose) cat("Adding PCs as covariates...\n")
    X_with_pcs <- cbind(X, pcs)
    colnames(X_with_pcs)[(ncol(X)+1):ncol(X_with_pcs)] <- paste0("PC", 1:pc.count)
    X <- X_with_pcs
  }
  
  # Set the marker subset if requested
  if (!is.null(marker.select)) {
    if (verbose) cat("Selecting subset of", length(marker.select), "markers...\n")
    markers <- markers[, marker.select, drop = FALSE]
    n_markers <- ncol(markers)
  }
  
  # Prepare result containers
  results <- data.frame(
    marker = colnames(markers),
    effect = numeric(n_markers),
    se = numeric(n_markers),
    p_value = numeric(n_markers),
    stringsAsFactors = FALSE
  )
  
  if (method == "mixed") {
    if (is.null(K)) {
      stop("Kinship matrix required for mixed model analysis")
    }
    
    # Check if K is identity matrix for optimization
    K_is_identity <- is_identity_matrix(K)
    
    if (K_is_identity && verbose) {
      cat("Kinship matrix is identity - simplifying calculations\n")
    }
    
    # Accelerated mode: fit the base model once
    if (accelerated) {
      if (verbose) cat("Fitting base model without markers...\n")
      
      base_model <- emmreml(
        y = y, X = X, Z = Z, K = K,
        varbetahat = TRUE, test = FALSE
      )
      
      # Record variance components for faster marker testing
      Ve <- base_model$Ve
      Vu <- base_model$Vu
      
      # Test each marker
      if (verbose) {
        pb <- txtProgressBar(min = 0, max = n_markers, style = 3)
        cat("Testing", n_markers, "markers in accelerated mode...\n")
      }
      
      for (i in 1:n_markers) {
        # Add marker to fixed effects
        X_with_marker <- cbind(X, markers[, i])
        colnames(X_with_marker)[ncol(X_with_marker)] <- colnames(markers)[i]
        
        # Construct mixed model equations with fixed variance components
        q <- ncol(X_with_marker)
        spI <- diag(n_samples)
        
        # Optimize when K is identity
        if (K_is_identity) {
          ZKZt <- tcrossprod(Z)  # Z*Z' when K is identity
        } else {
          ZK <- Z %*% K
          ZKZt <- tcrossprod(ZK, Z)
        }
        
        # Calculate V matrix with fixed variance components
        V <- Ve * spI + Vu * ZKZt
        Vinv <- solve(V)
        
        # Calculate fixed effects and their variance
        XtVinv <- crossprod(X_with_marker, Vinv)
        C_inv <- solve(XtVinv %*% X_with_marker)
        beta <- C_inv %*% XtVinv %*% y
        beta_var <- diag(C_inv)
        
        # Extract marker effect and test
        marker_idx <- ncol(X_with_marker)
        marker_effect <- beta[marker_idx]
        marker_se <- sqrt(beta_var[marker_idx])
        t_stat <- marker_effect / marker_se
        p_value <- 2 * pnorm(-abs(t_stat))
        
        # Store results
        results$effect[i] <- marker_effect
        results$se[i] <- marker_se
        results$p_value[i] <- p_value
        
        if (verbose) setTxtProgressBar(pb, i)
      }
      
      if (verbose) {
        close(pb)
        cat("\nGWAS analysis complete.\n")
      }
      
    } else {
      # Regular mode: fit full model for each marker
      if (verbose) {
        pb <- txtProgressBar(min = 0, max = n_markers, style = 3)
        cat("Testing", n_markers, "markers with full mixed model...\n")
      }
      
      for (i in 1:n_markers) {
        # Add marker to fixed effects
        X_with_marker <- cbind(X, markers[, i])
        colnames(X_with_marker)[ncol(X_with_marker)] <- colnames(markers)[i]
        
        # Fit full model
        model <- emmreml(
          y = y, X = X_with_marker, Z = Z, K = K,
          varbetahat = TRUE, test = TRUE
        )
        
        # Extract marker effect and test
        marker_idx <- ncol(X_with_marker)
        marker_effect <- model$betahat[marker_idx]
        marker_se <- sqrt(model$varbetahat[marker_idx])
        
        # Use test results directly if available
        if (!is.null(model$pvalbeta)) {
          p_value <- model$pvalbeta[marker_idx, "none"]  # Unadjusted p-value
        } else {
          t_stat <- marker_effect / marker_se
          p_value <- 2 * pnorm(-abs(t_stat))
        }
        
        # Store results
        results$effect[i] <- marker_effect
        results$se[i] <- marker_se
        results$p_value[i] <- p_value
        
        if (verbose) setTxtProgressBar(pb, i)
      }
      
      if (verbose) {
        close(pb)
        cat("\nGWAS analysis complete.\n")
      }
      
      # Base model is NULL in regular mode
      base_model <- NULL
    }
    
  } else if (method == "lm") {
    # Linear model approach
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = n_markers, style = 3)
      cat("Testing", n_markers, "markers with linear model...\n")
    }
    
    # Pre-compute X'X inverse
    XtX_inv <- solve(crossprod(X))
    Xty <- crossprod(X, y)
    
    for (i in 1:n_markers) {
      # Add marker to fixed effects
      X_with_marker <- cbind(X, markers[, i])
      
      # Fit linear model
      XtX_with_marker <- crossprod(X_with_marker)
      Xty_with_marker <- crossprod(X_with_marker, y)
      
      # Calculate coefficients and their variance
      beta <- solve(XtX_with_marker) %*% Xty_with_marker
      residual <- y - X_with_marker %*% beta
      sigma2 <- sum(residual^2) / (n_samples - ncol(X_with_marker))
      var_beta <- sigma2 * solve(XtX_with_marker)
      
      # Extract marker effect and test
      marker_idx <- ncol(X_with_marker)
      marker_effect <- beta[marker_idx]
      marker_se <- sqrt(diag(var_beta)[marker_idx])
      t_stat <- marker_effect / marker_se
      p_value <- 2 * pt(-abs(t_stat), df = n_samples - ncol(X_with_marker))
      
      # Store results
      results$effect[i] <- marker_effect
      results$se[i] <- marker_se
      results$p_value[i] <- p_value
      
      if (verbose) setTxtProgressBar(pb, i)
    }
    
    if (verbose) {
      close(pb)
      cat("\nGWAS analysis complete.\n")
    }
    
    base_model <- NULL
    
  } else if (method == "score") {
    # Score test approach
    if (verbose) cat("Score test not yet implemented, using mixed model instead\n")
    # TODO: Implement score test
    
    base_model <- NULL
  }
  
  # Calculate Bonferroni and FDR adjusted p-values
  results$p_bonf <- p.adjust(results$p_value, method = "bonferroni")
  results$p_fdr <- p.adjust(results$p_value, method = "fdr")
  
  # Create result object
  gwas_result <- list(
    results = results,
    method = method,
    accelerated = accelerated,
    base_model = base_model,
    call = call
  )
  
  # Generate plots if requested
  if (plots) {
    plot_list <- list()
    
    if ("manhattan" %in% plot.type) {
      if (verbose) cat("Generating Manhattan plot...\n")
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        # Try to get chromosome information from marker names
        if (!is.null(colnames(markers))) {
          # Check for patterns like "chr1_12345" or "1:12345"
          chr_pattern <- "^(chr)?([0-9XY]+)[_:-]([0-9]+)$"
          marker_names <- colnames(markers)
          chr <- rep(NA, length(marker_names))
          pos <- rep(NA, length(marker_names))
          
          for (i in 1:length(marker_names)) {
            match <- regexpr(chr_pattern, marker_names[i], perl = TRUE)
            if (match != -1) {
              parts <- strsplit(marker_names[i], "[_:-]")[[1]]
              chr[i] <- gsub("chr", "", parts[1])
              pos[i] <- as.numeric(parts[length(parts)])
            }
          }
          
          results$chr <- chr
          results$pos <- pos
        }
        
        # Create Manhattan plot
        if (all(is.na(results$chr)) || all(is.na(results$pos))) {
          # No chromosome/position info, use marker index
          results$chr <- "1"
          results$pos <- 1:nrow(results)
        }
        
        # Make sure chr is a factor with proper order
        results$chr <- factor(results$chr, levels = unique(results$chr[order(as.numeric(gsub("[^0-9]", "", results$chr)))]))
        
        # Add color column
        results$color <- as.integer(results$chr) %% 2 + 1
        
        # Bonferroni threshold
        bonf_threshold <- 0.05 / nrow(results)
        
        # Create Manhattan plot
        manhattan_plot <- ggplot2::ggplot(results, ggplot2::aes(x = pos, y = -log10(p_value), 
                                                   color = factor(color))) +
          ggplot2::geom_point(alpha = 0.8) +
          ggplot2::geom_hline(yintercept = -log10(bonf_threshold), linetype = "dashed", color = "red") +
          ggplot2::facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "none") +
          ggplot2::labs(x = "Position", y = "-log10(p-value)", 
               title = "Manhattan Plot")
        
        plot_list$manhattan <- manhattan_plot
      } else {
        warning("ggplot2 package is required for Manhattan plot")
      }
    }
    
    if ("qq" %in% plot.type) {
      if (verbose) cat("Generating QQ plot...\n")
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        # Calculate expected p-values
        observed <- -log10(sort(results$p_value))
        expected <- -log10(ppoints(length(observed)))
        qq_data <- data.frame(expected = expected, observed = observed)
        
        # Create QQ plot
        qq_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(x = expected, y = observed)) +
          ggplot2::geom_point(alpha = 0.8) +
          ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
          ggplot2::theme_minimal() +
          ggplot2::labs(x = "Expected -log10(p-value)", y = "Observed -log10(p-value)", 
               title = "Quantile-Quantile Plot")
        
        plot_list$qq <- qq_plot
      } else {
        warning("ggplot2 package is required for QQ plot")
      }
    }
    
    if ("pca" %in% plot.type && pc.count >= 2) {
      if (verbose) cat("Generating PCA plot...\n")
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        # Use PCs calculated earlier
        if (!exists("pca_result")) {
          markers_scaled <- scale(markers)
          pca_result <- prcomp(markers_scaled, center = TRUE, scale. = FALSE)
          pcs <- pca_result$x[, 1:min(pc.count, ncol(pca_result$x))]
        }
        
        # Create PCA plot with first two components
        pca_data <- data.frame(PC1 = pcs[,1], PC2 = pcs[,2], y = y)
        pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = y)) +
          ggplot2::geom_point(alpha = 0.8) +
          ggplot2::theme_minimal() +
          ggplot2::scale_color_gradient(low = "blue", high = "red") +
          ggplot2::labs(x = "PC1", y = "PC2", color = "Phenotype", 
               title = "PCA Plot of Genotypes")
        
        plot_list$pca <- pca_plot
      } else {
        warning("ggplot2 package is required for PCA plot")
      }
    }
    
    if ("effect" %in% plot.type) {
      if (verbose) cat("Generating effect size plot...\n")
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        # Create effect size plot
        effects_data <- results[order(abs(results$effect), decreasing = TRUE)[1:min(20, nrow(results))], ]
        
        effect_plot <- ggplot2::ggplot(effects_data, 
                                       ggplot2::aes(x = reorder(marker, abs(effect)), y = effect)) +
          ggplot2::geom_bar(stat = "identity", alpha = 0.8) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = effect - 1.96*se, ymax = effect + 1.96*se), width = 0.2) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(x = "Marker", y = "Effect Size", 
               title = "Top Marker Effects with 95% CI")
        
        plot_list$effect <- effect_plot
      } else {
        warning("ggplot2 package is required for effect size plot")
      }
    }
    
    gwas_result$plots <- plot_list
  }
  
  # Set class
  class(gwas_result) <- "emmreml_gwas"
  
  return(gwas_result)
}

#' Plot Method for emmreml_gwas Objects
#'
#' @param x An object of class 'emmreml_gwas'
#' @param type Type of plot to display: "manhattan", "qq", "pca", or "effect"
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns the plot object
#'
#' @method plot emmreml_gwas
#' @export
plot.emmreml_gwas <- function(x, type = c("manhattan", "qq", "pca", "effect"), ...) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for plotting")
  }
  
  if (is.null(x$plots) || is.null(x$plots[[type]])) {
    stop(paste("Plot of type '", type, "' not available in this GWAS result", sep = ""))
  }
  
  print(x$plots[[type]])
  invisible(x$plots[[type]])
}

#' Print Method for emmreml_gwas Objects
#'
#' @param x An object of class 'emmreml_gwas'
#' @param ... Additional arguments
#'
#' @return Invisibly returns the input object
#'
#' @method print emmreml_gwas
#' @export
print.emmreml_gwas <- function(x, ...) {
  cat("GWAS analysis with EMMREML\n")
  cat("Method:", x$method, ifelse(x$accelerated, "(accelerated)", ""), "\n")
  cat("Number of markers:", nrow(x$results), "\n")
  
  # Print top hits
  top_n <- min(5, nrow(x$results))
  cat("\nTop", top_n, "markers by p-value:\n")
  top_hits <- x$results[order(x$results$p_value)[1:top_n], c("marker", "effect", "se", "p_value")]
  print(top_hits)
  
  invisible(x)
}

#' Summary Method for emmreml_gwas Objects
#'
#' @param object An object of class 'emmreml_gwas'
#' @param ... Additional arguments
#'
#' @return A summary object
#'
#' @method summary emmreml_gwas
#' @export
summary.emmreml_gwas <- function(object, ...) {
  # Create a summary object
  result <- list(
    call = object$call,
    method = object$method,
    accelerated = object$accelerated,
    n_markers = nrow(object$results),
    top_hits = object$results[order(object$results$p_value)[1:min(10, nrow(object$results))], ]
  )
  
  # Calculate lambda (genomic inflation factor)
  observed_chisq <- qchisq(object$results$p_value, df = 1, lower.tail = FALSE)
  expected_median <- qchisq(0.5, df = 1)
  result$lambda <- median(observed_chisq) / expected_median
  
  # Count significant hits
  result$bonf_hits <- sum(object$results$p_bonf < 0.05)
  result$fdr_hits <- sum(object$results$p_fdr < 0.05)
  
  class(result) <- "summary.emmreml_gwas"
  return(result)
}

#' Print Method for summary.emmreml_gwas Objects
#'
#' @param x An object of class 'summary.emmreml_gwas'
#' @param ... Additional arguments
#'
#' @return Invisibly returns the input object
#'
#' @method print summary.emmreml_gwas
#' @export
print.summary.emmreml_gwas <- function(x, ...) {
  cat("GWAS analysis with EMMREML\n")
  cat("\nCall:\n")
  print(x$call)
  
  cat("\nMethod:", x$method, ifelse(x$accelerated, "(accelerated)", ""), "\n")
  cat("Number of markers:", x$n_markers, "\n")
  cat("Genomic inflation factor (lambda):", round(x$lambda, 3), "\n")
  cat("Significant markers (Bonferroni p < 0.05):", x$bonf_hits, "\n")
  cat("Significant markers (FDR < 0.05):", x$fdr_hits, "\n")
  
  cat("\nTop", nrow(x$top_hits), "markers by p-value:\n")
  print(x$top_hits)
  
  invisible(x)
}
