#' @title Kinship Matrix Calculation and Marker Data Processing
#' @description Functions for calculating relationship matrices from marker data,
#'   marker data preprocessing, filtering, imputation, and visualization.

#' Identify marker data coding type
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param check_samples Integer, number of markers to sample when checking (default: 100)
#' @param tolerance Numeric, tolerance for comparing values (default: 1e-8)
#'
#' @return A list with elements:
#'   \item{coding_type}{Character, one of "01", "012", "-101", or "custom"}
#'   \item{min_val}{Minimum value in the matrix}
#'   \item{max_val}{Maximum value in the matrix}
#'   \item{unique_values}{Vector of unique values found}
#'   \item{ploidy}{Estimated ploidy level based on unique values}
#'
#' @examples
#' \dontrun{
#' # Binary markers (0,1)
#' M_01 <- matrix(sample(c(0,1), 100*10, replace=TRUE), 100, 10)
#' identify_marker_coding(M_01)
#'
#' # SNP markers (0,1,2)
#' M_012 <- matrix(sample(c(0,1,2), 100*10, replace=TRUE), 100, 10)
#' identify_marker_coding(M_012)
#'
#' # Centered markers (-1,0,1)
#' M_n101 <- matrix(sample(c(-1,0,1), 100*10, replace=TRUE), 100, 10)
#' identify_marker_coding(M_n101)
#' }
#'
#' @export
identify_marker_coding <- function(M, check_samples = 100, tolerance = 1e-8) {
  # Get dimensions
  n_samples <- nrow(M)
  n_markers <- ncol(M)
  
  # Sample columns to check (or use all if fewer than check_samples)
  cols_to_check <- min(n_markers, check_samples)
  sampled_cols <- sample(n_markers, cols_to_check)
  
  # Get unique values across sampled columns
  unique_vals <- sort(unique(as.vector(M[, sampled_cols])))
  
  # Get min and max values
  min_val <- min(unique_vals)
  max_val <- max(unique_vals)
  
  # Identify coding type
  if (all(abs(unique_vals - round(unique_vals)) < tolerance)) {
    # Integer coding
    if (min_val >= 0) {
      if (max_val <= 1 && all(unique_vals %in% c(0, 1))) {
        coding_type <- "01"
        ploidy <- 1
      } else if (max_val <= 2 && all(unique_vals %in% c(0, 1, 2))) {
        coding_type <- "012"
        ploidy <- 2
      } else {
        # Check if it's a higher ploidy level
        if (length(unique_vals) == max_val + 1 && 
            all(unique_vals == seq(0, max_val))) {
          coding_type <- paste0("0-", max_val)
          ploidy <- max_val
        } else {
          coding_type <- "custom"
          ploidy <- NA
        }
      }
    } else if (min_val == -1 && max_val == 1 && 
               all(unique_vals %in% c(-1, 0, 1))) {
      coding_type <- "-101"
      ploidy <- 2
    } else {
      coding_type <- "custom"
      ploidy <- NA
    }
  } else {
    # Non-integer coding (possibly scaled or normalized)
    coding_type <- "custom"
    ploidy <- NA
  }
  
  return(list(
    coding_type = coding_type,
    min_val = min_val,
    max_val = max_val,
    unique_values = unique_vals,
    ploidy = ploidy
  ))
}

#' Convert marker coding between different types
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param from Character, the current coding type ("01", "012", "-101", or "custom")
#' @param to Character, the desired coding type ("01", "012", "-101")
#'
#' @return A marker matrix with the converted coding
#'
#' @examples
#' \dontrun{
#' # Convert 0,1,2 coding to -1,0,1
#' M_012 <- matrix(sample(c(0,1,2), 100*10, replace=TRUE), 100, 10)
#' M_n101 <- convert_marker_coding(M_012, from="012", to="-101")
#' }
#'
#' @export
convert_marker_coding <- function(M, from, to) {
  if (from == to) {
    return(M)
  }
  
  M_new <- M
  
  if (from == "01" && to == "012") {
    warning("Converting from binary (0,1) to SNP (0,1,2) format is ambiguous. ",
            "Treating 0 as 0/0 and 1 as 1/1 (skipping heterozygotes).")
    M_new <- M_new * 2
  } else if (from == "01" && to == "-101") {
    warning("Converting from binary (0,1) to centered (-1,0,1) format is ambiguous. ",
            "Treating 0 as -1 and 1 as 1 (centered around 0).")
    M_new <- M_new * 2 - 1
  } else if (from == "012" && to == "01") {
    warning("Converting from SNP (0,1,2) to binary (0,1) format loses information. ",
            "Treating 0 and 1 as 0, and 2 as 1.")
    M_new <- ifelse(M_new >= 2, 1, 0)
  } else if (from == "012" && to == "-101") {
    # Center around 1 (the mean for 0,1,2 coded markers)
    M_new <- M_new - 1
  } else if (from == "-101" && to == "01") {
    warning("Converting from centered (-1,0,1) to binary (0,1) format loses information. ",
            "Treating -1 as 0, and 0 and 1 as 1.")
    M_new <- ifelse(M_new >= 0, 1, 0)
  } else if (from == "-101" && to == "012") {
    # Shift to make minimum 0
    M_new <- M_new + 1
  } else {
    stop("Unsupported conversion from '", from, "' to '", to, "'")
  }
  
  return(M_new)
}

#' Calculate allele frequencies from marker data
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#' @param na.rm Logical, whether to remove NA values when calculating frequencies
#'
#' @return A vector of allele frequencies for each marker
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2), 100*10, replace=TRUE), 100, 10)
#' freq <- calculate_allele_frequencies(M_012, coding="012")
#' }
#'
#' @export
calculate_allele_frequencies <- function(M, coding = NULL, na.rm = TRUE) {
  # If coding not provided, try to identify it
  if (is.null(coding)) {
    coding_info <- identify_marker_coding(M)
    coding <- coding_info$coding_type
  }
  
  n_markers <- ncol(M)
  n_samples <- nrow(M)
  
  # Initialize frequencies vector
  frequencies <- numeric(n_markers)
  
  # Calculate based on coding type
  if (coding == "01") {
    # For binary markers, frequency is simply the mean
    frequencies <- colMeans(M, na.rm = na.rm)
  } else if (coding == "012") {
    # For 0,1,2 coding, divide by 2 to get allele frequency
    frequencies <- colMeans(M, na.rm = na.rm) / 2
  } else if (coding == "-101") {
    # For -1,0,1 coding, add 1, then divide by 2
    frequencies <- colMeans(M + 1, na.rm = na.rm) / 2
  } else if (grepl("^0-[0-9]+$", coding)) {
    # For 0-k coding, divide by ploidy (k)
    ploidy <- as.numeric(sub("0-", "", coding))
    frequencies <- colMeans(M, na.rm = na.rm) / ploidy
  } else {
    stop("Unsupported coding type: ", coding)
  }
  
  return(frequencies)
}

#' Calculate minor allele frequencies and return marker statistics
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#' @param na.rm Logical, whether to remove NA values when calculating frequencies
#'
#' @return A data frame with marker statistics:
#'   \item{marker}{Marker name or index}
#'   \item{freq}{Allele frequency}
#'   \item{maf}{Minor allele frequency}
#'   \item{missing_rate}{Proportion of missing values}
#'   \item{heterozygosity}{Heterozygosity rate (for diploid markers)}
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2), 100*10, replace=TRUE), 100, 10)
#' stats <- marker_statistics(M_012, coding="012")
#' }
#'
#' @export
marker_statistics <- function(M, coding = NULL, na.rm = TRUE) {
  # If coding not provided, try to identify it
  if (is.null(coding)) {
    coding_info <- identify_marker_coding(M)
    coding <- coding_info$coding_type
  }
  
  n_markers <- ncol(M)
  n_samples <- nrow(M)
  
  # Get marker names if available
  marker_names <- colnames(M)
  if (is.null(marker_names)) {
    marker_names <- paste0("Marker", 1:n_markers)
  }
  
  # Calculate allele frequencies
  frequencies <- calculate_allele_frequencies(M, coding, na.rm)
  
  # Calculate minor allele frequencies
  mafs <- pmin(frequencies, 1 - frequencies)
  
  # Calculate missing rates
  missing_rates <- colMeans(is.na(M))
  
  # Calculate heterozygosity (only relevant for diploid markers)
  heterozygosity <- rep(NA, n_markers)
  if (coding == "012") {
    heterozygosity <- colMeans(M == 1, na.rm = na.rm)
  } else if (coding == "-101") {
    heterozygosity <- colMeans(M == 0, na.rm = na.rm)
  }
  
  # Compile results
  stats <- data.frame(
    marker = marker_names,
    freq = frequencies,
    maf = mafs,
    missing_rate = missing_rates,
    heterozygosity = heterozygosity
  )
  
  return(stats)
}

#' Filter markers based on various criteria
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param min_maf Numeric, minimum minor allele frequency (default: 0.01)
#' @param max_missing Numeric, maximum proportion of missing values (default: 0.2)
#' @param min_samples Integer, minimum number of samples with non-missing values
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#'
#' @return A filtered marker matrix
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2, NA), 100*1000, replace=TRUE), 100, 1000)
#' M_filtered <- filter_markers(M_012, min_maf=0.05, max_missing=0.1)
#' }
#'
#' @export
filter_markers <- function(M, min_maf = 0.01, max_missing = 0.2, 
                           min_samples = NULL, coding = NULL) {
  # Calculate marker statistics
  stats <- marker_statistics(M, coding)
  
  # Define filter criteria
  keep <- rep(TRUE, ncol(M))
  
  # MAF filter
  if (!is.null(min_maf)) {
    keep <- keep & (stats$maf >= min_maf)
  }
  
  # Missing data filter
  if (!is.null(max_missing)) {
    keep <- keep & (stats$missing_rate <= max_missing)
  }
  
  # Minimum samples filter
  if (!is.null(min_samples)) {
    n_present <- nrow(M) * (1 - stats$missing_rate)
    keep <- keep & (n_present >= min_samples)
  }
  
  # Apply filter
  M_filtered <- M[, keep, drop = FALSE]
  
  # Print summary
  n_removed <- sum(!keep)
  cat(sprintf("Filtered %d/%d markers (%.1f%%)\n", 
              n_removed, ncol(M), 100 * n_removed / ncol(M)))
  cat(sprintf("- MAF < %g: %d markers\n", 
              min_maf, sum(stats$maf < min_maf)))
  cat(sprintf("- Missing rate > %g: %d markers\n", 
              max_missing, sum(stats$missing_rate > max_missing)))
  if (!is.null(min_samples)) {
    n_present <- nrow(M) * (1 - stats$missing_rate)
    cat(sprintf("- Samples < %d: %d markers\n", 
                min_samples, sum(n_present < min_samples)))
  }
  
  return(M_filtered)
}

#' Impute missing marker data
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param method Character, imputation method:
#'   - "mean": Impute with the mean (most common for quantitative data)
#'   - "mode": Impute with the most frequent value (most common for discrete data)
#'   - "median": Impute with the median value
#'   - "knn": K-nearest neighbor imputation
#'   - "random": Random imputation from empirical distribution
#' @param k Integer, number of neighbors for KNN imputation (default: 5)
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#'
#' @return A marker matrix with imputed values
#'
#' @examples
#' \dontrun{
#' M <- matrix(sample(c(0,1,2, NA), 100*10, replace=TRUE), 100, 10)
#' M_imp <- impute_markers(M, method="mean", coding="012")
#' }
#'
#' @export
impute_markers <- function(M, method = "mean", k = 5, coding = NULL) {
  # If coding not provided, try to identify it
  if (is.null(coding)) {
    coding_info <- identify_marker_coding(M)
    coding <- coding_info$coding_type
  }
  
  # Initialize result
  M_imp <- M
  
  # Check if there are any missing values
  if (!any(is.na(M))) {
    message("No missing values to impute")
    return(M)
  }
  
  # Apply imputation method
  if (method == "mean") {
    # Impute with column means
    col_means <- colMeans(M, na.rm = TRUE)
    
    # For specified coding types, round to nearest valid value
    if (coding %in% c("01", "012", "-101")) {
      for (j in 1:ncol(M)) {
        na_idx <- is.na(M[, j])
        if (any(na_idx)) {
          if (coding == "01") {
            M_imp[na_idx, j] <- round(col_means[j])
          } else if (coding == "012") {
            M_imp[na_idx, j] <- round(col_means[j])
          } else if (coding == "-101") {
            # For -1,0,1 coding, round to nearest of -1, 0, 1
            M_imp[na_idx, j] <- round(col_means[j])
          }
        }
      }
    } else {
      # For continuous data, use means directly
      for (j in 1:ncol(M)) {
        na_idx <- is.na(M[, j])
        if (any(na_idx)) {
          M_imp[na_idx, j] <- col_means[j]
        }
      }
    }
  } else if (method == "mode") {
    # Impute with most frequent value in each column
    for (j in 1:ncol(M)) {
      na_idx <- is.na(M[, j])
      if (any(na_idx)) {
        # Get non-NA values
        values <- M[!na_idx, j]
        # Find mode (most frequent value)
        tab <- table(values)
        mode_val <- as.numeric(names(tab)[which.max(tab)])
        # Impute
        M_imp[na_idx, j] <- mode_val
      }
    }
  } else if (method == "median") {
    # Impute with column medians
    for (j in 1:ncol(M)) {
      na_idx <- is.na(M[, j])
      if (any(na_idx)) {
        M_imp[na_idx, j] <- median(M[, j], na.rm = TRUE)
      }
    }
  } else if (method == "knn") {
    # K-nearest neighbor imputation
    if (requireNamespace("impute", quietly = TRUE)) {
      M_t <- t(M)
      M_t_imp <- impute::impute.knn(M_t, k = k)$data
      M_imp <- t(M_t_imp)
    } else {
      message("The 'impute' Bioconductor package is not available. ",
              "Falling back to mean imputation. ",
              "Install with: BiocManager::install('impute')")
      M_imp <- impute_markers(M, method = "mean", coding = coding)
    }
  } else if (method == "random") {
    # Random imputation from empirical distribution
    for (j in 1:ncol(M)) {
      na_idx <- is.na(M[, j])
      if (any(na_idx)) {
        # Get non-NA values
        values <- M[!na_idx, j]
        # Sample randomly from observed values
        M_imp[na_idx, j] <- sample(values, sum(na_idx), replace = TRUE)
      }
    }
  } else {
    stop("Unsupported imputation method: ", method)
  }
  
  # For integer coding, ensure values are integers
  if (coding %in% c("01", "012", "-101")) {
    M_imp <- round(M_imp)
  }
  
  # Print summary
  n_imputed <- sum(is.na(M))
  cat(sprintf("Imputed %d missing values (%.1f%%) using method '%s'\n",
              n_imputed, 100 * n_imputed / (nrow(M) * ncol(M)), method))
  
  return(M_imp)
}

#' Calculate a genomic relationship matrix (GRM/kinship matrix)
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param method Character, method for calculating kinship:
#'   - "vanraden1": First method from VanRaden (2008)
#'   - "vanraden2": Second method from VanRaden (2008)
#'   - "astle": Astle & Balding method
#'   - "ibs": Identity by state (IBS)
#'   - "ibd": Identity by descent (IBD)
#'   - "centered": Simple centered cross-product
#'   - "correlation": Correlation-based
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#' @param scale Logical, whether to scale the kinship matrix (default: TRUE)
#' @param na.rm Logical, whether to remove missing values (default: TRUE)
#' @param na.method Character, how to handle missing data ("omit", "mean", "zero")
#' @param return_scaled_M Logical, whether to also return the scaled marker matrix
#'
#' @return A kinship matrix or a list with the kinship matrix and scaled markers
#'
#' @references
#' VanRaden, P. M. (2008). Efficient methods to compute genomic predictions.
#' Journal of Dairy Science, 91(11), 4414-4423.
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2), 100*1000, replace=TRUE), 100, 1000)
#' K <- calculate_kinship(M_012, method="vanraden1", coding="012")
#' }
#'
#' @export
calculate_kinship <- function(M, method = "vanraden1", coding = NULL, 
                              scale = TRUE, na.rm = TRUE, 
                              na.method = "mean", return_scaled_M = FALSE) {
  # If coding not provided, try to identify it
  if (is.null(coding)) {
    coding_info <- identify_marker_coding(M)
    coding <- coding_info$coding_type
  }
  
  n_samples <- nrow(M)
  n_markers <- ncol(M)
  
  # Handle missing values
  if (any(is.na(M))) {
    if (na.method == "omit") {
      # Only use markers with no missing values
      complete_cols <- colSums(is.na(M)) == 0
      if (sum(complete_cols) == 0) {
        stop("No complete markers available with na.method='omit'")
      }
      M <- M[, complete_cols, drop = FALSE]
      cat(sprintf("Using %d/%d complete markers\n", 
                  sum(complete_cols), n_markers))
    } else if (na.method == "mean") {
      # Impute missing values with column means
      M <- impute_markers(M, method = "mean", coding = coding)
    } else if (na.method == "zero") {
      # Replace missing values with zero
      M[is.na(M)] <- 0
    } else {
      stop("Unsupported na.method: ", na.method)
    }
  }
  
  # Calculate allele frequencies
  p <- calculate_allele_frequencies(M, coding, na.rm)
  
  # Initialize scaled marker matrix
  M_scaled <- M
  
  # Apply method-specific scaling
  if (method == "vanraden1") {
    # VanRaden's first method (2008)
    if (coding == "012") {
      # Center the marker matrix
      for (j in 1:ncol(M)) {
        M_scaled[, j] <- M[, j] - 2 * p[j]
      }
      
      # Calculate kinship
      K <- tcrossprod(M_scaled) / (2 * sum(p * (1 - p)))
    } else if (coding == "-101") {
      # For centered markers, adjust scaling
      for (j in 1:ncol(M)) {
        M_scaled[, j] <- M[, j] - 0 # Already centered
      }
      
      # Calculate kinship
      K <- tcrossprod(M_scaled) / (2 * sum(p * (1 - p)))
    } else {
      stop("VanRaden method requires '012' or '-101' coding")
    }
  } else if (method == "vanraden2") {
    # VanRaden's second method (2008)
    if (coding == "012") {
      # Scale the marker matrix
      for (j in 1:ncol(M)) {
        M_scaled[, j] <- (M[, j] - 2 * p[j]) / sqrt(2 * p[j] * (1 - p[j]))
      }
      
      # Calculate kinship
      K <- tcrossprod(M_scaled) / n_markers
    } else if (coding == "-101") {
      # For centered markers, adjust scaling
      for (j in 1:ncol(M)) {
        M_scaled[, j] <- M[, j] / sqrt(2 * p[j] * (1 - p[j]))
      }
      
      # Calculate kinship
      K <- tcrossprod(M_scaled) / n_markers
    } else {
      stop("VanRaden method requires '012' or '-101' coding")
    }
  } else if (method == "astle") {
    # Astle & Balding method
    if (coding == "012") {
      # Scale the marker matrix
      for (j in 1:ncol(M)) {
        M_scaled[, j] <- (M[, j] - 2 * p[j]) / sqrt(2 * p[j] * (1 - p[j]))
      }
      
      # Calculate kinship
      K <- tcrossprod(M_scaled) / n_markers
    } else if (coding == "-101") {
      # For centered markers, adjust scaling
      for (j in 1:ncol(M)) {
        M_scaled[, j] <- M[, j] / sqrt(2 * p[j] * (1 - p[j]))
      }
      
      # Calculate kinship
      K <- tcrossprod(M_scaled) / n_markers
    } else {
      stop("Astle & Balding method requires '012' or '-101' coding")
    }
  } else if (method == "ibs") {
    # Identity by state
    if (coding != "012") {
      stop("IBS method requires '012' coding")
    }
    
    # Calculate IBS matrix directly
    K <- matrix(0, n_samples, n_samples)
    for (i in 1:(n_samples-1)) {
      for (j in (i+1):n_samples) {
        # Count matching alleles
        matches <- sum(M[i, ] == M[j, ], na.rm = na.rm)
        total <- sum(!is.na(M[i, ]) & !is.na(M[j, ]))
        K[i, j] <- K[j, i] <- matches / (2 * total)
      }
    }
    
    # Set diagonal to 1
    diag(K) <- 1
  } else if (method == "ibd") {
    # Identity by descent approximation
    if (coding != "012") {
      stop("IBD method requires '012' coding")
    }
    
    # Calculate IBD matrix
    K <- matrix(0, n_samples, n_samples)
    for (i in 1:(n_samples-1)) {
      for (j in (i+1):n_samples) {
        # Count identical genotypes
        identical <- sum(M[i, ] == M[j, ] & !is.na(M[i, ]) & !is.na(M[j, ]))
        total <- sum(!is.na(M[i, ]) & !is.na(M[j, ]))
        K[i, j] <- K[j, i] <- identical / total
      }
    }
    
    # Set diagonal to 1
    diag(K) <- 1
  } else if (method == "centered") {
    # Simple centered cross-product
    
    # Center the marker matrix
    for (j in 1:ncol(M)) {
      col_mean <- mean(M[, j], na.rm = na.rm)
      M_scaled[, j] <- M[, j] - col_mean
    }
    
    # Calculate kinship
    K <- tcrossprod(M_scaled) / n_markers
  } else if (method == "correlation") {
    # Correlation-based
    
    # Standardize the marker matrix
    for (j in 1:ncol(M)) {
      col_mean <- mean(M[, j], na.rm = na.rm)
      col_sd <- sd(M[, j], na.rm = na.rm)
      if (col_sd > 0) {
        M_scaled[, j] <- (M[, j] - col_mean) / col_sd
      } else {
        M_scaled[, j] <- 0
      }
    }
    
    # Calculate kinship
    K <- tcrossprod(M_scaled) / n_markers
  } else {
    stop("Unsupported kinship method: ", method)
  }
  
  # Scale the kinship matrix to have mean diagonal = 1
  if (scale) {
    K <- K / mean(diag(K))
  }
  
  # Check and ensure matrix is positive definite
  if (any(is.na(K)) || any(!is.finite(K))) {
    warning("Kinship matrix contains NA or infinite values")
  }
  
  # Return result
  if (return_scaled_M) {
    return(list(K = K, M_scaled = M_scaled))
  } else {
    return(K)
  }
}

#' Visualize a kinship matrix as a heatmap
#'
#' @param K A kinship (relationship) matrix
#' @param labels Character vector of sample labels (optional)
#' @param clustering Logical, whether to perform hierarchical clustering (default: TRUE)
#' @param cluster_method Character, clustering method (default: "complete")
#' @param dist_method Character, distance method for clustering (default: "euclidean")
#' @param palette Character, color palette to use (default: "viridis")
#' @param title Character, plot title (default: "Kinship Matrix")
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2), 100*1000, replace=TRUE), 100, 1000)
#' K <- calculate_kinship(M_012, method="vanraden1", coding="012")
#' visualize_kinship(K)
#' }
#'
#' @export
visualize_kinship <- function(K, labels = NULL, clustering = TRUE,
                             cluster_method = "complete", dist_method = "euclidean",
                             palette = "viridis", title = "Kinship Matrix") {
  # Check for ggplot2 package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required for visualization. ",
         "Please install it with install.packages('ggplot2')")
  }
  
  # Get dimensions
  n <- nrow(K)
  
  # Set default labels if not provided
  if (is.null(labels)) {
    if (!is.null(rownames(K))) {
      labels <- rownames(K)
    } else {
      labels <- 1:n
    }
  }
  
  # Apply clustering if requested
  if (clustering) {
    # Compute distance matrix
    d <- stats::dist(K, method = dist_method)
    
    # Perform hierarchical clustering
    hc <- stats::hclust(d, method = cluster_method)
    
    # Get the clustering order
    ordering <- hc$order
  } else {
    # Use original order
    ordering <- 1:n
  }
  
  # Reorder the matrix
  K_ordered <- K[ordering, ordering]
  labels_ordered <- labels[ordering]
  
  # Convert to long format for ggplot
  K_long <- data.frame(
    row = rep(1:n, n),
    col = rep(1:n, each = n),
    value = as.vector(K_ordered)
  )
  
  # Create plot with ggplot2
  p <- ggplot2::ggplot(K_long, ggplot2::aes(x = col, y = n - row + 1, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_continuous(breaks = 1:n, labels = labels_ordered, 
                               expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = 1:n, labels = rev(labels_ordered), 
                               expand = c(0, 0)) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "", y = "", fill = "Value") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  # Apply color palette
  if (palette == "viridis") {
    if (requireNamespace("viridis", quietly = TRUE)) {
      p <- p + viridis::scale_fill_viridis()
    } else {
      p <- p + ggplot2::scale_fill_gradient2(low = "blue", mid = "white", 
                                            high = "red", midpoint = 0.5)
    }
  } else if (palette == "heat") {
    p <- p + ggplot2::scale_fill_gradient2(low = "white", mid = "orange", 
                                          high = "red", midpoint = 0.5)
  } else {
    p <- p + ggplot2::scale_fill_gradient2(low = "blue", mid = "white", 
                                          high = "red", midpoint = 0.5)
  }
  
  return(p)
}

#' Plot marker data summary statistics
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#' @param plot_type Character, type of plot to create ("maf", "missing", "both")
#'
#' @return A ggplot2 object or list of plots
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2, NA), 100*1000, replace=TRUE, prob=c(0.3, 0.4, 0.2, 0.1)), 100, 1000)
#' p <- plot_marker_stats(M_012, coding="012", plot_type="both")
#' print(p$maf)
#' print(p$missing)
#' }
#'
#' @export
plot_marker_stats <- function(M, coding = NULL, plot_type = "both") {
  # Check for ggplot2 package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required for visualization. ",
         "Please install it with install.packages('ggplot2')")
  }
  
  # Calculate marker statistics
  stats <- marker_statistics(M, coding)
  
  # Prepare plots
  plots <- list()
  
  # MAF distribution
  if (plot_type %in% c("maf", "both")) {
    p_maf <- ggplot2::ggplot(stats, ggplot2::aes(x = maf)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
      ggplot2::labs(title = "Minor Allele Frequency Distribution",
                   x = "Minor Allele Frequency", y = "Count") +
      ggplot2::theme_minimal()
    
    plots$maf <- p_maf
  }
  
  # Missing data distribution
  if (plot_type %in% c("missing", "both")) {
    p_missing <- ggplot2::ggplot(stats, ggplot2::aes(x = missing_rate)) +
      ggplot2::geom_histogram(bins = 30, fill = "firebrick", color = "black", alpha = 0.7) +
      ggplot2::labs(title = "Missing Data Distribution",
                   x = "Missing Rate", y = "Count") +
      ggplot2::theme_minimal()
    
    plots$missing <- p_missing
  }
  
  # Return result
  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}

#' Create a summary table of marker data properties
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#'
#' @return A data.frame with summary statistics
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2, NA), 100*1000, replace=TRUE), 100, 1000)
#' summary_table <- summarize_marker_data(M_012, coding="012")
#' print(summary_table)
#' }
#'
#' @export
summarize_marker_data <- function(M, coding = NULL) {
  # If coding not provided, try to identify it
  if (is.null(coding)) {
    coding_info <- identify_marker_coding(M)
    coding <- coding_info$coding_type
  }
  
  # Calculate marker statistics
  stats <- marker_statistics(M, coding)
  
  # Create summary table
  summary <- data.frame(
    Metric = c(
      "Number of samples",
      "Number of markers",
      "Marker coding",
      "Missing rate (mean)",
      "Missing rate (min)",
      "Missing rate (max)",
      "MAF (mean)",
      "MAF (min)",
      "MAF (max)",
      "Heterozygosity (mean)",
      "Heterozygosity (min)",
      "Heterozygosity (max)"
    ),
    Value = c(
      nrow(M),
      ncol(M),
      coding,
      mean(stats$missing_rate),
      min(stats$missing_rate),
      max(stats$missing_rate),
      mean(stats$maf),
      min(stats$maf),
      max(stats$maf),
      mean(stats$heterozygosity, na.rm = TRUE),
      min(stats$heterozygosity, na.rm = TRUE),
      max(stats$heterozygosity, na.rm = TRUE)
    )
  )
  
  # Format numeric values
  for (i in 4:nrow(summary)) {
    summary$Value[i] <- sprintf("%.4f", as.numeric(summary$Value[i]))
  }
  
  return(summary)
}

#' Subset a marker matrix based on MAF buckets
#'
#' This function creates a stratified subset of markers based on MAF buckets,
#' which can be useful for testing different marker selection strategies.
#'
#' @param M A marker matrix with samples in rows and markers in columns
#' @param n_markers Number of markers to select in total
#' @param breaks Vector of MAF breakpoints for stratification (default: 5 equal buckets)
#' @param coding Character, the marker coding type ("01", "012", "-101", or "custom")
#'
#' @return A subset of the marker matrix with selected markers
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2), 100*1000, replace=TRUE), 100, 1000)
#' M_sub <- subset_markers_by_maf(M_012, n_markers=100, coding="012")
#' }
#'
#' @export
subset_markers_by_maf <- function(M, n_markers, breaks = NULL, coding = NULL) {
  # Calculate marker statistics
  stats <- marker_statistics(M, coding)
  
  # Set default breaks if not provided
  if (is.null(breaks)) {
    breaks <- seq(0, 0.5, length.out = 6)
  }
  
  # Create MAF buckets
  maf_buckets <- cut(stats$maf, breaks = breaks, include.lowest = TRUE)
  bucket_counts <- table(maf_buckets)
  
  # Calculate number of markers to select from each bucket
  n_per_bucket <- round(n_markers * bucket_counts / sum(bucket_counts))
  
  # Ensure total matches n_markers
  diff <- n_markers - sum(n_per_bucket)
  if (diff > 0) {
    n_per_bucket[which.max(bucket_counts)] <- n_per_bucket[which.max(bucket_counts)] + diff
  } else if (diff < 0) {
    i <- which.max(n_per_bucket)
    n_per_bucket[i] <- n_per_bucket[i] + diff
  }
  
  # Select markers from each bucket
  selected_markers <- integer(0)
  for (i in 1:length(levels(maf_buckets))) {
    bucket <- levels(maf_buckets)[i]
    markers_in_bucket <- which(maf_buckets == bucket)
    if (length(markers_in_bucket) > 0) {
      n_select <- min(n_per_bucket[i], length(markers_in_bucket))
      if (n_select > 0) {
        sampled <- sample(markers_in_bucket, n_select)
        selected_markers <- c(selected_markers, sampled)
      }
    }
  }
  
  # Create subset matrix
  M_subset <- M[, selected_markers, drop = FALSE]
  
  # Report info
  cat(sprintf("Selected %d markers with stratified MAF distribution:\n", ncol(M_subset)))
  for (i in 1:length(levels(maf_buckets))) {
    bucket <- levels(maf_buckets)[i]
    n_selected <- sum(maf_buckets[selected_markers] == bucket)
    pct <- 100 * n_selected / ncol(M_subset)
    cat(sprintf("  %s: %d markers (%.1f%%)\n", bucket, n_selected, pct))
  }
  
  return(M_subset)
}

#' Compare different kinship matrices
#'
#' @param K_list List of kinship matrices
#' @param names Character vector of names for the matrices
#' @param plot Logical, whether to create a plot (default: TRUE)
#'
#' @return If plot=TRUE, a ggplot2 object, otherwise correlation matrix
#'
#' @examples
#' \dontrun{
#' M_012 <- matrix(sample(c(0,1,2), 100*1000, replace=TRUE), 100, 1000)
#' K1 <- calculate_kinship(M_012, method="vanraden1", coding="012")
#' K2 <- calculate_kinship(M_012, method="vanraden2", coding="012")
#' K3 <- calculate_kinship(M_012, method="correlation", coding="012")
#' compare_kinship(list(K1, K2, K3), names=c("VanRaden1", "VanRaden2", "Correlation"))
#' }
#'
#' @export
compare_kinship <- function(K_list, names = NULL, plot = TRUE) {
  # Check if all matrices have the same dimension
  dims <- unique(sapply(K_list, nrow))
  if (length(dims) > 1) {
    stop("All kinship matrices must have the same dimensions")
  }
  
  # Set default names if not provided
  if (is.null(names)) {
    names <- paste0("K", 1:length(K_list))
  }
  
  # Convert matrices to vectors for correlation analysis
  K_vecs <- lapply(K_list, function(K) K[upper.tri(K, diag = FALSE)])
  K_df <- as.data.frame(do.call(cbind, K_vecs))
  colnames(K_df) <- names
  
  # Calculate correlation matrix
  cor_mat <- cor(K_df)
  
  if (!plot) {
    return(cor_mat)
  }
  
  # Check for ggplot2 package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("The 'ggplot2' package is required for plotting. Returning correlation matrix instead.")
    return(cor_mat)
  }
  
  # Convert correlation matrix to long format for plotting
  cor_long <- data.frame(
    Method1 = rep(names, each = length(names)),
    Method2 = rep(names, length(names)),
    Correlation = as.vector(cor_mat)
  )
  
  # Create heatmap
  p <- ggplot2::ggplot(cor_long, ggplot2::aes(x = Method1, y = Method2, fill = Correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                 midpoint = 0.5, limits = c(0, 1)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", Correlation)),
                      color = "black", size = 3) +
    ggplot2::labs(title = "Correlation Between Kinship Matrices",
                 x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}
