# S3 Methods for EMMREML objects

#' Print Method for emmreml Objects
#'
#' @param x An object of class 'emmreml'
#' @param ... Additional arguments
#'
#' @return Invisibly returns the input object
#'
#' @method print emmreml
#' @export
print.emmreml <- function(x, ...) {
  cat("Mixed Model with Known Covariance Structure\n\n")
  
  # Fixed effects
  cat("Fixed effects estimates (BLUEs):\n")
  fe <- x$betahat
  if (!is.null(names(fe))) {
    fe_df <- data.frame(Estimate = fe)
  } else {
    fe_df <- data.frame(Estimate = fe)
    rownames(fe_df) <- paste0("beta", seq_along(fe))
  }
  print(fe_df)
  cat("\n")
  
  # Variance components
  cat("Variance components:\n")
  vc <- c(x$Vu, x$Ve)
  names(vc) <- c("Vu", "Ve")
  print(vc)
  cat("\n")
  
  # Model fit
  cat("Model fit:\n")
  cat("Log-likelihood:", x$loglik, "\n")
  if (!is.null(x$AIC)) cat("AIC:", x$AIC, "\n")
  if (!is.null(x$BIC)) cat("BIC:", x$BIC, "\n")
  
  invisible(x)
}

#' Summary Method for emmreml Objects
#'
#' @param object An object of class 'emmreml'
#' @param ... Additional arguments
#'
#' @return A summary object
#'
#' @method summary emmreml
#' @export
summary.emmreml <- function(object, ...) {
  # Create a summary object
  result <- list(
    call = object$call,
    betahat = object$betahat,
    Vu = object$Vu,
    Ve = object$Ve,
    loglik = object$loglik,
    AIC = object$AIC,
    BIC = object$BIC
  )
  
  # Add standard errors if available
  if (!is.null(object$varbetahat)) {
    result$se_beta <- sqrt(object$varbetahat)
    
    # Add t-values and p-values
    result$t_value <- object$betahat / result$se_beta
    result$p_value <- 2 * pt(abs(result$t_value), df = length(object$uhat) - length(object$betahat), lower.tail = FALSE)
  }
  
  class(result) <- "summary.emmreml"
  return(result)
}

#' Print Method for summary.emmreml Objects
#'
#' @param x An object of class 'summary.emmreml'
#' @param ... Additional arguments
#'
#' @return Invisibly returns the input object
#'
#' @method print summary.emmreml
#' @export
print.summary.emmreml <- function(x, ...) {
  cat("Mixed Model with Known Covariance Structure\n\n")
  
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  
  # Fixed effects table
  cat("Fixed effects:\n")
  if (!is.null(x$se_beta)) {
    fe_table <- data.frame(
      Estimate = x$betahat,
      StdError = x$se_beta,
      tValue = x$t_value,
      pValue = x$p_value
    )
    # Add significance stars
    stars <- symnum(x$p_value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
    fe_table$Significance <- stars
  } else {
    fe_table <- data.frame(Estimate = x$betahat)
  }
  if (is.null(names(x$betahat))) {
    rownames(fe_table) <- paste0("beta", seq_along(x$betahat))
  }
  print(fe_table)
  
  if (!is.null(x$se_beta)) {
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  }
  
  # Variance components
  cat("Variance components:\n")
  vc <- c(x$Vu, x$Ve)
  names(vc) <- c("Vu", "Ve")
  print(vc)
  cat("\n")
  
  # Model fit statistics
  cat("Model fit:\n")
  cat("Log-likelihood:", x$loglik, "\n")
  if (!is.null(x$AIC)) cat("AIC:", x$AIC, "\n")
  if (!is.null(x$BIC)) cat("BIC:", x$BIC, "\n")
  
  invisible(x)
}

#' Predict Method for emmreml Objects
#'
#' @param object An object of class 'emmreml'
#' @param newdata Optional data frame for predictions
#' @param type Type of prediction: "fixed" for fixed effects only, "random" for random effects only, or "total" for both
#' @param se.pred Logical, whether to return standard errors
#' @param ... Additional arguments
#'
#' @return A list with predictions
#'
#' @method predict emmreml
#' @export
predict.emmreml <- function(object, newdata = NULL, type = "total", se.pred = FALSE, ...) {
  if (!is.null(newdata)) {
    stop("Prediction with new data not yet implemented")
  }
  
  # Return different types of predictions based on 'type'
  result <- list()
  
  if (type == "fixed" || type == "total") {
    result$fixed <- object$betahat
  }
  
  if (type == "random" || type == "total") {
    result$random <- object$uhat
  }
  
  if (type == "total") {
    result$fitted <- result$fixed + result$random
  }
  
  # Add standard errors if requested
  if (se.pred) {
    if ((type == "fixed" || type == "total") && !is.null(object$varbetahat)) {
      result$se.fixed <- sqrt(object$varbetahat)
    }
    
    if ((type == "random" || type == "total") && !is.null(object$varuhat)) {
      result$se.random <- sqrt(object$varuhat)
    }
    
    if (type == "total" && !is.null(object$varbetahat) && !is.null(object$varuhat)) {
      # For total, we need to account for covariance between fixed and random effects
      # This is a simplification assuming independence
      result$se.fitted <- sqrt(result$se.fixed^2 + result$se.random^2)
    }
  }
  
  return(result)
}

#' Print Method for emmreml_multi Objects
#'
#' @param x An object of class 'emmreml_multi'
#' @param ... Additional arguments
#'
#' @return Invisibly returns the input object
#'
#' @method print emmreml_multi
#' @export
print.emmreml_multi <- function(x, ...) {
  cat("Multi-Kernel Mixed Model with Known Covariance Structures\n\n")
  
  # Fixed effects
  cat("Fixed effects estimates (BLUEs):\n")
  fe <- x$betahat
  if (!is.null(names(fe))) {
    fe_df <- data.frame(Estimate = fe)
  } else {
    fe_df <- data.frame(Estimate = fe)
    rownames(fe_df) <- paste0("beta", seq_along(fe))
  }
  print(fe_df)
  cat("\n")
  
  # Variance components
  cat("Variance components:\n")
  vc <- c(x$Vu, x$Ve)
  names(vc) <- c("Vu", "Ve")
  print(vc)
  cat("\n")
  
  # Kernel weights
  if (!is.null(x$weights)) {
    cat("Kernel weights:\n")
    w_df <- data.frame(Weight = x$weights)
    rownames(w_df) <- paste0("K", seq_along(x$weights))
    print(w_df)
    cat("\n")
  }
  
  # Model fit
  cat("Model fit:\n")
  cat("Log-likelihood:", x$loglik, "\n")
  if (!is.null(x$AIC)) cat("AIC:", x$AIC, "\n")
  if (!is.null(x$BIC)) cat("BIC:", x$BIC, "\n")
  
  invisible(x)
}

#' Summary Method for emmreml_multi Objects
#'
#' @param object An object of class 'emmreml_multi'
#' @param ... Additional arguments
#'
#' @return A summary object
#'
#' @method summary emmreml_multi
#' @export
summary.emmreml_multi <- function(object, ...) {
  # Call the base summary method
  result <- summary.emmreml(object, ...)
  
  # Add multi-kernel specific information
  result$weights <- object$weights
  
  class(result) <- c("summary.emmreml_multi", "summary.emmreml")
  return(result)
}

#' Print Method for summary.emmreml_multi Objects
#'
#' @param x An object of class 'summary.emmreml_multi'
#' @param ... Additional arguments
#'
#' @return Invisibly returns the input object
#'
#' @method print summary.emmreml_multi
#' @export
print.summary.emmreml_multi <- function(x, ...) {
  cat("Multi-Kernel Mixed Model with Known Covariance Structures\n\n")
  
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  
  # Fixed effects table
  cat("Fixed effects:\n")
  if (!is.null(x$se_beta)) {
    fe_table <- data.frame(
      Estimate = x$betahat,
      StdError = x$se_beta,
      tValue = x$t_value,
      pValue = x$p_value
    )
    # Add significance stars
    stars <- symnum(x$p_value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
    fe_table$Significance <- stars
  } else {
    fe_table <- data.frame(Estimate = x$betahat)
  }
  if (is.null(names(x$betahat))) {
    rownames(fe_table) <- paste0("beta", seq_along(x$betahat))
  }
  print(fe_table)
  
  if (!is.null(x$se_beta)) {
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  }
  
  # Variance components
  cat("Variance components:\n")
  vc <- c(x$Vu, x$Ve)
  names(vc) <- c("Vu", "Ve")
  print(vc)
  cat("\n")
  
  # Kernel weights
  if (!is.null(x$weights)) {
    cat("Kernel weights:\n")
    w_df <- data.frame(Weight = x$weights)
    rownames(w_df) <- paste0("K", seq_along(x$weights))
    print(w_df)
    cat("\n")
  }
  
  # Model fit statistics
  cat("Model fit:\n")
  cat("Log-likelihood:", x$loglik, "\n")
  if (!is.null(x$AIC)) cat("AIC:", x$AIC, "\n")
  if (!is.null(x$BIC)) cat("BIC:", x$BIC, "\n")
  
  invisible(x)
}
