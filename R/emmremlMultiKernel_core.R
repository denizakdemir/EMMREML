# Continuation of the Multi-Kernel Core function

#' Core Solver for Multi-Kernel EMMA Algorithm
#'
#' @description
#' Core function that implements the multi-kernel EMMA algorithm.
#' This is an internal function called by the main emmremlMultiKernel function.
#'
#' @param y Response vector
#' @param X Fixed effects design matrix
#' @param Zlist List of random effects design matrices
#' @param Klist List of covariance matrices for random effects
#' @param varbetahat Whether to calculate variance of beta estimates
#' @param varuhat Whether to calculate variance of BLUPs
#' @param PEVuhat Whether to calculate prediction error variance
#' @param test Whether to perform significance tests
#'
#' @return A list containing model results
#'
#' @keywords internal
emmremlMultiKernel_core <- function(y, X, Zlist, Klist, varbetahat = FALSE, 
                                    varuhat = FALSE, PEVuhat = FALSE, test = FALSE) {
  q <- dim(X)[2]
  n <- length(y)
  lz <- length(Zlist)
  
  spI <- diag(n)
  S <- spI - tcrossprod(X %*% solve(crossprod(X)), X)
  
  K_is_identity <- vapply(seq_len(lz), function(i) is_identity_matrix(Klist[[i]]), logical(1))
  
  Z <- c()
  for (i in 1:lz) {
    Z <- cbind(Z, Zlist[[i]])
  }
  
  # Define the function to minimize
  minimfunctionouter <- function(weights = rep(1/lz, lz)) {
    weights <- weights / sum(weights)
    ZKZt <- matrix(0, nrow = n, ncol = n)
    for (i in 1:lz) {
      if (K_is_identity[i]) {
        # Simplified computation when K is identity
        ZKZt <- ZKZt + weights[i] * tcrossprod(Zlist[[i]])
      } else {
        ZKZt <- ZKZt + weights[i] * tcrossprod(Zlist[[i]] %*% Klist[[i]], Zlist[[i]])
      }
    }
    
    offset <- log(n)
    ZKZtandoffset <- ZKZt + offset * spI
    SZKZtSandoffset <- {S %*% ZKZtandoffset} %*% S
    
    svdSZKZtSandspI <- eigen(SZKZtSandoffset, symmetric = TRUE)
    Ur <- svdSZKZtSandspI$vectors[, 1:(n-q)]
    
    lambda <- svdSZKZtSandspI$values[1:(n-q)] - offset
    eta <- crossprod(Ur, y)
    minimfunc <- function(delta) {
      (n-q) * log(sum(eta^2 / {lambda + delta})) + sum(log(lambda + delta))
    }
    optimout <- optimize(minimfunc, lower = 0, upper = 10000, tol = 1e-9)
    return(optimout$objective)
  }
  
  # Optimize to find weights
  weights <- optim(
    par = rep(1/lz, lz), 
    fn = minimfunctionouter,
    method = "L-BFGS-B", 
    lower = rep(0, lz), 
    upper = rep(1, lz)
  )$par
  
  # Normalize weights
  weights <- weights / sum(weights)
  
  # Recombine Z matrices
  Z <- c()
  for (i in 1:lz) {
    Z <- cbind(Z, Zlist[[i]])
  }
  
  # Combine weighted K matrices
  ZKZt <- matrix(0, nrow = n, ncol = n)
  Klistweighted <- Klist
  for (i in 1:lz) {
    Klistweighted[[i]] <- weights[i] * Klist[[i]]
    
    if (K_is_identity[i]) {
      # Simplified computation when K is identity
      ZKZt <- ZKZt + weights[i] * tcrossprod(Zlist[[i]])
    } else {
      ZKZt <- ZKZt + weights[i] * tcrossprod(Zlist[[i]] %*% Klist[[i]], Zlist[[i]])
    }
  }
  
  # Create block diagonal K matrix
  K <- Matrix::.bdiag(Klistweighted)
  ZK <- as.matrix(Z %*% K)
  
  # Continue with EMMA algorithm
  offset <- log(n)
  ZKZtandoffset <- ZKZt + offset * spI
  SZKZtSandoffset <- {S %*% ZKZtandoffset} %*% S
  
  svdSZKZtSandspI <- eigen(SZKZtSandoffset, symmetric = TRUE)
  Ur <- svdSZKZtSandspI$vectors[, 1:(n-q)]
  
  lambda <- svdSZKZtSandspI$values[1:(n-q)] - offset
  eta <- crossprod(Ur, y)
  minimfunc <- function(delta) {
    (n-q) * log(sum(eta^2 / {lambda + delta})) + sum(log(lambda + delta))
  }
  optimout <- optimize(minimfunc, lower = 0, upper = 10000, tol = 1e-9)
  deltahat <- optimout$minimum
  
  # Calculate estimates and predictions
  Hinvhat <- solve(ZKZt + deltahat * spI)
  XtHinvhat <- crossprod(X, Hinvhat)
  betahat <- solve(XtHinvhat %*% X, XtHinvhat %*% y)
  ehat <- (y - {X %*% betahat})
  Hinvhatehat <- Hinvhat %*% ehat
  
  sigmausqhat <- sum(eta^2 / {lambda + deltahat}) / (n-q)
  sigmaesqhat <- deltahat * sigmausqhat
  uhat <- crossprod(ZK, Hinvhatehat)
  Vinv <- (1/sigmausqhat) * Hinvhat
  
  # Create names for random effects
  namesuhat <- c()
  for (i in 1:length(Klist)) {
    if (!is.null(colnames(Klist[[i]]))) {
      namesuhat <- c(namesuhat, paste(paste("K", i, sep = "_"), colnames(Klist[[i]]), sep = ""))
    } else {
      namesuhat <- c(namesuhat, paste("K", i, "_", 1:ncol(Klist[[i]]), sep = ""))
    }
  }
  
  # Calculate log-likelihood
  df <- n - q
  loglik <- -0.5 * (optimout$objective + df + df * log(2 * pi/df))
  
  # Initialize empty variables for returns
  varuhat_diag <- NULL
  varbetahat_diag <- NULL
  PEVuhat_diag <- NULL
  Xsqtestbeta <- NULL
  Xsqtestu <- NULL
  p_adjbetahat <- NULL
  p_adjuhat <- NULL
  
  # Calculate variances and tests if requested
  if (varuhat) {
    P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
    varuhat_mat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK
    varuhat_diag <- diag(varuhat_mat)
  }
  
  # Calculate prediction error variance if requested
  if (PEVuhat) {
    if (!exists("P")) {
      P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
    }
    if (!exists("varuhat_mat")) {
      varuhat_mat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK
    }
    PEVuhat_mat <- sigmausqhat * K - varuhat_mat
    PEVuhat_diag <- diag(PEVuhat_mat)
  }
  
  # Calculate variance of fixed effects if requested
  if (varbetahat) {
    varbetahat_mat <- solve(crossprod(X, Vinv %*% X))
    varbetahat_diag <- diag(varbetahat_mat)
  }
  
  # Perform significance tests if requested
  if (test) {
    if (!exists("varuhat_mat")) {
      P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
      varuhat_mat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK
      varuhat_diag <- diag(varuhat_mat)
    }
    if (!exists("varbetahat_mat")) {
      varbetahat_mat <- solve(crossprod(X, Vinv %*% X))
      varbetahat_diag <- diag(varbetahat_mat)
    }
    
    # Test for random effects
    Xsqtestu <- uhat^2 / varuhat_diag
    puhat <- pchisq(Xsqtestu, df = 1, lower.tail = FALSE, log.p = FALSE)
    p.adjust.M <- p.adjust.methods
    p_adjuhat <- sapply(p.adjust.M, function(meth) p.adjust(puhat, meth))
    
    # Test for fixed effects
    Xsqtestbeta <- betahat^2 / varbetahat_diag
    pbetahat <- pchisq(Xsqtestbeta, df = 1, lower.tail = FALSE, log.p = FALSE)
    p_adjbetahat <- sapply(p.adjust.M, function(meth) p.adjust(pbetahat, meth))
  }
  
  # Set names for the random effects and return results
  names(uhat) <- namesuhat
  
  return(list(
    Vu = sigmausqhat,
    Ve = sigmaesqhat,
    betahat = betahat,
    uhat = uhat,
    weights = weights,
    Xsqtestbeta = Xsqtestbeta,
    pvalbeta = p_adjbetahat,
    Xsqtestu = Xsqtestu,
    pvalu = p_adjuhat,
    varuhat = varuhat_diag,
    varbetahat = varbetahat_diag,
    PEVuhat = PEVuhat_diag,
    loglik = loglik
  ))
}
