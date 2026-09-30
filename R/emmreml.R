emmreml <- function(y, X, Z, K, varbetahat = FALSE, varuhat = FALSE, PEVuhat = FALSE, test = FALSE) {
  q <- dim(X)[2]
  n <- length(y)

  K_is_identity <- is_identity_matrix(K)

  spI <- diag(n)
  S <- spI - tcrossprod(X %*% solve(crossprod(X)), X)

  if (K_is_identity) {
    ZK <- Z
    ZKZt <- tcrossprod(Z)
  } else {
    ZK <- Z %*% K
    ZKZt <- tcrossprod(ZK, Z)
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
  
  optimout <- optimize(minimfunc, lower = 9^(-9), upper = 9^9, tol = .000001)
  deltahat <- optimout$minimum
 
  Hinvhat <- solve(ZKZt + deltahat * spI)
 
  XtHinvhat <- crossprod(X, Hinvhat)
  betahat <- solve(XtHinvhat %*% X, XtHinvhat %*% y)
  ehat <- (y - {X %*% betahat})
  Hinvhatehat <- Hinvhat %*% ehat
  
  sigmausqhat <- sum(eta^2 / {lambda + deltahat}) / (n-q)
  Vinv <- (1/sigmausqhat) * Hinvhat
  sigmaesqhat <- deltahat * sigmausqhat
  
  uhat <- crossprod(ZK, Hinvhatehat)
  
  df <- n - q
  loglik <- -0.5 * (optimout$objective + df + df * log(2 * pi/df))
  
  # Initialize empty variables
  varuhat_diag <- c()
  varbetahat_diag <- c()
  PEVuhat_diag <- c()
  
  ####VAR U
  if (varuhat) {
    P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
    varuhat_mat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK
    varuhat_diag <- diag(varuhat_mat)
  }
  
  if (PEVuhat) {
    if (!exists("P")) {
      P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
    }
    if (!exists("varuhat_mat")) {
      varuhat_mat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK
    }
    
    # Simplify when K is identity
    if(K_is_identity) {
      PEVuhat_mat <- sigmausqhat * diag(ncol(Z)) - varuhat_mat
    } else {
      PEVuhat_mat <- sigmausqhat * K - varuhat_mat
    }
    
    PEVuhat_diag <- diag(PEVuhat_mat)
  }

  # Var beta
  if (varbetahat) {
    varbetahat_mat <- solve(crossprod(X, Vinv %*% X))
    varbetahat_diag <- diag(varbetahat_mat)
  }
  
  # Test statistics
  Xsqtestbeta <- c()
  Xsqtestu <- c()
  p_adjuhat <- c()
  p_adjbetahat <- c()
  
  if (test) {
    if (length(varuhat_diag) == 0) {
      P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
      varuhat_mat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK
      varuhat_diag <- diag(varuhat_mat)
    }
    
    if (length(varbetahat_diag) == 0) {
      varbetahat_mat <- solve(crossprod(X, Vinv %*% X))
      varbetahat_diag <- diag(varbetahat_mat)
    }
    
    Xsqtestu <- uhat^2 / varuhat_diag
    puhat <- pchisq(Xsqtestu, df = 1, lower.tail = FALSE, log.p = FALSE)
    p.adjust.M <- p.adjust.methods
    p_adjuhat <- sapply(p.adjust.M, function(meth) p.adjust(puhat, meth))
    
    Xsqtestbeta <- betahat^2 / varbetahat_diag
    pbetahat <- pchisq(Xsqtestbeta, df = 1, lower.tail = FALSE, log.p = FALSE)
    p_adjbetahat <- sapply(p.adjust.M, function(meth) p.adjust(pbetahat, meth))
  }

  return(list(
    Vu = sigmausqhat,
    Ve = sigmaesqhat,
    betahat = betahat,
    uhat = uhat, 
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
