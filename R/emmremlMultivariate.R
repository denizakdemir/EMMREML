#' Fit a Multivariate Gaussian Mixed Model with Known Covariance Structure
#'
#' @description
#' Fits a multivariate Gaussian mixed model of the form
#' \eqn{Y = XB + ZG + E}, where \eqn{vec(G) \sim N(0, K \otimes V_g)} and
#' \eqn{vec(E) \sim N(0, I_n \otimes V_e)}, using the ECM algorithm.
#'
#' This is the primary (matrix-based) interface, preserved for backward
#' compatibility with v3.1. The formula-based wrapper \code{emmremlMultivariate}
#' calls this function internally.
#'
#' @param Y Response matrix (d x n): d traits, n individuals. NOTE: traits in
#'   rows, individuals in columns — same convention as v3.1.
#' @param X Fixed effects design matrix (m x n): m effects, n individuals.
#' @param Z Random effects incidence matrix (n x N): n individuals, N levels.
#' @param K Covariance/kernel matrix (N x N) for random effects.
#' @param varBhat Logical, whether to return variance of B estimates.
#' @param varGhat Logical, whether to return variance of G predictions.
#' @param PEVGhat Logical, whether to return prediction error variance of G.
#' @param test Logical, whether to perform chi-square significance tests.
#' @param tolpar Tolerance for convergence of variance components.
#' @param tolparinv Regularisation added to matrices before inversion.
#' @param maxIter Maximum number of ECM iterations.
#'
#' @return A list (class \code{emmreml_mv}) containing:
#'   \item{Bhat}{Fixed effects estimates (m x n eigenvector-space; same as v3.1 Bt)}
#'   \item{Vg}{Estimated genetic covariance matrix (d x d)}
#'   \item{Ve}{Estimated residual covariance matrix (d x d)}
#'   \item{Gpred}{BLUPs for G (d x N matrix, same as v3.1)}
#'   \item{iterations}{Number of ECM iterations run}
#'   \item{convergence}{Logical: did the algorithm converge?}
#'   \item{varBhat}{(optional) variance of vec(B)}
#'   \item{varGhat}{(optional) variance of vec(G)}
#'   \item{PEVGhat}{(optional) prediction error variance of vec(G)}
#'   \item{XsqtestB, pvalB, XsqtestG, pvalG}{(optional) test statistics}
#'
#' @export
emmremlMultivariate_core <- function(Y, X, Z, K,
                                     varBhat = FALSE, varGhat = FALSE,
                                     PEVGhat = FALSE, test = FALSE,
                                     tolpar = 1e-6, tolparinv = 1e-6,
                                     maxIter = 100) {

  if (sum(is.na(Y)) > 0) stop("Missing values in Y are not supported")

  # --- dimensions (v3.1 convention: traits=rows, individuals=cols) ---
  d <- nrow(Y)   # number of traits
  n <- ncol(Y)   # number of individuals (= ncol(X) = nrow(ZKZt))

  # Precompute ZKZ' and its eigen-decomposition (core of v3.1 ECM)
  Zt <- t(Z)                        # N x n  ->  t(Z) is n x N ... wait, Z is n x N
  KZt <- tcrossprod(K, Z)           # N x n
  ZKZt <- Z %*% KZt                 # n x n

  eigZKZt <- eigen(ZKZt)
  # Project Y and X into eigenvector space (v3.1 convention)
  Yt <- Y %*% eigZKZt$vectors       # d x n
  Xt <- X %*% eigZKZt$vectors       # m x n

  # --- initialise variance components ---
  Vgt <- cov(t(Y)) / 2
  Vet <- cov(t(Y)) / 2

  # Initial B in eigenvector space (OLS)
  XttinvXtXtt <- t(Xt) %*% solve(tcrossprod(Xt))
  Bt <- Yt %*% XttinvXtXtt

  # --- ECM inner functions (identical to v3.1) ---
  ECM1 <- function(ytl, xtl, Vgt, Vet, Bt, deltal) {
    Vlt <- deltal * Vgt + Vet
    invVlt <- solve(Vlt + tolparinv * diag(d))
    gtl <- deltal * Vgt %*% invVlt %*% (ytl - Bt %*% xtl)
    Sigmalt <- deltal * Vgt - deltal * Vgt %*% invVlt %*% (deltal * Vgt)
    list(Vlt = Vlt, gtl = gtl, Sigmalt = Sigmalt)
  }

  Vgfunc <- function(l) {
    Vgl <- tcrossprod(outfromECM1[[l]]$gtl)
    (1/n) * (1/eigZKZt$values[l]) * (Vgl + outfromECM1[[l]]$Sigmalt)
  }

  Vefunc <- function(l) {
    etl <- Yt[, l] - Bt %*% Xt[, l] - outfromECM1[[l]]$gtl
    (1/n) * (tcrossprod(etl) + outfromECM1[[l]]$Sigmalt)
  }

  # --- main ECM loop ---
  Vetm1 <- Vet
  convergence <- FALSE
  iterations <- 0L

  repeat {
    iterations <- iterations + 1L
    Vetm1 <- Vet

    outfromECM1 <- lapply(seq_len(n), function(l) {
      ECM1(ytl = Yt[, l], xtl = Xt[, l],
           Vgt = Vgt, Vet = Vet, Bt = Bt,
           deltal = eigZKZt$values[l])
    })

    Gt <- sapply(outfromECM1, function(x) cbind(x$gtl))
    Bt <- (Yt - Gt) %*% XttinvXtXtt

    listVgts <- lapply(seq_len(n), Vgfunc)
    Vgt <- Reduce("+", listVgts)

    listVets <- lapply(seq_len(n), Vefunc)
    Vet <- Reduce("+", listVets)

    convnum <- tryCatch(
      abs(sum(diag(Vet - Vetm1))) / abs(sum(diag(Vetm1))),
      error = function(e) Inf
    )

    if (convnum < tolpar) { convergence <- TRUE; break }
    if (iterations >= maxIter) break
  }

  # --- BLUPs for G (Gpred), same as v3.1 ---
  HobsInv <- solve(
    kronecker(ZKZt, Vgt) + kronecker(diag(n), Vet) + tolparinv * diag(d * n)
  )
  ehat <- matrix(Y - Bt %*% X, ncol = 1, byrow = FALSE)
  HobsInve <- HobsInv %*% ehat
  varvecG <- kronecker(K, Vgt)
  gpred <- varvecG %*% (kronecker(t(Z), diag(d))) %*% HobsInve
  Gpred <- matrix(gpred, nrow = d, byrow = FALSE)
  colnames(Gpred) <- rownames(K)

  # --- optional variance / test quantities ---
  Xforvec <- kronecker(t(X), diag(d))
  Zforvec <- kronecker(Z, diag(d))
  ZKforvec <- Zforvec %*% varvecG

  P_mat <- NULL
  varGhat_out <- NULL
  varBhat_out <- NULL
  PEVGhat_out <- NULL
  XsqtestB <- NULL; p_adjBhat <- NULL
  XsqtestG <- NULL; p_adjGhat <- NULL

  if (varGhat || PEVGhat || test) {
    P_mat <- HobsInv - HobsInv %*% Xforvec %*%
      solve(crossprod(Xforvec, HobsInv %*% Xforvec), crossprod(Xforvec, HobsInv))
    varGhat_mat <- crossprod(ZKforvec, P_mat) %*% ZKforvec
    varGhat_out <- diag(varGhat_mat)
  }

  if (PEVGhat) {
    PEVGhat_out <- diag(varvecG - varGhat_mat)
  }

  if (varBhat || test) {
    if (is.null(P_mat)) {
      P_mat <- HobsInv - HobsInv %*% Xforvec %*%
        solve(crossprod(Xforvec, HobsInv %*% Xforvec), crossprod(Xforvec, HobsInv))
    }
    varBhat_mat <- solve(crossprod(Xforvec, HobsInv %*% Xforvec))
    varBhat_out <- diag(varBhat_mat)
  }

  if (test) {
    XsqtestG <- matrix(gpred, ncol = 1)^2 / varGhat_out
    XsqtestG <- matrix(matrix(XsqtestG, ncol = d, byrow = TRUE), ncol = 1, byrow = TRUE)
    pGhat <- pchisq(XsqtestG, df = 1, lower.tail = FALSE, log.p = FALSE)
    p_adjGhat <- sapply(p.adjust.methods, function(m) p.adjust(pGhat, m))

    XsqtestB <- matrix(Bt, ncol = 1)^2 / varBhat_out
    pBhat <- pchisq(XsqtestB, df = 1, lower.tail = FALSE, log.p = FALSE)
    p_adjBhat <- sapply(p.adjust.methods, function(m) p.adjust(pBhat, m))
  }

  result <- list(
    Bhat       = Bt,
    Vg         = Vgt,
    Ve         = Vet,
    Gpred      = Gpred,
    iterations = iterations,
    convergence = convergence,
    varBhat    = varBhat_out,
    varGhat    = varGhat_out,
    PEVGhat    = PEVGhat_out,
    XsqtestB   = XsqtestB,
    pvalB      = p_adjBhat,
    XsqtestG   = XsqtestG,
    pvalG      = p_adjGhat
  )
  class(result) <- c("emmreml_mv", "emmreml")
  result
}


#' Fit a Multivariate Gaussian Mixed Model — formula or matrix interface
#'
#' @description
#' Two calling conventions are supported:
#'
#' \strong{Matrix interface (v3.1 compatible):}
#' \preformatted{emmremlMultivariate(Y, X, Z, K, varBhat, varGhat, PEVGhat,
#'                    test, tolpar, tolparinv)}
#' where Y, X are (traits x individuals) matrices and Z is (individuals x levels).
#'
#' \strong{Formula interface (v5.x):}
#' \preformatted{emmremlMultivariate(formula, random, data, K, verbose,
#'                    tolParConv, tolParInv, maxIter, ...)}
#'
#' Both interfaces call \code{emmremlMultivariate_core} and return the same
#' result structure.
#'
#' @param Y_or_formula Either a (d x n) response matrix \strong{or} a formula
#'   such as \code{cbind(y1, y2) ~ x}.
#' @param X_or_random Either a (m x n) fixed-effects matrix \strong{or} an
#'   optional random-effects formula (or \code{NULL}).
#' @param Z_or_data Either a (n x N) incidence matrix \strong{or} a data frame.
#' @param K Covariance/kernel matrix (N x N).
#' @param varBhat Logical (matrix interface). Default \code{FALSE}.
#' @param varGhat Logical (matrix interface). Default \code{FALSE}.
#' @param PEVGhat Logical (matrix interface). Default \code{FALSE}.
#' @param test Logical (matrix interface). Default \code{FALSE}.
#' @param tolpar Convergence tolerance (matrix interface). Default \code{1e-6}.
#' @param tolparinv Regularisation tolerance (matrix interface). Default \code{1e-6}.
#' @param verbose Logical (formula interface). Default \code{TRUE}.
#' @param tolParConv Convergence tolerance (formula interface). Default \code{1e-6}.
#' @param tolParInv Regularisation tolerance (formula interface). Default \code{1e-6}.
#' @param maxIter Maximum ECM iterations. Default \code{100}.
#' @param ... Ignored (for forward compatibility).
#'
#' @return A list of class \code{emmreml_mv} — see \code{emmremlMultivariate_core}
#'   for full description of return values.
#'
#' @export
emmremlMultivariate <- function(Y_or_formula,
                                X_or_random = NULL,
                                Z_or_data,
                                K,
                                varBhat    = FALSE,
                                varGhat    = FALSE,
                                PEVGhat    = FALSE,
                                test       = FALSE,
                                tolpar     = 1e-6,
                                tolparinv  = 1e-6,
                                verbose    = TRUE,
                                tolParConv = 1e-6,
                                tolParInv  = 1e-6,
                                maxIter    = 100,
                                ...) {

  # ---- detect calling convention ----
  if (inherits(Y_or_formula, "formula")) {
    # ---- formula interface ----
    formula <- Y_or_formula
    random  <- X_or_random
    data    <- Z_or_data

    if (!is.data.frame(data)) stop("'data' must be a data frame")
    if (!is.matrix(K))        stop("K must be a matrix")

    if (verbose) cat("Preparing model matrices...\n")

    mf <- model.frame(formula, data, na.action = na.pass)
    X_mat <- model.matrix(formula, mf)
    Y_mat <- model.response(mf)

    if (!is.matrix(Y_mat)) stop("Response must be multivariate (use cbind)")

    n <- nrow(Y_mat)
    Z_mat <- diag(n)

    if (!is.null(random) && inherits(random, "formula")) {
      Z_mat <- diag(n)   # extend here if needed
    } else if (!is.null(random)) {
      stop("'random' must be a formula or NULL")
    }

    if (ncol(Z_mat) != nrow(K))
      stop("ncol(Z) must equal nrow(K)")

    # Convert to v3.1 convention: traits=rows, individuals=cols
    Y_core <- t(Y_mat)   # d x n
    X_core <- t(X_mat)   # m x n

    if (verbose) cat("Fitting multivariate model using ECM algorithm...\n")

    result <- emmremlMultivariate_core(
      Y = Y_core, X = X_core, Z = Z_mat, K = K,
      tolpar = tolParConv, tolparinv = tolParInv,
      maxIter = maxIter
    )
    result$call <- match.call()
    if (verbose) cat("Model fitting complete\n")
    return(result)

  } else {
    # ---- matrix interface (v3.1 compatible) ----
    Y <- Y_or_formula
    X <- X_or_random
    Z <- Z_or_data

    if (is.null(X)) stop("X must be provided")
    if (!is.matrix(Y)) stop("Y must be a matrix (traits x individuals)")
    if (!is.matrix(X)) stop("X must be a matrix (fixed effects x individuals)")
    if (!is.matrix(Z)) stop("Z must be a matrix (individuals x levels)")
    if (!is.matrix(K)) stop("K must be a matrix")

    return(emmremlMultivariate_core(
      Y = Y, X = X, Z = Z, K = K,
      varBhat = varBhat, varGhat = varGhat,
      PEVGhat = PEVGhat, test = test,
      tolpar = tolpar, tolparinv = tolparinv,
      maxIter = maxIter
    ))
  }
}
