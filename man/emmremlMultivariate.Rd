\name{emmremlMultivariate}
\alias{emmremlMultivariate}
\title{
Fit a Multivariate Gaussian Mixed Model with Known Covariance Structure
}
\description{
Fits a multivariate Gaussian mixed model of the form \eqn{Y = XB + ZG + E},
where \eqn{vec(G) \sim N(0, K \otimes V_g)} and
\eqn{vec(E) \sim N(0, I_n \otimes V_e)}, using the ECM algorithm.

Two calling conventions are supported:

\strong{Matrix interface (v3.1 compatible):}
\preformatted{emmremlMultivariate(Y, X, Z, K, varBhat, varGhat, PEVGhat,
                   test, tolpar, tolparinv)}
where \code{Y} and \code{X} are (traits x individuals) matrices and
\code{Z} is an (individuals x levels) incidence matrix.

\strong{Formula interface (v5.x):}
\preformatted{emmremlMultivariate(formula, random, data, K, verbose,
                   tolParConv, tolParInv, maxIter)}
}
\usage{
emmremlMultivariate(Y_or_formula, X_or_random = NULL, Z_or_data, K,
                    varBhat = FALSE, varGhat = FALSE, PEVGhat = FALSE,
                    test = FALSE, tolpar = 1e-6, tolparinv = 1e-6,
                    verbose = TRUE, tolParConv = 1e-6, tolParInv = 1e-6,
                    maxIter = 100, ...)
}
\arguments{
  \item{Y_or_formula}{
    Either a (d x n) response matrix (traits in rows, individuals in columns)
    for the matrix interface, or a formula such as \code{cbind(y1, y2) ~ x}
    for the formula interface.
  }
  \item{X_or_random}{
    Either a (m x n) fixed-effects design matrix for the matrix interface,
    or an optional random-effects formula (or \code{NULL}) for the formula interface.
  }
  \item{Z_or_data}{
    Either a (n x N) incidence matrix for the matrix interface,
    or a data frame for the formula interface.
  }
  \item{K}{(N x N) covariance/kernel matrix for random effects.}
  \item{varBhat}{Logical, whether to return variance of B estimates. Default \code{FALSE}.}
  \item{varGhat}{Logical, whether to return variance of G predictions. Default \code{FALSE}.}
  \item{PEVGhat}{Logical, whether to return prediction error variance of G. Default \code{FALSE}.}
  \item{test}{Logical, whether to perform chi-square significance tests. Default \code{FALSE}.}
  \item{tolpar}{Convergence tolerance for variance components (matrix interface). Default \code{1e-6}.}
  \item{tolparinv}{Regularisation added to matrices before inversion (matrix interface). Default \code{1e-6}.}
  \item{verbose}{Logical, whether to print progress (formula interface). Default \code{TRUE}.}
  \item{tolParConv}{Convergence tolerance (formula interface). Default \code{1e-6}.}
  \item{tolParInv}{Regularisation tolerance (formula interface). Default \code{1e-6}.}
  \item{maxIter}{Maximum number of ECM iterations. Default \code{100}.}
  \item{...}{Ignored (for forward compatibility).}
}
\value{
A list of class \code{emmreml_mv} containing:
\item{Bhat}{Fixed effects estimates in eigenvector space (m x n matrix).}
\item{Vg}{Estimated genetic covariance matrix (d x d).}
\item{Ve}{Estimated residual covariance matrix (d x d).}
\item{Gpred}{BLUPs for G (d x N matrix).}
\item{iterations}{Number of ECM iterations run.}
\item{convergence}{Logical: did the algorithm converge?}
\item{varBhat}{(if requested) Diagonal of variance matrix for vec(B).}
\item{varGhat}{(if requested) Diagonal of variance matrix for vec(G).}
\item{PEVGhat}{(if requested) Diagonal of prediction error variance for vec(G).}
\item{XsqtestB}{(if test=TRUE) Chi-square statistics for fixed effects.}
\item{pvalB}{(if test=TRUE) Adjusted p-values for fixed effects.}
\item{XsqtestG}{(if test=TRUE) Chi-square statistics for G predictions.}
\item{pvalG}{(if test=TRUE) Adjusted p-values for G predictions.}
}
\details{
The ECM algorithm projects Y and X into the eigenvector space of ZKZ',
then iterates per-eigenvalue E-steps and M-steps to update the genetic
(Vg) and residual (Ve) covariance matrices. This is the same algorithm
as in v3.1, with the addition of a \code{maxIter} guard and a robust
convergence check. See \code{\link{emmremlMultivariate_core}} for the
underlying implementation.
}
\seealso{
\code{\link{emmremlMultivariate_core}}, \code{\link{emmreml}}, \code{\link{emmremlMultiKernel}}
}
\examples{
set.seed(123)
n <- 15   # individuals
d <- 2    # traits

# --- matrix interface (v3.1 style) ---
Y <- matrix(rnorm(d * n), d, n)   # traits x individuals
X <- rbind(rep(1, n), rnorm(n))   # fixed effects x individuals
Z <- diag(n)
K_raw <- matrix(rnorm(n * n), n, n)
K <- tcrossprod(K_raw)
diag(K) <- diag(K) + 1e-4
K <- K / mean(diag(K))

result <- emmremlMultivariate(Y, X, Z, K)
print(result$Vg)
print(result$Ve)
dim(result$Gpred)   # d x n

\dontrun{
# --- formula interface (v5.x style) ---
df <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x1 = rnorm(n))
result2 <- emmremlMultivariate(
  formula = cbind(y1, y2) ~ x1,
  data = df,
  K = K,
  verbose = FALSE
)
print(result2$Vg)
}
}
