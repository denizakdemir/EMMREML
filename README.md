# EMMREML: Fitting Mixed Models with Known Covariance Structures

## Overview

The EMMREML package provides efficient algorithms for fitting mixed models with known covariance structures. It is particularly useful for genomic prediction, quantitative genetics, and other applications where the genetic relationships between individuals or items are known.

## Key Features

- **Formula interface**: Easy to use R-style formula interface for model specification
- **Multiple model types**:
  - Single-kernel mixed models (EMMA algorithm)
  - Multi-kernel mixed models (weighted combination of kernels)
  - Multivariate mixed models (ECM algorithm)
- **Comprehensive output**: Variance components, fixed and random effects, standard errors, p-values
- **Diagnostics**: AIC, BIC, log-likelihood for model comparison
- **S3 methods**: Consistent print, summary, and predict methods

## Installation

```r
# Install from CRAN
install.packages("EMMREML")

# Or the development version from GitHub:
# install.packages("devtools")
# devtools::install_github("dakdemir/EMMREML")
```

## Basic Usage

### Single Kernel Model

```r
library(EMMREML)

# Example data
n <- 200
M1 <- matrix(rnorm(n*300), nrow=n)
K1 <- cov(t(M1))
K1 <- K1/mean(diag(K1))  # Normalize

# Generate phenotype data
covY <- 2*K1 + 1*diag(n)  # Genetic + residual covariance
Y <- 10 + crossprod(chol(covY), rnorm(n))  # Phenotype with mean 10

# Create data frame
df <- data.frame(y = Y, x = rnorm(n), id = 1:n)

# Fit model using formula interface
mod <- emmreml(
  formula = y ~ x,
  data = df,
  K = K1,
  getVarBeta = TRUE,
  verbose = TRUE
)

# View results
print(mod)
summary(mod)

# Predictions
pred <- predict(mod)
```

### Multi-Kernel Model

```r
# Fit model with multiple kernels/random effects
mod_multi <- emmremlMultiKernel(
  formula = y ~ x,
  random = list(id = ~ id, group = ~ group - 1),
  data = df_multi,
  K = list(K1, K2),
  getPEV = TRUE
)
```

### Multivariate Model

```r
# Fit multivariate model for multiple traits
mod_mv <- emmremlMultivariate(
  formula = cbind(y1, y2, y3) ~ x,
  random = ~ id,
  data = df_mv,
  K = K1,
  getVarB = TRUE
)
```

## References

For more details on the underlying algorithms:

- Kang HM, et al. (2008). Efficient control of population structure in model organism association mapping. *Genetics*, 178(3), 1709-1723.
- Meng X-L, Rubin DB (1993). Maximum likelihood estimation via the ECM algorithm: A general framework. *Biometrika*, 80(2), 267-278.

## Version 5.1 Update

The 5.1 release includes:

- Fixed multivariate implementation with corrected matrix dimensions and eigenvector projection
- Streamlined codebase with removal of duplicate/debug files
- Enhanced GWAS functionality with accelerated computation options
- Updated dependencies and compatibility
- Additional bug fixes and stability improvements

See the TODO.md file for upcoming improvements.

## Citation

If you use EMMREML in your research, please cite:

Akdemir D, Godfrey OU (2025). EMMREML: Fitting Mixed Models with Known Covariance Structures. R package version 5.1, https://github.com/dakdemir/EMMREML.
