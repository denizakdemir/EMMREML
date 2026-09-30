# EMMREML Package ToDo List

This file contains a list of issues that should be addressed in future versions of the package.

## Resolved in v5.2

1. `emmremlMultivariate` backward compatibility restored:
   - Matrix interface `(Y, X, Z, K, varBhat, varGhat, PEVGhat, test, tolpar, tolparinv)` is the primary API again
   - `Gpred` restored in all return values
   - All test outputs (`XsqtestB`, `pvalB`, `XsqtestG`, `pvalG`, `varBhat`, `varGhat`, `PEVGhat`) restored
   - Formula interface kept as a thin wrapper
2. ECM convergence fixed: reverted to v3.1 ECM update equations (correct per-eigenvalue updates)
3. Identity-matrix shortcut now uses `is_identity_matrix()` from `utils.R` everywhere (emmreml, emmremlMultiKernel_core, gwas)
4. `impute` (Bioconductor) dependency: graceful fallback to mean imputation with informative message
5. `Matrix::.bdiag` call in `emmremlMultiKernel` already uses `Matrix::` prefix (correct)

## Documentation Issues

1. Add proper roxygen documentation for all exported functions:
   - `calculate_allele_frequencies` 
   - `calculate_kinship`
   - `compare_kinship`
   - `convert_marker_coding` 
   - `demo_K_optimization`
   - `demo_gwas`
   - `demo_kinship`
   - `demo_multivariate`
   - `filter_markers`
   - `generate_example_data`
   - `gwas`
   - `identify_marker_coding` 
   - `impute_markers`
   - `marker_statistics`
   - `plot_marker_stats`
   - `subset_markers_by_maf`
   - `summarize_marker_data`
   - `visualize_kinship`

2. Fix Rd cross-references in kinship.Rd documentation
3. Update `emmremlMultivariate.Rd` to document both calling conventions

## Code Issues

1. Address R CMD check warnings and notes
2. Update vignette to include examples of all main functions

## Dependencies

1. Add proper imports for visualizations in ggplot2 and viridis packages

## Testing

1. Create comprehensive test suite for all main functions
2. Add regression tests that verify v3.1 matrix-interface outputs match v5.2 outputs numerically
3. Add more test cases for the GWAS functionality