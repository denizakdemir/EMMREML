\name{kinship}
\alias{kinship}
\title{Kinship Matrix Calculation and Marker Data Functions}
\description{
  A collection of functions for kinship matrix calculation, marker data preprocessing, filtering, imputation, and visualization. 
  These functions provide comprehensive tools for handling genomic marker data and calculating various types of relationship matrices.
}
\details{
  The functions in this module fall into several categories:
  
  \strong{Marker Data Inspection:}
  \itemize{
    \item \code{\link{identify_marker_coding}} - Identify the coding system used in marker data (0/1, 0/1/2, -1/0/1, etc.)
    \item \code{\link{marker_statistics}} - Calculate various statistics for marker data (MAF, heterozygosity, etc.)
    \item \code{\link{summarize_marker_data}} - Create summary table of marker data properties
  }
  
  \strong{Marker Data Preprocessing:}
  \itemize{
    \item \code{\link{convert_marker_coding}} - Convert marker data between different coding systems
    \item \code{\link{filter_markers}} - Filter markers based on MAF, missing data, etc.
    \item \code{\link{impute_markers}} - Impute missing values in marker data
    \item \code{\link{subset_markers_by_maf}} - Select markers based on MAF distribution
  }
  
  \strong{Kinship Matrix Calculation:}
  \itemize{
    \item \code{\link{calculate_allele_frequencies}} - Calculate allele frequencies from marker data
    \item \code{\link{calculate_kinship}} - Calculate kinship matrices using various methods
    \item \code{\link{compare_kinship}} - Compare different kinship matrices
  }
  
  \strong{Visualization:}
  \itemize{
    \item \code{\link{visualize_kinship}} - Create heatmap visualizations of kinship matrices
    \item \code{\link{plot_marker_stats}} - Plot marker statistics (MAF, missing data, etc.)
  }
}
\references{
  VanRaden, P. M. (2008). Efficient methods to compute genomic predictions.
  Journal of Dairy Science, 91(11), 4414-4423.
  
  Astle, W. & Balding, D. J. (2009). Population structure and cryptic relatedness in genetic association studies.
  Statistical Science, 24(4), 451-471.
  
  Endelman, J. B. & Jannink, J. L. (2012). Shrinkage estimation of the realized relationship matrix.
  G3: Genes, Genomes, Genetics, 2(11), 1405-1413.
}
\seealso{
  \code{\link{emmreml}} for mixed model analysis using kinship matrices
  \code{\link{emmremlMultiKernel}} for mixed models with multiple kinship matrices
  \code{\link{gwas}} for genome-wide association studies using kinship matrices
}
\examples{
\dontrun{
# Generate example marker data
set.seed(123)
n_samples <- 100
n_markers <- 1000
M <- matrix(sample(c(0,1,2), n_samples*n_markers, replace=TRUE), 
            nrow=n_samples, ncol=n_markers)
colnames(M) <- paste0("SNP", 1:n_markers)
rownames(M) <- paste0("Sample", 1:n_samples)

# Identify marker coding
coding_info <- identify_marker_coding(M)
print(coding_info)

# Calculate marker statistics
stats <- marker_statistics(M, coding="012")
head(stats)

# Filter markers by MAF
M_filtered <- filter_markers(M, min_maf=0.1, coding="012")

# Calculate VanRaden kinship matrix
K_vanraden <- calculate_kinship(M_filtered, method="vanraden1", coding="012")

# Calculate correlation-based kinship matrix
K_corr <- calculate_kinship(M_filtered, method="correlation", coding="012")

# Compare different kinship matrices
comparison <- compare_kinship(
  list(K_vanraden, K_corr), 
  names=c("VanRaden", "Correlation")
)
print(comparison)

# Use in mixed model
y <- rnorm(n_samples)  # Example phenotypes
X <- cbind(1, rnorm(n_samples))  # Example fixed effects
colnames(X) <- c("intercept", "covariate")

# Fit model with kinship
model <- emmreml(y ~ X, K=K_vanraden)
summary(model)
}
}
\author{
  Deniz Akdemir
}
