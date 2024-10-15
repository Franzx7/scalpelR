\name{FindIsoforms}
\alias{FindIsoforms}
\title{Identify Isoform Expression with Differential Expression Across Conditions}
\description{
  The \code{FindIsoforms} function identifies isoform expression with differential expression across specified conditions in single-cell RNA-seq data.
}
\usage{
FindIsoforms(seurat_obj, group.by = "orig.ident", split.by = NULL, assay = "RNA",
  threshold_pval = 0.05, threshold_abund = 0.10, min.cells.fract = 0.10)
}
\arguments{
  \item{seurat_obj}{The input Seurat object containing single-cell RNA-seq data.}
  \item{group.by}{The condition to compare isoforms (default: "orig.ident").}
  \item{split.by}{Optional parameter to split the analysis by an additional condition.}
  \item{assay}{The Seurat assay to use (default: "RNA").}
  \item{threshold_pval}{Adjusted P-value threshold for significance (default: 0.05).}
  \item{threshold_abund}{Minimum isoform abundance percentage (default: 0.10).}
  \item{min.cells.fract}{Minimum fraction of cells expressing isoforms (default: 0.10).}
}
\details{
  This function performs differential expression analysis of isoforms across specified conditions within the input Seurat object. It includes checks for the validity of the input parameters, calculates isoform expression fractions, and performs a Chi-squared test to identify significant isoforms.
}
\value{
  A data frame with identified isoforms, including P-values and adjusted P-values for significant isoforms.
}
\examples{
\dontrun{
# Example usage of FindIsoforms
result <- FindIsoforms(seurat_obj = my_seurat_obj,
                       group.by = "Cluster",
                       split.by = "Treatment",
                       assay = "RNA",
                       threshold_pval = 0.05,
                       threshold_abund = 0.10,
                       min.cells.fract = 0.10)
}
}
\seealso{
  \code{\link[Seurat]{Seurat}}, \code{\link[stats]{chisq.test}}
}
\author{
  Franz, University of Barcelona, Bioinformatics Researcher
}
\keyword{isoform}
\keyword{differential expression}
\keyword{single-cell RNA-seq}
