\name{FindIsoforms}
\title{Identify Differentially Expressed Isoforms Across Conditions}
\description{
Identifies isoform expression with differential expression across specified conditions in single-cell RNA-seq data using a Seurat object.
}
\usage{
FindIsoforms(
  seurat_obj,
  group_by,
  split_by = NULL,
  assay = "RNA",
  threshold_pval = 0.05,
  threshold_abund = 0.1,
  threshold_sd = 0.05
)
}
\arguments{
  \item{seurat_obj}{The input Seurat object containing scRNA-seq data.}
  \item{group_by}{The metadata variable used to group cells for isoform expression comparison.}
  \item{split_by}{Optional. Metadata variable used to split the Seurat object before analysis. If \code{NULL}, no split is performed. Default is \code{NULL}.}
  \item{assay}{Character. The Seurat assay to use for analysis. Default is \code{"RNA"}.}
  \item{threshold_pval}{Numeric. Adjusted p-value threshold for significance. Default is \code{0.05}.}
  \item{threshold_abund}{Numeric. Minimum isoform abundance proportion per group. Default is \code{0.1}.}
  \item{threshold_sd}{Numeric. Minimum standard deviation of counts across groups for an isoform to be tested. Default is \code{0.05}.}
}
\details{
The function aggregates isoform (transcript) expression per group, filters isoforms based on abundance and variability, and applies a chi-squared test to assess differential isoform usage. Adjusted p-values (Benjamini-Hochberg) are reported. If \code{split_by} is specified and different from \code{group_by}, the analysis is performed separately for each split.
}
\value{
A data.frame (or data.table) containing isoforms passing the significance threshold, with statistics for differential expression. If splitting, results from each split are combined.
}
\examples{
\dontrun{
FindIsoforms(
  seurat_obj = my_seurat,
  group_by = "condition",
  split_by = NULL,
  assay = "RNA",
  threshold_pval = 0.05,
  threshold_abund = 0.1,
  threshold_sd = 0.05
)
}
}
\author{}