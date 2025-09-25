\name{plotExpression}
\title{Plot Relative Isoform or Feature Expression Across Groups}
\description{
Generates boxplots of relative expression for selected features (e.g., isoforms or genes) across groups in a Seurat object.
}
\usage{
plotExpression(
  seurat_obj,
  features,
  group_var,
  plot = "group",
  assay = "RNA",
  log_normalization = TRUE
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing single-cell RNA-seq data.}
  \item{features}{Character vector of features (e.g., genes or isoforms) to plot.}
  \item{group_var}{The metadata column in the Seurat object to group cells by.}
  \item{plot}{Character. Whether to plot by group (\code{"group"}) or by features (\code{"features"}). Default is \code{"group"}.}
  \item{assay}{Character. The Seurat assay to use. Default is \code{"RNA"}.}
  \item{log_normalization}{Logical. Whether to log-normalize the expression values (\code{TRUE}) or plot raw counts (\code{FALSE}). Default is \code{TRUE}.}
}
\details{
The function extracts counts for the specified features from the specified assay in the Seurat object. It plots boxplots of expression values, grouped either by metadata group or by features, and can display either raw or log-normalized counts.
}
\value{
A ggplot2 object showing boxplots of expression values for the selected features across the specified groups.
}
\examples{
\dontrun{
plotExpression(
  seurat_obj = my_seurat_obj,
  features = c("Isoform1", "Isoform2"),
  group_var = "Cluster",
  plot = "group",
  assay = "RNA",
  log_normalization = TRUE
)
}
}
\author{Franz AKE}