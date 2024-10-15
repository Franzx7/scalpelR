\name{plotDIU}
\alias{plotDIU}
\title{Plot Differential Isoform Usage (DIU)}
\description{
  The \code{plotDIU} function visualizes isoform usage across groups from a Seurat object based on differential isoform usage (DIU) data. It calculates isoform fractions and generates a violin plot of the isoform expression across clusters.
}
\usage{
plotDIU(seurat_obj, DIU_table, group.by = NULL, gene, extended=F, type="fraction")
}
\arguments{
  \item{seurat_obj}{A Seurat object containing single-cell RNA-seq data.}
  \item{DIU_table}{A data frame containing differential isoform usage information.
  It must include a column \code{gene_name} and \code{gene_tr}.}
  \item{group.by}{The metadata column in the Seurat object to group cells by.}
  \item{gene}{The gene for which isoform usage will be plotted.}
  \item{extended}{To add the gene expression plot for each cluster (T) or not(F)}
  \item{type}{To show the fraction of Isoform or normalized expression in clusters}
}
\details{
  This function computes the isoform fractions for the specified gene and
  visualizes them using a violin plot. It relies on the differential isoform usage
  data to summarize the expression of each isoform across different groups or
  clusters defined in the Seurat object.
}
\value{
  A \code{ggplot} object showing the isoform fractions in the specified groups.
}
\examples{
\dontrun{
# Example usage of plotDIU
plotDIU(seurat_obj = my_seurat_obj,
        DIU_table = my_DIU_table,
        group.by = "Cluster",
        gene = "GeneA")
}
}
\seealso{
  \code{\link[ggplot2]{ggplot}}, \code{\link[Seurat]{Seurat}}
}
\author{
  Franz AKE, IDIBELL-PCMRC, Bioinformatics Researcher (fake@idibell.cat)
}
\keyword{visualization}
\keyword{isoform usage}
