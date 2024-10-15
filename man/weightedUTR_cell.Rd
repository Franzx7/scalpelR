\name{weightedUTR_cell}
\alias{weightedUTR_cell}
\title{Calculate Weighted UTR Length per Cell}
\description{
  The \code{weightedUTR_cell} function computes the weighted UTR length per cell
  using transcript information from a Seurat object. It requires a gene annotation
  file with UTR and CDS information to calculate UTR lengths across isoforms.
}
\usage{
weightedUTR_cell(seurat_obj, features = NULL, annotation_gr, group.by = NULL,
assay = "RNA", min.counts = 3)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing single-cell RNA-seq data.}
  \item{features}{A vector of features (genes/transcripts) to subset the Seurat
  object by. Default is \code{NULL}.}
  \item{annotation_gr}{A GRanges object containing gene annotation, including
  UTR and CDS regions.}
  \item{group.by}{The metadata column in the Seurat object to group cells by.}
  \item{assay}{The assay to use for expression values. Default is "RNA".}
  \item{min.counts}{The minimum number of counts required to consider a
  transcript for further analysis. Default is 3.}
}
\details{
  This function calculates the weighted UTR length for each cell by combining
  gene annotations (including UTR and CDS regions) and isoform expression levels
  from the Seurat object. The UTR length is weighted by the relative expression
  of each isoform in a given cell.
}
\value{
  A \code{data.table} object containing the weighted UTR length for each cell,
  along with cluster annotations.
}
\examples{
\dontrun{
# Example usage of weightedUTR_cell
weightedUTR_cell(seurat_obj = my_seurat_obj,
                 features = NULL,
                 annotation_gr = my_annotation,
                 group.by = "Cluster",
                 assay = "RNA",
                 min.counts = 3)
}
}
\seealso{
  \code{\link[Seurat]{Seurat}}, \code{\link[GenomicRanges]{GRanges}}
}
\author{
  Franz, University of Barcelona, Bioinformatics Researcher
}
\keyword{RNA-seq}
\keyword{Seurat}
