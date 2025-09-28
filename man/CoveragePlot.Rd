\name{CoveragePlot}
\title{Plot Read Coverage and Gene/Transcript Structure from BAM Files}
\description{
Generates a visualization of gene and transcript annotation alongside read coverage from one or more BAM files using the Gviz package.
}
\usage{
CoveragePlot(
  genome_gr,
  bamfiles,
  gene_in,
  features = NULL,
  sizes = NULL,
  genome_sp = "mm10",
  distZOOM = NULL,
  annot_tab = NULL,
  filter_trs = FALSE,
  extend.left = 1000,
  extend.right = 1000,
  samtools.bin = "samtools"
)
}
\arguments{
  \item{genome_gr}{A GRanges object containing genome annotation (e.g., from GTF).}
  \item{bamfiles}{A named vector of BAM file paths. Names will be used as track labels.}
  \item{gene_in}{The gene name of interest (must be present in \code{genome_gr}).}
  \item{features}{Optional. A vector of transcript names to be displayed for the gene (\code{NULL} for all transcripts).}
  \item{sizes}{Optional. A numeric vector specifying the number of cells in each BAM file (for normalization purposes).}
  \item{genome_sp}{The genome build identifier, either \code{"mm10"} or \code{"hg38"}. Default is \code{"mm10"}.}
  \item{distZOOM}{Optional. Zoom to a subregion (numeric, in bp) from gene boundary. If \code{NULL}, plot the entire gene region.}
  \item{annot_tab}{Optional. Annotation table. Not used internally but reserved for future use.}
  \item{filter_trs}{Logical. Whether to filter transcripts. Default is \code{FALSE}.}
  \item{extend.left}{Numeric. Number of bases to extend to the left of the gene region. Default is 1000.}
  \item{extend.right}{Numeric. Number of bases to extend to the right of the gene region. Default is 1000.}
  \item{samtools.bin}{Character. Path to the \code{samtools} executable. Default is \code{"samtools"}.}
}
\details{
For each BAM file, coverage is computed over the gene region (optionally including selected transcripts and extensions)
using \code{samtools view} and \code{samtools depth}. The gene region and transcript structure are shown using Gviz tracks.
The function can zoom to either end of the gene, depending on strand and the \code{distZOOM} parameter.
}
\value{
Invisibly returns the Gviz plot object (plot is drawn as a side effect).
}
\examples{
\dontrun{
CoveragePlot(
  genome_gr = gtf_gr,
  bamfiles = c(Sample1 = "sample1.bam", Sample2 = "sample2.bam"),
  gene_in = "GeneA",
  features = c("Transcript1", "Transcript2"),
  sizes = c(500, 600),
  genome_sp = "mm10",
  extend.left = 2000,
  extend.right = 2000
)
}
}
\author{Franz AKE}