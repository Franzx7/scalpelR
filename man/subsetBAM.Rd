\name{subsetBAM}
\title{Subset BAM File by Barcodes}
\description{
Subsets a BAM file to include only reads with specified barcode tags, writing the result to a new BAM file.
}
\usage{
subsetBAM(
  bam_file,
  output_file,
  barcodes,
  tag,
  samtools_bin = "samtools",
  cores = 4
)
}
\arguments{
  \item{bam_file}{Path to the input BAM file to be subset.}
  \item{output_file}{Path to the output BAM file to be written.}
  \item{barcodes}{Character vector of cell barcodes to retain.}
  \item{tag}{Character. The tag in the BAM file identifying the barcode (e.g., \code{"BC"}).}
  \item{samtools_bin}{Character. Path to the \code{samtools} binary. Default is \code{"samtools"}.}
  \item{cores}{Integer. Number of cores to use for filtering. Default is \code{4}.}
}
\details{
The function writes the selected barcodes to a temporary file and uses \code{samtools view} to extract reads with matching barcode tags into a new BAM file, followed by indexing.
}
\value{
NULL. The function is called for its side effects (writing a subset BAM and its index file).
}
\examples{
\dontrun{
subsetBAM(
  bam_file = "input.bam",
  output_file = "subset.bam",
  barcodes = c("AAACCTGAGC", "AAACCTGCTT"),
  tag = "BC",
  samtools_bin = "samtools",
  cores = 4
)
}
}
\author{Franz AKE}