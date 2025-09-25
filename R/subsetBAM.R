#' Subset BAM File by Barcodes
#'
#' Subsets a BAM file to include only reads with specified barcode tags, writing the result to a new BAM file.
#'
#' @param bam_file Path to the input BAM file to be subset.
#' @param output_file Path to the output BAM file to be written.
#' @param barcodes Character vector of cell barcodes to retain.
#' @param tag Character. The tag in the BAM file identifying the barcode (e.g., \code{"BC"}).
#' @param samtools_bin Character. Path to the \code{samtools} binary. Default is \code{"samtools"}.
#' @param cores Integer. Number of cores to use for filtering. Default is \code{4}.
#'
#' @details
#' The function writes the selected barcodes to a temporary file and uses \code{samtools view} to extract reads with matching barcode tags into a new BAM file, followed by indexing.
#'
#' @return
#' NULL. The function is called for its side effects (writing a subset BAM and its index file).
#'
#' @examples
#' \dontrun{
#' subsetBAM(
#'   bam_file = "input.bam",
#'   output_file = "subset.bam",
#'   barcodes = c("AAACCTGAGC", "AAACCTGCTT"),
#'   tag = "BC",
#'   samtools_bin = "samtools",
#'   cores = 4
#' )
#' }
#'
#' @export
subsetBAM <- function(
    bam_file,
    output_file,
    barcodes,
    tag,
    samtools_bin = "samtools",
    cores = 4) {
  # Function to subset BAM file based on barcode tags...
  data.frame(bcs = barcodes) %>% data.table::fwrite(file = "bc.tmp", col.names = F)
  samtools.exp <- paste0(samtools_bin, " view -b -@ ", cores, " -D ", tag, ":bc.tmp ", bam_file, " > ", output_file)
  samtools.ind <- paste0(samtools_bin, " index -@ ", cores, " ", output_file)
  message("Filtering BAM file...")
  system(samtools.exp)
  system(samtools.ind)
  system("rm -r ./bc.tmp")
}
