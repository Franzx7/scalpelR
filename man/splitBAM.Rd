\name{splitBAM}
\alias{splitBAM}
\title{Split BAM Files by Condition}
\description{
  This function splits a BAM file according to the cell barcodes associated with
  a condition level defined in the Seurat object. It generates one BAM file for
  each condition level.
}
\usage{
splitBAM(seurat_obj, barcode_tag = "CB", bam_file, samtools_bin = 'samtools',
  group.by, ncores = 4, output_dir)
}
\arguments{
  \item{seurat_obj}{A \code{Seurat} object containing cell metadata.}
  \item{barcode_tag}{A string indicating the barcode tag (default is "BC").}
  \item{bam_file}{Path to the BAM file to be split.}
  \item{samtools_bin}{Path to the samtools binary (default is 'samtools').}
  \item{group.by}{Column name in the Seurat object metadata to group cells.}
  \item{ncores}{Number of cores to use for filtering (default is 4).}
  \item{output_dir}{Directory where the output BAM files will be saved.}
}
\value{
  Returns \code{NULL}. BAM files are generated in the specified output directory.
}
\details{
  The function retrieves barcodes from the Seurat object and formats them
  according to the specified condition levels. Ensure
  that the cell barcodes in the Seurat object are formatted similarly to those
  in the BAM file !
}
\examples{
# Load the necessary libraries
library(Seurat)
library(data.table)

# Example usage (assuming you have a Seurat object and a BAM file)
# seurat_obj <- readRDS("path_to_your_seurat_object.rds")
# splitBAM(seurat_obj = seurat_obj,
#          bam_file = "path_to_your_file.bam",
#          group.by = "condition_column",
#          output_dir = "output_directory")
}

\author{
  Franz, University of Barcelona, Bioinformatics Researcher
}
