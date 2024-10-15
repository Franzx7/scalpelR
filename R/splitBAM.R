#' Split BAM Files by Condition
#'
#' This function splits a BAM file according to the cell barcodes associated with
#' a condition level defined in the Seurat object. It generates one BAM file for
#' each condition level.
#'
#' @param seurat_obj A Seurat object containing cell metadata.
#' @param barcode_tag A string indicating the barcode tag (default is "BC").
#' @param bam_file Path to the BAM file to be split.
#' @param samtools_bin Path to the samtools binary (default is 'samtools').
#' @param group.by Column name in the Seurat object metadata to group cells.
#' @param ncores Number of cores to use for filtering (default is 4).
#' @param output_dir Directory where the output BAM files will be saved.
#'
#' @return NULL
#' @export
splitBAM <- function(seurat_obj,
                     barcode_tag = "CB",
                     bam_file,
                     samtools_bin='samtools',
                     group.by,
                     ncores = 4,
                     output_dir) {

    # Check if group.by is a valid column in the Seurat object
    if (!(group.by %in% colnames(seurat_obj@meta.data))) {
        stop(paste0("Error: '", group.by, "' not found in the Seurat object!"))
    }

    # Check if output_dir exists, create if it doesn't
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    curr.dir <- getwd()

    # Retrieve barcodes from the Seurat object and format them
    message('Retrieve barcodes from the Seurat object and format them...')
    metadata.inf <- data.table(Barcodes = rownames(seurat_obj@meta.data),
                               condition_defined = unlist(seurat_obj[[group.by]]))
    barcodes <- split(metadata.inf, metadata.inf$condition_defined)

    # Print a warning about barcode formatting
    warning("Ensure that cell barcodes in the Seurat object are formatted similarly to those in the BAM file.")

    # Iterate through each condition level
    VOID <- lapply(names(barcodes), function(x) {
        message(paste0("Processing ", x, " ..."))
        # Write barcodes tags into a txt file
        data.table::fwrite(dplyr::select(barcodes[[x]], Barcodes),
                           file = paste0(x, '.txt'), col.names = F)

        # Filter bamfile & indexing
        filt.exp <- paste0(samtools_bin, ' view -b -D ', barcode_tag, ':', curr.dir, '/', x,
                           '.txt -@ ', ncores, ' ', bam_file,' > ',
                           output_dir, '/', x, '.bam')
        indexing.exp <- paste0(samtools_bin, ' index ',
                               output_dir, '/', x, '.bam')
        system(filt.exp)
        system(indexing.exp)
        system(paste0('rm ', x, '.txt'))
    })
}
