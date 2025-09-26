#' Calculate Weighted UTR Length per Cell
#'
#' This function computes the weighted UTR length per cell using transcript information from a Seurat object.
#' It requires a gene annotation file with UTR and CDS information to calculate UTR lengths across isoforms.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param features A vector of features (genes/transcripts) to subset the Seurat object by. Default is NULL.
#' @param annotation_gr A GRanges object containing gene annotation, including UTR and CDS regions.
#' @param group.by The metadata column in the Seurat object to group cells by.
#' @param assay The assay to use for expression values, default is "RNA".
#' @param min.counts The minimum number of counts required to consider a transcript for further analysis. Default is 3.
#'
#' @return A data.table object containing the weighted UTR length for each cell, along with cluster annotations.
#'
#' @examples
#' \dontrun{
#' weightedUTR_cell(seurat_obj = my_seurat_obj, features = NULL,
#'                  annotation_gr = my_annotation, group.by = "Cluster",
#'                  assay = "RNA", min.counts = 3)
#' }
#'
#' @export
weightedUTR_cell <- function(
        seurat_obj,
        features = NULL,
        annotation_gr,
        group.by = NULL,
        assay = "RNA",
        min.counts = 3) {
    # Checks
    if (!is.null(features)) {seurat_obj <- subset(seurat_obj, features = features)}
    if (is.null(group.by)) stop("Error! Provide group.by argument.")
    if (!(group.by %in% colnames(seurat_obj[[]]))) stop(paste0("Error! ", group.by, " not found in Seurat object."))

    message("Processing annotation file...(1/3)")
    # a. Subset isoforms from annotation
    my_annot <- annotation_gr[annotation_gr$type=="CDS" | annotation_gr$type=="UTR"][,c("type","gene_name","transcript_name")] %>%
        data.frame() %>%
        dplyr::distinct(seqnames, start, end, strand, type, width, gene_name, transcript_name) %>%
        dplyr::group_by(transcript_name) %>%
        dplyr::mutate(end_CDS = dplyr::if_else(strand == "-", min(start[which(type == "CDS")]), max(end[which(type == "CDS")]))) %>%
        dplyr::filter((strand == "-" & end < end_CDS) | (strand == "+" & start > end_CDS)) %>%
        dplyr::mutate(UTRlength = sum(width)) %>%
        dplyr::distinct(gene_name, transcript_name, UTRlength) %>%
        data.table::data.table()

    message("Processing isoform expression...(2/3)")
    # Extract isoform expression from Seurat object
    my_expr <- tibble::as_tibble(data.frame(seurat_obj[[assay]]$counts), rownames = "gene_tr") %>%
        tidyr::separate(col = gene_tr, into = c("gene_name", "transcript_name"), sep = "\\*\\*\\*", remove = FALSE) %>%
        reshape2::melt(id.vars = c('gene_tr', 'gene_name', 'transcript_name')) %>% suppressWarnings() %>%
        dplyr::filter(value > min.counts & transcript_name %in% my_annot$transcript_name) %>%
        dplyr::mutate(variable = stringr::str_replace(variable, "\\.", "\\-")) %>%
        data.table::data.table() %>%
        dplyr::left_join(my_annot, by = c('gene_name', 'transcript_name')) %>%
        dplyr::group_by(variable, gene_name) %>%
        dplyr::mutate(value.norm = value / sum(value),
                      wUTRlength.isof = value.norm * UTRlength,
                      wUTRlength.gene = mean(wUTRlength.isof)) %>%
        data.table::data.table() %>%
        dplyr::distinct(gene_tr, gene_name, transcript_name, variable,
                        value, UTRlength, wUTRlength.gene)

    message("Integrating Seurat metadata...(3/3)")
    my_expr = tibble::tibble(group=seurat_obj[[]][,group.by], variable=stringr::str_replace(colnames(seurat_obj),"\\.","\\-")) %>%
        dplyr::left_join(my_expr, by="variable") %>%
        dplyr::group_by(variable) %>%
        dplyr::mutate(wUTRlength.cell = mean(wUTRlength.gene)) %>%
        dplyr::distinct(variable, group, wUTRlength.cell) %>%
        dplyr::rename(cells=variable) %>%
        data.table::data.table()
    return(my_expr)
}
