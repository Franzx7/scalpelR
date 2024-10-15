#' Plot Differential Isoform Usage (DIU)
#'
#' This function visualizes isoform usage across groups from a Seurat object
#' based on differential isoform usage (DIU) data. It calculates isoform
#' fractions and generates a violin plot of the isoform expression across clusters.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param DIU_table A data frame containing differential isoform usage information.
#' It must include a column `gene_name` and `gene_tr`.
#' @param group.by The metadata column in the Seurat object to group cells by.
#' @param gene The gene for which isoform usage will be plotted.
#' @param extended To add the gene expression plot for each cluster (T) or not(F)
#' @param type To show the fraction of Isoform or normalized expression in clusters
#'
#' @return A ggplot object showing the isoform fractions in the specified groups.
#'
#' @examples
#' \dontrun{
#' plotDIU(seurat_obj = my_seurat_obj, DIU_table = my_DIU_table, group.by = "Cluster", gene = "GeneA")
#' }
#'
#' @export
plotDIU <- function(
        seurat_obj,
        DIU_table,
        group.by = NULL,
        gene,
        extended=F,
        type='fraction') {

    # Checks
    if (is.null(group.by)) stop("Provide group.by argument!")
    if (!(group.by %in% colnames(seurat_obj@meta.data))) {
        stop("group.by argument not found in the provided Seurat object!")
    }

    # Extract differential isoforms from DIU table (a)
    tmp_tab <- dplyr::filter(DIU_table, gene_name == gene)
    seurat_obj$Barcodes <- stringr::str_replace(colnames(seurat_obj), "\\.", "\\-")

    # Extract Differential Isoform expression from Seurat object (b)
    tmp_tab <- subset(seurat_obj, features = tmp_tab$gene_tr)[["RNA"]]$counts %>%
        data.frame() %>%
        tibble::as_tibble(rownames = "gene_tr") %>%
        reshape2::melt(id.vars = "gene_tr") %>%
        dplyr::rename(Barcodes = 'variable') %>%
        dplyr::mutate(Barcodes = stringr::str_replace(Barcodes, "\\.", "\\-")) %>%
        data.table::data.table() %>%
        dplyr::left_join(
            dplyr::distinct(seurat_obj[[]][, c('Barcodes', group.by)]),
            by = "Barcodes"
        ) %>%
        dplyr::filter(value>0) %>%
        dplyr::group_by(Barcodes) %>%
        dplyr::mutate(IF = value / sum(value), value.tot=sum(value)) %>%
        data.table::data.table() %>%
        tidyr::separate(col='gene_tr', into=c('gene_name','transcript_name'),
                        remove=F) %>%
        dplyr::mutate(value.log=log(value),value.tot.log=log(value.tot))


    # Plotting (c)
    OUT.extended = ggplot2::ggplot(tmp_tab, ggplot2::aes_string(group.by, 'value.tot.log',
                                                                color=group.by, fill=group.by)) +
        ggplot2::geom_point(size=3, color='gray10') +
        ggplot2::geom_violin(color = "black", alpha = 0.5, trim = T) +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::facet_wrap(~gene_name) +
        ggplot2::ylab("Log_normalized expression") +
        ggplot2::ggtitle(paste0(gene, ' expression in defined cluster'))


    OUT.a = ggplot2::ggplot(tmp_tab, ggplot2::aes_string(group.by, 'value.log',
                                                         color = group.by,
                                                         fill = group.by)) +
        ggplot2::geom_point(size = 3, color = "gray10") +
        ggplot2::geom_violin(color = "black", alpha = 0.5, trim = FALSE) +
        ggplot2::facet_wrap(~gene_tr) +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::ylab("Isoform expresion (log_normalized) ") +
        ggplot2::ggtitle("Isoform expression in defined cluster for each cell") %>%
        suppressWarnings()


    OUT.b = ggplot2::ggplot(tmp_tab, ggplot2::aes_string(group.by, 'IF',
                                                         color = group.by,
                                                         fill = group.by)) +
        ggplot2::geom_point(size = 3, color = "gray10") +
        ggplot2::geom_violin(color = "black", alpha = 0.5, trim = FALSE) +
        ggplot2::facet_wrap(~gene_tr) +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::ylab("Isoform fraction ") +
        ggplot2::ggtitle("Isoform fraction in defined cluster for each cell") %>%
        suppressWarnings()

    #return
    if(type=="expression"){
        if(extended){OUT.extended + OUT.a}else{OUT.a}
    }else if(type=="fraction"){
        if(extended){OUT.extended + OUT.b}else{OUT.b}
    }
}
