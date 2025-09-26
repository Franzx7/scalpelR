#' Plot Relative Isoform or Feature Expression Across Groups
#'
#' Generates boxplots of relative expression for selected features (e.g., isoforms or genes) across groups in a Seurat object.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param features Character vector of features (e.g., genes or isoforms) to plot.
#' @param group_var The metadata column in the Seurat object to group cells by.
#' @param plot Character. Whether to plot by group (\code{"group"}) or by features (\code{"features"}). Default is \code{"group"}.
#' @param assay Character. The Seurat assay to use. Default is \code{"RNA"}.
#' @param log_normalization Logical. Whether to log-normalize the expression values (\code{TRUE}) or plot raw counts (\code{FALSE}). Default is \code{TRUE}.
#'
#' @details
#' The function extracts counts for the specified features from the specified assay in the Seurat object. It plots boxplots of expression values, grouped either by metadata group or by features, and can display either raw or log-normalized counts.
#'
#' @return
#' A ggplot2 object showing boxplots of expression values for the selected features across the specified groups.
#'
#' @examples
#' \dontrun{
#' plotExpression(
#'   seurat_obj = my_seurat_obj,
#'   features = c("Isoform1", "Isoform2"),
#'   group_var = "Cluster",
#'   plot = "group",
#'   assay = "RNA",
#'   log_normalization = TRUE
#' )
#' }
#'
#' @export
plotExpression <- function(
    seurat_obj,
    features,
    group_var,
    plot = "group",
    assay = "RNA",
    layer = "counts",
    log_normalization = TRUE,
    only_expressed_cells = FALSE) {
  # Function to plot with geom_boxplot relative isoform expression

  if (!plot %in% c("group", "features")) {
    stop("Error plot args: group or features")
  }

  # extract meta.data
  meta.infos <- data.table::data.table(group = seurat_obj@meta.data[, group_var], cells = colnames(seurat_obj))

  # compute counts tab
  tab <- seurat_obj[[assay]][layer][features, ] %>%
    as_tibble(rownames = "features") %>%
    reshape2::melt(variable.name = "cells") %>%
    dplyr::left_join(meta.infos, by = "cells") %>%
    dplyr::mutate(value.logNorm = log1p(value)) %>%
    dplyr::ungroup()

  # discard all the cells not expressed cells if required
  if (only_expressed_cells) {
    tab <- dplyr::filter(tab, value != 0)
  }

  # plot
  p1 <- ggplot2::ggplot(tab, ggplot2::aes(group, value, fill = features, color = features)) +
    ggplot2::xlab(group_var) +
    ggplot2::ylab("counts") +
    ggplot2::geom_violin(width = .2, trim = T, scale = "width", position = ggplot2::position_dodge(width = 0.3), alpha = .6) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge(width = 0.3), alpha = .5) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.3), width = .05, color = "black", outlier.shape = NA, size = .2, outlier.size = .2) +
    ggplot2::ggtitle(paste0("Features expression in ", group_var)) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(legend.position = "top", axis.text = ggplot2::element_text(size = 17, face = "bold"))

  p2 <- ggplot2::ggplot(tab, ggplot2::aes(group, value.logNorm, fill = features, color = features)) +
    ggplot2::xlab(group_var) +
    ggplot2::ylab("counts") +
    ggplot2::geom_violin(width = .2, trim = T, scale = "width", position = ggplot2::position_dodge(width = 0.3), alpha = .7, color = "black", size = .2) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge(width = 0.3), alpha = .5) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.3), width = .05, color = "black", outlier.shape = NA, size = .2, outlier.size = .2) +
    ggplot2::ggtitle(paste0("Features expression in ", group_var)) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(legend.position = "top", axis.text = ggplot2::element_text(size = 17, face = "bold"))

  p3 <- ggplot2::ggplot(tab, ggplot2::aes(features, value, fill = group, color = group)) +
    ggplot2::xlab(group_var) +
    ggplot2::ylab("counts") +
    ggplot2::geom_violin(width = .2, trim = T, scale = "width", position = ggplot2::position_dodge(width = 0.3), alpha = .7, color = "black", size = .2) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge(width = 0.3), alpha = .5) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.3), width = .05, color = "black", outlier.shape = NA, size = .2, outlier.size = .2) +
    ggplot2::ggtitle(paste0("Features expression in ", group_var)) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(legend.position = "top", axis.text = ggplot2::element_text(size = 17, face = "bold"))

  p4 <- ggplot2::ggplot(tab, ggplot2::aes(features, value.logNorm, fill = group, color = group)) +
    ggplot2::xlab(group_var) +
    ggplot2::ylab("counts") +
    ggplot2::geom_violin(width = .2, trim = T, scale = "width", position = ggplot2::position_dodge(width = 0.3), alpha = .7, color = "black", size = .2) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge(width = 0.3), alpha = .5) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.3), width = .05, color = "black", outlier.shape = NA, size = .2, outlier.size = .2) +
    ggplot2::ggtitle(paste0("Features expression in ", group_var)) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(legend.position = "top", axis.text = ggplot2::element_text(size = 17, face = "bold"))

  # return
  if (log_normalization & plot == "group") {
    return(p2)
  } else if (log_normalization & plot == "features") {
    return(p4)
  } else if (!log_normalization & plot == "group") {
    return(p1)
  } else if (!log_normalization & plot == "features") {
    return(p3)
  }
}
