#' Identify Differentially Expressed Isoforms Across Conditions
#'
#' Identifies isoform expression with differential expression across specified conditions in single-cell RNA-seq data using a Seurat object.
#'
#' @param seurat_obj The input Seurat object containing scRNA-seq data.
#' @param group_by The metadata variable used to group cells for isoform expression comparison.
#' @param split_by Optional. Metadata variable used to split the Seurat object before analysis. If \code{NULL}, no split is performed. Default is \code{NULL}.
#' @param assay Character. The Seurat assay to use for analysis. Default is \code{"RNA"}.
#' @param threshold_pval Numeric. Adjusted p-value threshold for significance. Default is \code{0.05}.
#' @param threshold_abund Numeric. Minimum isoform abundance proportion per group. Default is \code{0.1}.
#' @param threshold_sd Numeric. Minimum standard deviation of counts across groups for an isoform to be tested. Default is \code{0.05}.
#'
#' @details
#' The function aggregates isoform (transcript) expression per group, filters isoforms based on abundance and variability, and applies a chi-squared test to assess differential isoform usage. Adjusted p-values (Benjamini-Hochberg) are reported. If \code{split_by} is specified and different from \code{group_by}, the analysis is performed separately for each split.
#'
#' @return
#' A data.frame (or data.table) containing isoforms passing the significance threshold, with statistics for differential expression. If splitting, results from each split are combined.
#'
#' @examples
#' \dontrun{
#' FindIsoforms(
#'   seurat_obj = my_seurat,
#'   group_by = "condition",
#'   split_by = NULL,
#'   assay = "RNA",
#'   threshold_pval = 0.05,
#'   threshold_abund = 0.1,
#'   threshold_sd = 0.05
#' )
#' }
#'
#' @export
FindIsoforms <- function(
    seurat_obj,
    group_by,
    split_by = NULL,
    assay = "RNA",
    threshold_pval = 0.05,
    threshold_abund = 0.1,
    threshold_sd = 0.05) {
  # check split.by
  if (!is.null(split_by)) {
    if (split_by != group_by) {
      message("Splitting Seurat object...")
      seurat_obj <- Seurat::SplitObject(seurat_obj, split.by = split_by)
    } else {
      stop("split_by & group_by args are identical !")
    }
  }

  # custom functions
  extract_aggtab <- function(seurat_obj, group, assay) {
    # Get aggregated counts table
    tmp <- Seurat::AggregateExpression(
      seurat_obj,
      group.by = group,
      assays = assay
    )[[1]]

    col_names <- colnames(tmp)

    tmp <- tmp %>%
      tibble::as_tibble(rownames = "gene_tr") %>%
      tidyr::separate(
        gene_tr,
        c("gene_name", "transcript_name"),
        sep = "\\*\\*\\*",
        remove = FALSE
      ) %>%
      dplyr::arrange(gene_name, transcript_name)

    # Compute number of isoforms per gene
    tmp_stats <- dplyr::distinct(tmp, gene_name, transcript_name) %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarise(nb_isoforms = dplyr::n_distinct(transcript_name))

    tmp <- dplyr::filter(tmp, gene_name %in% tmp_stats$gene_name) %>%
      data.table::as.data.table()

    colnames(tmp) <- c("gene_tr", "gene_name", "transcript_name", col_names)

    tmp <- tmp %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)))

    return(tmp)
  }

  processing_tab <- function(seurat.obj, threshold.pval, threshold.abund, threshold.sd) {
    # Aggregate counts
    message("Aggregating counts in conditions...")
    counts_tab <- extract_aggtab(seurat.obj, group = group_by, assay = assay) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.))) %>%
      dplyr::filter(!if_all(where(is.numeric), ~ . < 5)) %>%
      dplyr::group_by(gene_name) %>%
      dplyr::filter(dplyr::if_all(dplyr::where(is.numeric), ~ (. / sum(.)) > threshold.abund)) %>%
      dplyr::filter(dplyr::if_all(dplyr::where(is.numeric), ~ stats::sd(.) > threshold.sd)) %>%
      dplyr::filter(dplyr::n() >= 2) %>%
      ungroup()

    # Filtering and tests
    message("Performing filtering and comparison test...")
    all_results <- pbapply::pblapply(
      split(counts_tab, counts_tab$gene_name),
      function(gene_tab) {
        chi2_res <- stats::chisq.test(
          as.matrix(gene_tab %>% dplyr::select(dplyr::where(is.numeric))),
          simulate.p.value = T
        ) %>%
          suppressWarnings()
        return(gene_tab %>% dplyr::mutate(p_value = chi2_res$p.value, statistic = chi2_res$statistic))
      }
    ) %>%
      data.table::rbindlist(fill = T) %>%
      dplyr::filter(!is.na(p_value))
    adjusted.pvals <- stats::p.adjust(all_results$p_value, method = "BH")
    all_results$p_value_adjusted <- adjusted.pvals
    all_results <- dplyr::filter(all_results, p_value_adjusted < threshold.pval) %>%
      dplyr::arrange(p_value_adjusted, gene_tr)

    return(all_results)
  }

  # __main__

  if (is.null(split_by)) {
    return(processing_tab(seurat_obj, threshold_pval, threshold_abund, threshold_sd))
  } else {
    return(lapply(names(seurat_obj), function(x) {
      res <- processing_tab(seurat_obj[[x]], threshold_pval, threshold_abund, threshold_sd)
      res$split_by <- x
      return(res)
    }) %>% data.table::rbindlist())
  }
}
