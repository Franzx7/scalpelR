#' FindIsoforms
#'
#' This function identifies isoform expression with differential expression across conditions.
#'
#' @param seurat_obj The input scRNA-seq data
#' @param group.by The condition to compare isoforms
#' @param split.by description
#' @param assay Seurat Assay to use
#' @param threshold_pval Adjusted P.value threshold
#' @param threshold_abund Mininum Isoform abundance %
#' @param min.cells.fract Minimum fraction of cells expressing isoforms
#' @return A data frame with identified isoforms
#' @export
FindIsoforms <- function(
        seurat_obj,
        group.by = "orig.ident",
        split.by = NULL,
        assay = "RNA",
        threshold_pval = 0.05,
        threshold_abund = 0.10,
        min.cells.fract = 0.10) {

    #Perform checks
    #is group.by provided ?
    if(is.null(group.by)){
        stop("Error: Provide group.by !")}
    #is group.by %in% Seurat_obj ?
    if(!(group.by %in% colnames(seurat_obj@meta.data))){
        stop(paste0("Error: ", group.by, "not found in Seurat object !"))}
    #if group.by levels OK?
    if(length(unique(unlist(seurat_obj[[group.by]]))) == 1){
        stop(paste0(group.by, " contains an unique level ! ERROR"))}

    #is split.by %in% Seurat.obj ?
    if(!is.null(split.by)){
        if(!(split.by %in% colnames(seurat_obj@meta.data))){
            stop(paste0("Error: ", split.by, "not found in Seurat object !"))
        }else{
            message("splitting seurat object...")
            seurat_obj@meta.data = seurat_obj@meta.data %>%
                dplyr::rename(split=split.by)

            isoform_counts = lapply(unique(seurat_obj$split), function(sp_level){
                #a. subset Seurat object
                tmp_seurat = subset(seurat_obj, split==sp_level)
                #b. Processing
                message("Isoform filtering in cells...(1/5)")
                tmp_seurat@meta.data = seurat_obj@meta.data %>%
                    dplyr::rename('condition_defined'=group.by)

                isof_frac = lapply(unique(tmp_seurat$condition_defined), function(x){
                    tmp_obj = subset(tmp_seurat, condition_defined==x)
                    #get isoform expression
                    tmp_tab = reshape2::melt(data.frame(tmp_obj[["RNA"]]$counts) %>%
                                                 tibble::as_tibble(rownames = "gene_tr"), id.vars="gene_tr")
                    tmp_tab$cluster = x
                    tmp_tab = tmp_tab %>%
                        dplyr::group_by(gene_tr) %>%
                        dplyr::mutate(unexpressed=sum(value==0),
                                      expressed=sum(value>0),
                                      frac.expressed = expressed/(unexpressed+expressed)) %>%
                        dplyr::ungroup() %>%
                        dplyr::distinct(gene_tr,cluster,frac.expressed)
                }) %>%
                    data.table::rbindlist() %>%
                    dplyr::group_by(gene_tr) %>%
                    dplyr::reframe(check.frac=sum(frac.expressed<min.cells.fract),
                                   nclusters=dplyr::n_distinct(cluster)) %>%
                    dplyr::filter(check.frac!=nclusters)
                tmp_seurat = subset(tmp_seurat, features = isof_frac$gene_tr)

                message("Find genes with at least 2 isoforms...(2/5)")
                metainfos = tibble::tibble(gene_tr = rownames(tmp_seurat)) %>%
                    tidyr::separate(col = gene_tr,
                                    into = c("gene_name", "transcript_name"),
                                    remove = F, sep="\\*\\*\\*") %>%
                    dplyr::group_by(gene_name) %>%
                    dplyr::mutate(nb_trs = dplyr::n_distinct(transcript_name)) %>%
                    dplyr::filter(nb_trs>1) %>%
                    dplyr::distinct(gene_tr, gene_name, transcript_name)

                #check
                if(nrow(metainfos)==0) stop("No genes with > 1 isoforms... ! ERROR")

                message("Arrange selected isoforms in the defined conditions...(3/5)")
                isoforms_counts = Seurat::AggregateExpression(
                    subset(tmp_seurat, features = metainfos$gene_tr),
                    assays = assay,
                    group.by = 'condition_defined',
                    layers = "counts",
                    verbose = F,
                    return.seurat = F)[[assay]] %>%
                    data.frame() %>%
                    tibble::as_tibble(rownames = "gene_tr") %>%
                    dplyr::left_join(metainfos, by = "gene_tr") %>%
                    suppressMessages() %>%
                    suppressWarnings()

                message("Perform X2 test per gene across the defined condition...(4/5)")
                conds = colnames(dplyr::select(
                    isoforms_counts, !c(gene_name, transcript_name, gene_tr)))
                isoforms_counts = pbapply::pblapply(
                    split(isoforms_counts, isoforms_counts$gene_name),
                    function(x){
                        tmp = apply(x[,conds], 2, function(y) round(y/sum(y),2))
                        tmp = x[which(rowSums(tmp<threshold_abund)!=ncol(tmp)),]
                        tmp = tmp %>% dplyr::mutate_if(is.numeric, round)
                        if(nrow(tmp)>1){
                            #perform X2 statistical test
                            tmp.test = stats::chisq.test(tmp[,conds]) %>%
                                suppressWarnings()
                            tmp = tmp %>% dplyr::mutate(p_value = tmp.test$p.value)
                            return(tmp)
                        }else{
                            return(NULL)
                        }
                    }) %>% data.table::rbindlist()

                message("Pvalue adjustment...(5/5)")
                isoforms_counts = isoforms_counts %>%
                    dplyr::filter(p_value<threshold_pval) %>%
                    dplyr::mutate(
                        p_value.adjusted = stats::p.adjust(p_value, method="BH"), split=sp_level) %>%
                    dplyr::filter(p_value.adjusted < threshold_pval) %>%
                    dplyr::arrange(p_value.adjusted)

                #return
                return(isoforms_counts)
            }) %>% data.table::rbindlist()
            return(isoform_counts)
        }
    }else{
        message("Isoform filtering in cells...(1/5)")
        seurat_obj@meta.data = seurat_obj@meta.data %>%
            dplyr::rename(condition_defined=group.by)

        isof_frac = lapply(unique(seurat_obj$condition_defined), function(x){
            tmp_obj = subset(seurat_obj, condition_defined==x)
            #get isoform expression
            tmp_tab = reshape2::melt(data.frame(tmp_obj[["RNA"]]$counts) %>%
                                         tibble::as_tibble(rownames = "gene_tr"), id.vars="gene_tr")
            tmp_tab$cluster = x
            tmp_tab = tmp_tab %>%
                dplyr::group_by(gene_tr) %>%
                dplyr::mutate(unexpressed=sum(value==0),
                              expressed=sum(value>0),
                              frac.expressed = expressed/(unexpressed+expressed)) %>%
                dplyr::ungroup() %>%
                dplyr::distinct(gene_tr,cluster,frac.expressed)
        }) %>%
            data.table::rbindlist() %>%
            dplyr::group_by(gene_tr) %>%
            dplyr::reframe(check.frac=sum(frac.expressed<min.cells.fract),
                           nclusters=dplyr::n_distinct(cluster)) %>%
            dplyr::filter(check.frac!=nclusters)
        seurat_obj = subset(seurat_obj, features = isof_frac$gene_tr)

        message("Find genes with at least 2 isoforms...(2/5)")
        metainfos = tibble::tibble(gene_tr = rownames(seurat_obj)) %>%
            tidyr::separate(col = gene_tr,
                            into = c("gene_name", "transcript_name"),
                            remove = F, sep="\\*\\*\\*") %>%
            dplyr::group_by(gene_name) %>%
            dplyr::mutate(nb_trs = dplyr::n_distinct(transcript_name)) %>%
            dplyr::filter(nb_trs>1) %>%
            dplyr::distinct(gene_tr, gene_name, transcript_name)

        #check
        if(nrow(metainfos)==0) stop("No genes with > 1 isoforms... ! ERROR")

        message("Arrange selected isoforms in the defined conditions...(3/5)")
        isoforms_counts = Seurat::AggregateExpression(
            subset(seurat_obj, features = metainfos$gene_tr),
            assays = assay,
            group.by = 'condition_defined',
            layers = "counts",
            verbose = F,
            return.seurat = F)[[assay]] %>%
            data.frame() %>%
            tibble::as_tibble(rownames = "gene_tr") %>%
            dplyr::left_join(metainfos, by = "gene_tr") %>%
            suppressMessages() %>%
            suppressWarnings()

        message("Perform X2 test per gene across the defined condition...(4/5)")
        conds = colnames(dplyr::select(
            isoforms_counts, !c(gene_name, transcript_name, gene_tr)))
        isoforms_counts = pbapply::pblapply(
            split(isoforms_counts, isoforms_counts$gene_name),
            function(x){
                tmp = apply(x[,conds], 2, function(y) round(y/sum(y),2))
                tmp = x[which(rowSums(tmp<threshold_abund)!=ncol(tmp)),]
                tmp = tmp %>% dplyr::mutate_if(is.numeric, round)
                if(nrow(tmp)>1){
                    #perform X2 statistical test
                    tmp.test = stats::chisq.test(tmp[,conds]) %>%
                        suppressWarnings()
                    tmp = tmp %>% dplyr::mutate(p_value = tmp.test$p.value)
                    return(tmp)
                }else{
                    return(NULL)
                }
            }) %>% data.table::rbindlist()

        message("Pvalue adjustment...(5/5)")
        isoforms_counts = isoforms_counts %>%
            dplyr::filter(p_value<threshold_pval) %>%
            dplyr::mutate(
                p_value.adjusted = stats::p.adjust(p_value, method="BH")) %>%
            dplyr::filter(p_value.adjusted < threshold_pval) %>%
            dplyr::arrange(p_value.adjusted)

        #return
        return(isoforms_counts)
    }
}
