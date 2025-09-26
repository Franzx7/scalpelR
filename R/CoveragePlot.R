#' Plot Read Coverage and Gene/Transcript Structure from BAM Files
#'
#' Generates a visualization of gene and transcript annotation alongside read coverage from one or more BAM files using the Gviz package.
#'
#' @param genome_gr A GRanges object containing genome annotation (e.g., from GTF).
#' @param bamfiles A named vector of BAM file paths. Names will be used as track labels.
#' @param gene_in The gene name of interest (must be present in \code{genome_gr}).
#' @param features Optional. A vector of transcript names to be displayed for the gene (\code{NULL} for all transcripts).
#' @param genome_sp The genome build identifier, either \code{"mm10"} or \code{"hg38"}. Default is \code{"mm10"}.
#' @param distZOOM Optional. Zoom to a subregion (numeric, in bp) from gene boundary. If \code{NULL}, plot the entire gene region.
#' @param annot_tab Optional. Annotation table. Not used internally but reserved for future use.
#' @param filter_trs Logical. Whether to filter transcripts. Default is \code{FALSE}.
#' @param extend.left Numeric. Number of bases to extend to the left of the gene region. Default is 1000.
#' @param extend.right Numeric. Number of bases to extend to the right of the gene region. Default is 1000.
#' @param samtools.bin Character. Path to the \code{samtools} executable. Default is \code{"samtools"}.
#'
#' @details
#' For each BAM file, coverage is computed over the gene region (optionally including selected transcripts and extensions)
#' using \code{samtools view} and \code{samtools depth}. The gene region and transcript structure are shown using Gviz tracks.
#' The function can zoom to either end of the gene, depending on strand and the \code{distZOOM} parameter.
#'
#' @return
#' Invisibly returns the Gviz plot object (plot is drawn as a side effect).
#'
#' @examples
#' \dontrun{
#' CoveragePlot(
#'   genome_gr = gtf_gr,
#'   bamfiles = c(Sample1 = "sample1.bam", Sample2 = "sample2.bam"),
#'   gene_in = "GeneA",
#'   features = c("Transcript1", "Transcript2"),
#'   genome_sp = "mm10",
#'   extend.left = 2000,
#'   extend.right = 2000
#' )
#' }
#'
#' @export
CoveragePlot <- function(
    genome_gr,
    bamfiles,
    gene_in,
    features = NULL,
    genome_sp = "mm10",
    distZOOM = NULL,
    extend.left = 1000,
    extend.right = 1000,
    samtools.bin = "samtools") {
  # Function to plot read coverage on genome

  # checks
  if (!genome_sp %in% c("mm10", "hg38")) {
    stop("! Enter genome species: mm10 or hg38")
  }

  # resume input args...
  inputs_in <- list(
    "BAM file:" = bamfiles,
    "Gene:" = gene_in,
    "Transcripts:" = features,
    "Genome annotation org:" = genome_sp
  )
  print(inputs_in)
  options(ucscChromosomeNames = FALSE)

  # 0. set axis track
  axisTrack <- Gviz::GenomeAxisTrack(genome = genome_sp, name = gene_in, cex = 1, col = "black")

  # 1. Set Gene Tracks
  message("Building Gene Track...")
  gene.tab <- data.frame(genome_gr[genome_gr$gene_name == gene_in]) %>%
    tibble::tibble() %>%
    dplyr::filter(type %in% c("UTR", "CDS", "exon")) %>%
    dplyr::distinct(seqnames, start, end, width, strand, gene_id, exon_id, transcript_id, gene_name, transcript_name, type) %>%
    dplyr::rename(Chromosome = "seqnames", feature = "type", gene = "gene_id", transcript = "transcript_id") %>%
    dplyr::group_by(transcript_name) %>%
    dplyr::mutate(check = dplyr::if_else(any(feature == "UTR"), "wUTR", "noUTR")) %>%
    dplyr::filter(!(feature == "exon" & check == "wUTR")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(feature = dplyr::if_else(feature == "UTR", "utr", feature), transcript = transcript_name) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::arrange(transcript_name, start)

  if (!is.null(features)) {
    gene.tab <- dplyr::filter(gene.tab, transcript_name %in% features)
  }

  jgroup <- NULL
  if (gene.tab$strand[1] == "-") {
    jgroup <- "left"
  } else {
    jgroup <- "right"
  }
  gene.track <- Gviz::GeneRegionTrack(GenomicRanges::makeGRangesFromDataFrame(gene.tab, keep.extra.columns = T),
    chromosome = gene.tab$Chromosome[1], name = gene_in, transcriptAnnotation = "transcript",
    thinBoxFeature = c("utr", "exon"), just.group = jgroup, genome = genome_sp, fontcolor.group = "black",
    fill = "#38B6A3", color = "black", col.title = "black", col.line = "black",
    lwd = 3, cex.title = 2, fontsize.group = 25, background.title = "cornflowerblue", size = 0.3
  )

  # 2. set Data track
  message("Building Alignment Track...")
  res <- lapply(1:length(bamfiles), function(x) {
    # get coverage
    data.table::fwrite(dplyr::distinct(gene.tab, Chromosome, start, end), file = "./coords.txt", col.names = F, sep = "\t")
    system(paste0(samtools.bin, " view -b --region-file ./coords.txt ", bamfiles[x], " > current.bam"))
    cov.exp <- system(paste0(samtools.bin, " depth -b coords.txt current.bam > current.cov"))
    cov.tab <- data.table::fread("current.cov", col.names = c("seqnames", "start", "depth")) %>% dplyr::filter(depth >= 0)
    # dataTrack
    curr.track <- Gviz::DataTrack(
      start = cov.tab$start, width = 1,
      data = cov.tab$depth, chromosome = gene.tab$Chromosome[1], genome = genome_sp,
      cex.title = 2, cex = 4, cex.axis = .9, col.line = "black",
      type = c("hist"), background.title = "gray20", name = names(bamfiles)[x], col.histogram = "#ede32b"
    )
    return(list(max.depth = max(cov.tab$depth), track = curr.track))
  })
  data.track <- lapply(res, function(x) x[[2]])
  YMAX <- max(unlist(lapply(res, function(x) x[[1]])))
  system("rm current.cov")
  system("rm current.bam")

  # Plot
  if (gene.tab$strand[1] == "+") {
    if (is.null(distZOOM)) {
      Gviz::plotTracks(c(axisTrack, gene.track, data.track),
        from = min(gene.tab$start), to = max(gene.tab$end) + extend.right,
        sizes = c(0.1, 1 / 3, rep(1 / 4, length(bamfiles))),
        ylim = c(0, YMAX),
        fig.width = 20, fig.height = 5
      )
    } else {
      Gviz::plotTracks(c(axisTrack, gene.track, data.track),
        from = max(gene.tab$end) - distZOOM, to = max(gene.tab$end) + extend.right,
        sizes = c(0.1, 1 / 3, rep(1 / 4, length(bamfiles))),
        ylim = c(0, YMAX),
        fig.width = 20, fig.height = 5
      )
    }
  } else {
    if (is.null(distZOOM)) {
      Gviz::plotTracks(c(axisTrack, gene.track, data.track),
        from = min(gene.tab$start) - extend.left, to = max(gene.tab$end),
        sizes = c(0.1, 1 / 3, rep(1 / 4, length(bamfiles))),
        ylim = c(0, YMAX),
        fig.width = 20, fig.height = 5
      )
    } else {
      Gviz::plotTracks(c(axisTrack, gene.track, data.track),
        from = min(gene.tab$start) - extend.left, to = min(gene.tab$start) + distZOOM,
        sizes = c(0.1, 1 / 3, rep(1 / 4, length(bamfiles))),
        ylim = c(0, YMAX),
        fig.width = 20, fig.height = 5
      )
    }
  }
}
