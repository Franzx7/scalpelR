#' Coverage Plot for Gene and Transcript Regions using BAM files
#'
#' This function generates a gene region track and corresponding read coverage plots
#' from BAM files using the Gviz package.
#'
#' @param bamfiles A named vector of BAM file paths. Names will be used as labels.
#' @param gtf A GRanges object containing the GTF annotation.
#' @param gene_name The gene name of interest (must be in GTF).
#' @param transcript_names A vector of transcript names to be plotted.
#' @param extend.range Numeric value specifying the extension of the plotting region
#'        beyond the gene's boundaries. Default is 0 (no extension).
#' @param genome The genome build identifier (e.g., "mm10", "hg19"). Default is "mm10".
#' @param zoom Logical. Should the plot zoom into a particular region? Default is FALSE.
#' @param zoom_focus Numeric. Percentage of the gene region to focus on when zooming (default: 0.1).
#' @return A Gviz plot displaying the gene and transcript annotations and coverage tracks from BAM files.
#' @examples
#' CoveragePlot(bamfiles = c("sample1.bam", "sample2.bam"),
#'              gtf = gtf_data,
#'              gene_name = "GeneA",
#'              transcript_names = c("Transcript1", "Transcript2"),
#'              extend.range = 1000)
#' @export
CoveragePlot <- function(bamfiles,
                         gtf,
                         gene_name,
                         transcript_names,
                         extend.range = 0,
                         genome = "mm10",
                         zoom = FALSE,
                         zoom_focus = 0.1) {

    options(ucscChromosomeNames=FALSE)

    # Validate inputs
    if (!is.vector(bamfiles) || length(bamfiles) == 0) {
        stop("Please provide a named vector of BAM file paths.")
    }
    if (!gene_name %in% gtf$gene_name) {
        stop(paste0("Gene name '", gene_name, "' not found in GTF annotation."))
    }
    if (length(transcript_names) == 0) {
        stop("Please provide a vector of transcript names.")
    }
    if (!all(file.exists(bamfiles))) {
        missing_bams <- bamfiles[!file.exists(bamfiles)]
        stop(paste("The following BAM files were not found:", paste(missing_bams, collapse = ", ")))
    }

    # Subset the GTF annotation to get the gene and transcripts of interest
    message('Subsetting the GTF annotation... (1/4)')
    gene_gtf <- gtf[gtf$gene_name == gene_name &
                        gtf$type %in% c('UTR', 'CDS', 'exon') &
                        gtf$transcript_name %in% transcript_names]

    # Modify annotations for easier plotting
    gene_gtf$transcript <- paste(gene_gtf$gene_name, gene_gtf$transcript_name, sep = " - ")
    gene_gtf$gene <- gene_gtf$gene_id
    gene_gtf$feature <- stringr::str_replace(gene_gtf$type, 'UTR', 'utr')

    # Check for UTR presence and modify exon handling
    gene_gtf$check <- lapply(GenomicRanges::split(gene_gtf, gene_gtf$transcript), function(x) {
        if ('UTR' %in% x$type) return(TRUE) else return(FALSE)
    }) %>% unlist()

    # Remove redundant exons in transcripts with UTR
    gene_gtf <- gene_gtf[!(gene_gtf$type == "exon" & gene_gtf$check == TRUE)]

    if (length(gene_gtf) == 0) {
        stop(paste0("No matching transcripts found for the gene '", gene_name, "'."))
    }

    # Define the genomic region for plotting
    plot_region <- range(gene_gtf)

    # Handle zooming logic
    if (zoom) {
        dist_zoom <- round(GenomicRanges::width(plot_region) * zoom_focus)  # Default zoom to 10% of the region size
    }

    # Create the GeneRegionTrack for plotting
    message('Creating the Gene Region Track... (2/4)')
    chromosome_in <- as.character(GenomicRanges::seqnames(gene_gtf))[1]
    strand_in <- as.character(GenomicRanges::strand(gene_gtf))[1]
    gene_name.pos <- ifelse(strand_in == "-", 'left', 'right')

    gene_track <- Gviz::GeneRegionTrack(gene_gtf, genome = genome,
                                        chromosome = chromosome_in,
                                        name = gene_name,
                                        transcriptAnnotation = "transcript",
                                        fill = "orange", col = "black",
                                        background.title = "black",
                                        fontsize.group = 20,
                                        just.group = 'below')

    # Create alignment tracks for each BAM file
    message('Creating alignment tracks for BAM files... (3/4)')
    cov_tracks <- lapply(names(bamfiles), function(bam) {
        Gviz::AlignmentsTrack(bamfiles[[bam]], genome = genome,
                              isPaired = TRUE, name = bam,
                              chromosome = chromosome_in,
                              type = c("coverage"),
                              col = "black", fill = "blue",
                              background.title = "coral4")
    })

    # Plot the gene and coverage tracks
    message('Plotting the gene and coverage tracks... (4/4)')
    if (strand_in == "+") {
        if (!zoom) {
            Gviz::plotTracks(
                c(gene_track, cov_tracks),
                sizes = c(1/3, rep(2/3 / length(cov_tracks), length(cov_tracks))),
                from = GenomicRanges::start(plot_region),
                to = GenomicRanges::end(plot_region) + extend.range,
                transcriptAnnotation = "transcript",
                main = paste("Coverage Plot for", gene_name)
            )
        } else {
            Gviz::plotTracks(
                c(gene_track, cov_tracks),
                sizes = c(1/3, rep(2/3 / length(cov_tracks), length(cov_tracks))),
                from = GenomicRanges::end(plot_region) - dist_zoom,
                to = GenomicRanges::end(plot_region) + extend.range,
                transcriptAnnotation = "transcript",
                main = paste("Coverage Plot for", gene_name)
            )
        }
    } else {
        if (!zoom) {
            Gviz::plotTracks(
                c(gene_track, cov_tracks),
                sizes = c(1/3, rep(2/3 / length(cov_tracks), length(cov_tracks))),
                from = GenomicRanges::start(plot_region) - extend.range,
                to = GenomicRanges::end(plot_region),
                transcriptAnnotation = "transcript",
                main = paste("Coverage Plot for", gene_name)
            )
        } else {
            Gviz::plotTracks(
                c(gene_track, cov_tracks),
                sizes = c(1/3, rep(2/3 / length(cov_tracks), length(cov_tracks))),
                from = GenomicRanges::start(plot_region) - extend.range,
                to = GenomicRanges::start(plot_region) + dist_zoom,
                transcriptAnnotation = "transcript",
                main = paste("Coverage Plot for", gene_name)
            )
        }
    }
}
