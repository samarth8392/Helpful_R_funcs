#' Visualize mutations on ideograms for specified genes from a MAF file
#'
#' @param maf_file Path to MAF file (character string)
#' @param genes Character vector of gene symbols to visualize
#' @param genome_build Genome build to use (default: "hg38")
#' @param output_dir Directory to save plots (NULL = don't save, just display)
#' @param region_padding Number of base pairs to add around gene region (default: 50,000)
#' @param show_all_genes Boolean to show all genes in region, not just target gene (default: TRUE)
#' @param height Plot height in inches when saving to PDF (default: 8)
#' @param width Plot width in inches when saving to PDF (default: 10)
#'
#' @return Invisibly returns a list of plot data
#' @export
#'
visualize_maf_ideogram <- function(maf_file, 
                                  genes, 
                                  genome_build = "hg38",
                                  output_dir = NULL,
                                  region_padding = 5e4,
                                  show_all_genes = TRUE,
                                  height = 8,
                                  width = 10) {
  
  # Check required packages
  required_packages <- c("Gviz", "GenomicRanges", "data.table", "utils")
  for(pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed. ",
                  "Please install it with: ", 
                  ifelse(pkg %in% c("Gviz", "GenomicRanges"),
                         "BiocManager::install('", 
                         "install.packages('"), 
                  pkg, "')"))
    }
  }
  
  # Check if TxDb is available
  txdb_pkg <- paste0("TxDb.Hsapiens.UCSC.", genome_build, ".knownGene")
  if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
    stop(paste0("Package '", txdb_pkg, "' is required but not installed. ",
                "Please install it with: BiocManager::install('", txdb_pkg, "')"))
  }
  
  # Load required packages without attaching
  txdb <- eval(parse(text = paste0(txdb_pkg, "::", txdb_pkg)))
  
  # Read MAF file
  message("Reading MAF file...")
  maf_data <- maf_file@data
  
  # Check if required columns exist
  required_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", 
                    "End_Position", "Variant_Classification", "Tumor_Sample_Barcode")
  missing_cols <- required_cols[!required_cols %in% colnames(maf_data)]
  if (length(missing_cols) > 0) {
    stop(paste0("MAF file is missing required columns: ", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Setup variant classification colors
  variant_colors <- c(
    "Missense_Mutation" = "red",
    "Nonsense_Mutation" = "blue",
    "Silent" = "grey",
    "Frame_Shift_Del" = "purple", 
    "Frame_Shift_Ins" = "orange",
    "In_Frame_Del" = "darkred",
    "In_Frame_Ins" = "darkgreen",
    "Splice_Site" = "darkblue",
    "Translation_Start_Site" = "gold",
    "Nonstop_Mutation" = "brown"
  )
  # Add default for any variant type not in our list
  setdiff_variants <- setdiff(unique(maf_data$Variant_Classification), names(variant_colors))
  if (length(setdiff_variants) > 0) {
    for (v in setdiff_variants) {
      variant_colors[v] <- "black"
    }
  }
  
  # Create output directory if needed
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # Process each gene
  result_list <- list()
  for (gene in genes) {
    message(paste0("Processing gene: ", gene))
    
    # Filter MAF for this gene
    gene_maf <- maf_data[maf_data$Hugo_Symbol == gene, ]
    
    if (nrow(gene_maf) == 0) {
      message(paste0("No mutations found for gene: ", gene))
      next
    }
    
    # Create GRanges object for mutations
    gr_gene <- GenomicRanges::GRanges(
      seqnames = paste0(gene_maf$Chromosome),
      ranges = IRanges::IRanges(
        start = gene_maf$Start_Position,
        end = gene_maf$End_Position
      ),
      gene = gene_maf$Hugo_Symbol,
      variant = gene_maf$Variant_Classification,
      sample = gene_maf$Tumor_Sample_Barcode
    )
    
    # Get chromosome and gene region boundaries
    chr <- as.character(GenomicRanges::seqnames(gr_gene)[1])
    gene_start <- max(min(GenomicRanges::start(gr_gene)) - region_padding, 1)
    gene_end <- min(max(GenomicRanges::end(gr_gene)) + region_padding, 3e9) # Safe upper limit
    
    # Create Gviz tracks
    message(paste0("Creating tracks for ", gene, " on ", chr))
    
    # 1. Ideogram track
    itrack <- Gviz::IdeogramTrack(genome = genome_build, chromosome = chr)
    
    # 2. Genome axis track
    gtrack <- Gviz::GenomeAxisTrack()
    
    # 3. Gene region track
    grtrack <- Gviz::GeneRegionTrack(
      txdb,
      chromosome = chr,
      start = gene_start,
      end = gene_end,
      name = "Genes",
      geneSymbol = TRUE,
      showId = TRUE
    )
    
    # If we only want to show the target gene, filter the track
    if (!show_all_genes) {
      grtrack <- Gviz::subset(grtrack, symbol == gene)
    }
    
    # 4. Annotation track for mutations (primary - boxes)
    colors_for_mut <- variant_colors[gr_gene$variant]
    colors_for_mut[is.na(colors_for_mut)] <- "black"
    
    atrack <- Gviz::AnnotationTrack(
      range = gr_gene,
      chromosome = chr,
      name = paste0(gene, " Mutations"),
      group = gr_gene$variant,
      stacking = "squish",
      shape = "box",
      fill = colors_for_mut,
      col = "black",
      alpha = 0.7,
      showFeatureId = FALSE,
      featureAnnotation = "group"
    )
    
    # 5. Additional Data Track for mutations (secondary - points)
    # This ensures mutations are visible as points even if boxes are hard to see
    dtrack <- Gviz::DataTrack(
      range = gr_gene,
      chromosome = chr,
      name = "Mutation Points",
      data = rep(1, length(gr_gene)),  # all points at y=1
      type = "p",
      cex = 1.2,
      pch = 19,  # solid circle
      col = colors_for_mut
    )
    
    # 6. Sample distribution track (shows which samples have mutations)
    # Get counts per sample
    sample_counts <- table(gr_gene$sample)
    unique_samples <- names(sample_counts)
    
    # Create sample ranges
    sample_gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(
        start = rep(gene_start, length(unique_samples)),
        end = rep(gene_end, length(unique_samples))
      ),
      sample = unique_samples,
      count = as.numeric(sample_counts)
    )
    
    # Sample bar chart
    strack <- Gviz::AnnotationTrack(
      range = sample_gr,
      chromosome = chr,
      name = "Samples",
      group = sample_gr$sample,
      stacking = "squish",
      fill = rainbow(length(unique_samples)),
      background.title = "darkblue"
    )

    sample_levels <- unique(gr_gene$sample)
    sample_map <- setNames(seq_along(sample_levels), sample_levels)
    gr_gene$y <- sample_map[gr_gene$sample]

    ggstyle_dtrack <- Gviz::DataTrack(
      start = GenomicRanges::start(gr_gene),
      end = GenomicRanges::end(gr_gene),
      chromosome = chr,
      genome = genome_build,
      name = "Mutations (by Sample)",
      type = "p",
      data = gr_gene$y,
      pch = 21,
      cex = 1.5,
      col = "black",  # border color
      fill = variant_colors[gr_gene$variant],  # fill color per point
      legend = FALSE  # suppress legend since grouping isn't used here
    )
    
    # Plot tracks
    if (!is.null(output_dir)) {
  pdf_file <- file.path(output_dir, paste0(gene, "_mutation_ideogram.pdf"))
  message(paste0("Saving plot to: ", pdf_file))
  grDevices::pdf(pdf_file, width = width, height = height * 1.5)  # Increased height to fit two panels
}

# Panel 1: Ideogram Plot (capture it in a grid object)
gviz_plot <- grid::grid.grabExpr({
  Gviz::plotTracks(
    list(itrack, gtrack, grtrack, atrack, ggstyle_dtrack),
    from = gene_start,
    to = gene_end,
    chromosome = chr,
    background.title = "darkgrey",
    col.axis = "black",
    col.title = "black",
    cex.title = 0.8,
    margin = 20,
    main = paste0("Mutations on ", gene, " locus")
  )
})

plot_df <- data.frame(
  sample = gr_gene$sample,
  start = GenomicRanges::start(gr_gene),
  variant_class = gr_gene$variant,
  variant_type = ifelse(width(gr_gene) == 1, "SNP", "INDEL")
)

# Order samples for vertical axis
plot_df$sample <- factor(plot_df$sample, levels = rev(unique(plot_df$sample)))

# Plot
gg_mut_plot <- ggplot(plot_df, aes(x = start, y = sample, color = variant_class, shape = variant_type)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = variant_colors) +
  labs(
    title = paste0("Mutations in ", gene, " on Chromosome ", chr),
    x = "Genomic Position", y = "Sample"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 6)
  )
    gridExtra::grid.arrange(gviz_plot, gg_mut_plot, ncol = 1, heights = c(2, 1))

    if (!is.null(output_dir)) {
      dev.off()
    }
    
    # Store result data
    result_list[[gene]] <- list(
      gene = gene,
      chromosome = chr,
      start = gene_start,
      end = gene_end,
      mutations = gr_gene,
      samples = unique_samples,
      sample_counts = sample_counts
    )
  }
  
  message("Completed gene mutation ideogram visualization")
  invisible(result_list)
}
