#' Visualize mutations on a grid for specified genes from a MAF file
#'
#' @param maf_input Either a MAF object from maftools or a data frame with MAF format data
#' @param genes Character vector of gene symbols to visualize
#' @param output_dir Directory to save plots (NULL = don't save, just display)
#' @param height Plot height in inches when saving to PDF (default: 8)
#' @param width Plot width in inches when saving to PDF (default: 10)
#' @param max_genes_per_row Maximum number of genes to display in each row of the grid (default: 3)
#' @param text_size Size of the HGVSp label text (default: 2.5)
#' @param point_size Size of mutation points (default: 3)
#'
#' @return Invisibly returns the ggplot object
#' @export
#'
visualize_maf_aa_grid <- function(maf_input, 
                                 genes, 
                                 output_dir = NULL,
                                 height = 8,
                                 width = 10,
                                 max_genes_per_row = 3,
                                 text_size = 2.5,
                                 point_size = 3) {
  
  # Check required packages
  required_packages <- c("ggplot2", "data.table", "gridExtra", "dplyr", "stringr", "maftools", "ggrepel")
  for(pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed. ",
                  "Please install it with: install.packages('", pkg, "')"))
    }
  }
  
  # Access the MAF data
  message("Reading MAF data...")
  if (inherits(maf_input, "MAF")) {
    # If input is a MAF object from maftools
    maf_data <- maf_input@data
  } else if (is.data.frame(maf_input)) {
    # If input is already a data frame
    maf_data <- maf_input
  } else {
    stop("Input must be either a MAF object from maftools or a data frame in MAF format")
  }
  
  # Check if required columns exist
  required_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "HGVSp_Short")
  missing_cols <- required_cols[!required_cols %in% colnames(maf_data)]
  if (length(missing_cols) > 0) {
    stop(paste0("MAF data is missing required columns: ", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Setup custom variant classification colors
  variant_colors <- c(
    "Missense_Mutation" = "#FF5733",     # Orange-red
    "Nonsense_Mutation" = "#3498DB",     # Blue
    "Silent" = "#95A5A6",                # Gray
    "Frame_Shift_Del" = "#9B59B6",       # Purple
    "Frame_Shift_Ins" = "#F1C40F",       # Yellow
    "In_Frame_Del" = "#E74C3C",          # Red
    "In_Frame_Ins" = "#2ECC71",          # Green
    "Splice_Site" = "#34495E",           # Dark blue
    "Translation_Start_Site" = "#D4AC0D", # Gold
    "Nonstop_Mutation" = "#A04000"       # Brown
  )
  
  # Create output directory if needed
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # Filter MAF for the specified genes
  gene_maf <- maf_data[maf_data$Hugo_Symbol %in% genes, ]
  
  if (nrow(gene_maf) == 0) {
    stop("No mutations found for any of the specified genes")
  }
  
  # Extract amino acid positions from HGVSp_Short
  gene_maf$aa_pos <- NA
  
  # Extract positions using regex
  pattern <- "p\\.[A-Za-z]+([0-9]+)[A-Za-z\\*]+"
  matches <- regmatches(gene_maf$HGVSp_Short, 
                        regexec(pattern, gene_maf$HGVSp_Short))
  
  # Extract the amino acid position for matching entries
  for (i in seq_along(matches)) {
    if (length(matches[[i]]) > 1) {
      gene_maf$aa_pos[i] <- as.numeric(matches[[i]][2])
    }
  }
  
  # Handle cases where regex didn't work (like frameshifts or splicing variants)
  frameshifts <- grep("^p\\.[A-Za-z]+[0-9]+fs", gene_maf$HGVSp_Short)
  if (length(frameshifts) > 0) {
    # Extract positions from frameshifts (e.g., "p.G542fs" -> 542)
    fs_pattern <- "p\\.[A-Za-z]+([0-9]+)fs"
    fs_matches <- regmatches(gene_maf$HGVSp_Short[frameshifts], 
                             regexec(fs_pattern, gene_maf$HGVSp_Short[frameshifts]))
    
    for (i in seq_along(fs_matches)) {
      if (length(fs_matches[[i]]) > 1) {
        gene_maf$aa_pos[frameshifts[i]] <- as.numeric(fs_matches[[i]][2])
      }
    }
  }
  
  # For splice sites or other variants without amino acid positions, 
  # set position to -1 (will be handled specially in the plot)
  gene_maf$aa_pos[is.na(gene_maf$aa_pos)] <- -1
  
  # Add a shortened version of HGVSp_Short for labeling (without the 'p.' prefix)
  gene_maf$hgvsp_label <- gsub("^p\\.", "", gene_maf$HGVSp_Short)
  
  # For any variant type not in our predefined colors, add a default color
  unique_variants <- unique(gene_maf$Variant_Classification)
  missing_variants <- setdiff(unique_variants, names(variant_colors))
  if (length(missing_variants) > 0) {
    # Generate some additional colors for any missing variant types
    additional_colors <- rainbow(length(missing_variants))
    variant_colors <- c(variant_colors, setNames(additional_colors, missing_variants))
  }
  
  # Calculate number of genes to process
  gene_count <- length(unique(gene_maf$Hugo_Symbol))
  
  # Calculate grid layout
  n_rows <- ceiling(gene_count / max_genes_per_row)
  n_cols <- min(gene_count, max_genes_per_row)
  
  # Create a plot for each gene
  plot_list <- list()
  unique_genes <- unique(gene_maf$Hugo_Symbol)
  
  for (i in seq_along(unique_genes)) {
    gene <- unique_genes[i]
    message(paste0("Processing gene: ", gene))
    
    # Filter data for this gene
    gene_data <- gene_maf[gene_maf$Hugo_Symbol == gene, ]
    
    if (nrow(gene_data) == 0) {
      message(paste0("No mutations found for gene: ", gene))
      next
    }
    
    # Order samples for vertical axis
    gene_data$Tumor_Sample_Barcode <- factor(gene_data$Tumor_Sample_Barcode, 
                                           levels = rev(unique(gene_data$Tumor_Sample_Barcode)))
    
    # Get variant types (SNV vs INDEL)
    # Assume anything not a missense/nonsense/silent is an INDEL
    snv_types <- c("Missense_Mutation", "Nonsense_Mutation", "Silent")
    gene_data$variant_type <- ifelse(gene_data$Variant_Classification %in% snv_types, "SNV", "INDEL")
    
    # Create the plot
    p <- ggplot2::ggplot(gene_data, 
                        ggplot2::aes(x = aa_pos, 
                                     y = Tumor_Sample_Barcode, 
                                     color = Variant_Classification,
                                     shape = variant_type)) +
      ggplot2::geom_point(size = point_size, alpha = 0.8) +
      ggplot2::scale_color_manual(values = variant_colors) +
      ggplot2::scale_shape_manual(values = c("SNV" = 16, "INDEL" = 17)) + # Circle for SNV, triangle for INDEL
      ggrepel::geom_text_repel(
        ggplot2::aes(label = hgvsp_label),
        size = text_size,
        color = "black",
        box.padding = 0.35,
        point.padding = 0.5,
        segment.color = "grey50",
        max.overlaps = 15
      ) +
      ggplot2::labs(
        title = paste0(gene),  # Simplified title
        x = "Amino Acid Position", 
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        legend.position = "right",
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
    
    # If we have special cases with no amino acid position, add a special section
    if (any(gene_data$aa_pos == -1)) {
      # Get the maximum AA position to determine where to place special cases
      max_aa <- max(gene_data$aa_pos[gene_data$aa_pos > 0], na.rm = TRUE)
      
      # Update positions for special cases to be at the rightmost end of the plot
      gene_data$aa_pos[gene_data$aa_pos == -1] <- max_aa * 1.1
      
      # Add an annotation for these special cases
      p <- p + ggplot2::annotate("text", 
                               x = max_aa * 1.1, 
                               y = 0, 
                               label = "Special\nvariants", 
                               size = 3,
                               hjust = 0.5)
    }
    
    plot_list[[i]] <- p
  }
  
  # Combine plots in a grid
  final_plot <- do.call(gridExtra::grid.arrange, 
                       c(plot_list, 
                         ncol = n_cols, 
                         nrow = n_rows))
  
  # Save the plot if an output directory is specified
  if (!is.null(output_dir)) {
    pdf_file <- file.path(output_dir, "maf_amino_acid_grid.pdf")
    message(paste0("Saving plot to: ", pdf_file))
    
    # Calculate appropriate dimensions based on grid size
    width_adjusted <- width * min(n_cols, max_genes_per_row) / max_genes_per_row
    height_adjusted <- height * n_rows / 2
    
    ggplot2::ggsave(
      pdf_file,
      plot = final_plot,
      width = width_adjusted,
      height = height_adjusted
    )
    
    # Also save individual gene plots if there are many
    if (gene_count > 4) {
      message("Saving individual gene plots as well...")
      for (i in seq_along(plot_list)) {
        gene <- unique_genes[i]
        gene_file <- file.path(output_dir, paste0(gene, "_aa_mutations.pdf"))
        ggplot2::ggsave(
          gene_file,
          plot = plot_list[[i]],
          width = width / 2,
          height = height / 2
        )
      }
    }
  }
  
  # Return the plot objects invisibly
  invisible(list(
    combined_plot = final_plot,
    individual_plots = plot_list
  ))
}

#' Alternative implementation using facets instead of grid for genes
#'
#' @param maf_input Either a MAF object from maftools or a data frame with MAF format data
#' @param genes Character vector of gene symbols to visualize
#' @param output_dir Directory to save plots (NULL = don't save, just display)
#' @param height Plot height in inches when saving to PDF (default: 10)
#' @param width Plot width in inches when saving to PDF (default: 12)
#' @param text_size Size of the HGVSp label text (default: 2.5)
#' @param point_size Size of mutation points (default: 3)
#'
#' @return Invisibly returns the ggplot object
#' @export
#'
visualize_maf_aa_facet <- function(maf_input, 
                                  genes, 
                                  title="",
                                  output_dir = NULL,
                                  height = 10,
                                  width = 12,
                                  text_size = 2.5,
                                  point_size = 3) {
  
  # Check required packages
  required_packages <- c("ggplot2", "data.table", "dplyr", "stringr", "maftools", "ggrepel")
  for(pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed. ",
                  "Please install it with: install.packages('", pkg, "')"))
    }
  }
  
  # Access the MAF data
  message("Reading MAF data...")
  if (inherits(maf_input, "MAF")) {
    # If input is a MAF object from maftools
    maf_data <- maf_input@data
  } else if (is.data.frame(maf_input)) {
    # If input is already a data frame
    maf_data <- maf_input
  } else {
    stop("Input must be either a MAF object from maftools or a data frame in MAF format")
  }
  
  # Check if required columns exist
  required_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "HGVSp_Short")
  missing_cols <- required_cols[!required_cols %in% colnames(maf_data)]
  if (length(missing_cols) > 0) {
    stop(paste0("MAF data is missing required columns: ", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Setup custom variant classification colors
  variant_colors <- c(
    "Missense_Mutation" = "#FF5733",     # Orange-red
    "Nonsense_Mutation" = "#3498DB",     # Blue
    "Silent" = "#95A5A6",                # Gray
    "Frame_Shift_Del" = "#9B59B6",       # Purple
    "Frame_Shift_Ins" = "#F1C40F",       # Yellow
    "In_Frame_Del" = "#E74C3C",          # Red
    "In_Frame_Ins" = "#2ECC71",          # Green
    "Splice_Site" = "#34495E",           # Dark blue
    "Translation_Start_Site" = "#D4AC0D", # Gold
    "Nonstop_Mutation" = "#A04000"       # Brown
  )
  
  # Create output directory if needed
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # Filter MAF for the specified genes
  gene_maf <- maf_data[maf_data$Hugo_Symbol %in% genes, ]
  
  if (nrow(gene_maf) == 0) {
    stop("No mutations found for any of the specified genes")
  }
  
  # Extract amino acid positions from HGVSp_Short
  gene_maf$aa_pos <- NA
  
  # Extract positions using regex
  pattern <- "p\\.[A-Za-z]+([0-9]+)[A-Za-z\\*]+"
  matches <- regmatches(gene_maf$HGVSp_Short, 
                        regexec(pattern, gene_maf$HGVSp_Short))
  
  # Extract the amino acid position for matching entries
  for (i in seq_along(matches)) {
    if (length(matches[[i]]) > 1) {
      gene_maf$aa_pos[i] <- as.numeric(matches[[i]][2])
    }
  }
  
  # Handle cases where regex didn't work (like frameshifts or splicing variants)
  frameshifts <- grep("^p\\.[A-Za-z]+[0-9]+fs", gene_maf$HGVSp_Short)
  if (length(frameshifts) > 0) {
    # Extract positions from frameshifts (e.g., "p.G542fs" -> 542)
    fs_pattern <- "p\\.[A-Za-z]+([0-9]+)fs"
    fs_matches <- regmatches(gene_maf$HGVSp_Short[frameshifts], 
                             regexec(fs_pattern, gene_maf$HGVSp_Short[frameshifts]))
    
    for (i in seq_along(fs_matches)) {
      if (length(fs_matches[[i]]) > 1) {
        gene_maf$aa_pos[frameshifts[i]] <- as.numeric(fs_matches[[i]][2])
      }
    }
  }
  
  # For splice sites or other variants without amino acid positions, 
  # set position to -1 (will be handled specially in the plot)
  gene_maf$aa_pos[is.na(gene_maf$aa_pos)] <- -1
  
  # Add a shortened version of HGVSp_Short for labeling (without the 'p.' prefix)
  gene_maf$hgvsp_label <- gsub("^p\\.", "", gene_maf$HGVSp_Short)
  
  # For any variant type not in our predefined colors, add a default color
  unique_variants <- unique(gene_maf$Variant_Classification)
  missing_variants <- setdiff(unique_variants, names(variant_colors))
  if (length(missing_variants) > 0) {
    # Generate some additional colors for any missing variant types
    additional_colors <- rainbow(length(missing_variants))
    variant_colors <- c(variant_colors, setNames(additional_colors, missing_variants))
  }
  
  # Get variant types (SNV vs INDEL)
  # Assume anything not a missense/nonsense/silent is an INDEL
  snv_types <- c("Missense_Mutation", "Nonsense_Mutation", "Silent")
  gene_maf$variant_type <- ifelse(gene_maf$Variant_Classification %in% snv_types, "SNV", "INDEL")
  
  # Order samples for vertical axis
  # First, count mutations per sample to rank them
  sample_counts <- table(gene_maf$Tumor_Sample_Barcode)
  sample_order <- names(sort(sample_counts, decreasing = TRUE))
  gene_maf$Tumor_Sample_Barcode <- factor(gene_maf$Tumor_Sample_Barcode, 
                                          levels = rev(sample_order))
  
  # Make Hugo_Symbol a factor with specified order
  gene_maf$Hugo_Symbol <- factor(gene_maf$Hugo_Symbol, levels = genes)
  
  # Create the faceted plot
  p <- ggplot2::ggplot(gene_maf, 
                      ggplot2::aes(x = aa_pos, 
                                   y = Tumor_Sample_Barcode, 
                                   color = Variant_Classification,
                                   shape = variant_type)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_color_manual(values = variant_colors) +
    ggplot2::scale_shape_manual(values = c("SNV" = 16, "INDEL" = 17)) + # Circle for SNV, triangle for INDEL
    ggrepel::geom_text_repel(
      ggplot2::aes(label = hgvsp_label),
      size = text_size,
      color = "black",
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = 15
    ) +
    ggplot2::facet_wrap(~ Hugo_Symbol, scales = "free_x") +
    ggplot2::labs(
      title = title,
      x = "Amino Acid Position", 
      y = "Sample"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      axis.text.y = ggplot2::element_text(size = 8),
      strip.background = ggplot2::element_rect(fill = "lightblue"),
      strip.text = ggplot2::element_text(face = "bold")
    )
  
  # Save the plot if an output directory is specified
  if (!is.null(output_dir)) {
    pdf_file <- file.path(output_dir, "maf_amino_acid_facet.pdf")
    message(paste0("Saving plot to: ", pdf_file))
    
    ggplot2::ggsave(
      pdf_file,
      plot = p,
      width = width,
      height = height
    )
  }
  
  # Return the plot object invisibly
  invisible(p)
}

# Example usage:
# library(maftools)
# 
# # Using a MAF object:
# maf <- read.maf("example.maf")
# genes_of_interest <- c("TP53", "PIK3CA", "ARID1A")
# visualize_maf_aa_grid(maf, genes_of_interest, output_dir = "plots")
# 
# # Or using a data frame directly:
# maf_df <- maf@data
# visualize_maf_aa_facet(maf_df, genes_of_interest, output_dir = "plots")
