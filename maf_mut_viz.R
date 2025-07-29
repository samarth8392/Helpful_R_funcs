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
#' @param sample_order Optional character vector specifying the order of tumor samples from top to bottom
#' @param title Plot title (default: "")
#' @param output_dir Directory to save plots (NULL = don't save, just display)
#' @param file_type Type of output file ("pdf" or "png", default: "pdf")
#' @param height Plot height in inches when saving (default: 10)
#' @param width Plot width in inches when saving (default: 12)
#' @param dpi Resolution for PNG output in dots per inch (default: 300)
#' @param text_size Size of the HGVSp label text (default: 2.5)
#' @param point_size Size of mutation points (default: 3)
#'
#' @return Invisibly returns the ggplot object
#' @export
#'
visualize_maf_aa_facet <- function(maf_input, 
                                  genes, 
                                  sample_order = NULL,
                                  patient_order = NULL,  # New parameter for patient ordering
                                  title = "",
                                  file_name = "maf_muts",
                                  output_dir = NULL,
                                  file_type = "pdf",
                                  height = 10,
                                  width = 12,
                                  dpi = 300,
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
  required_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "HGVSp_Short", "PatientID")
  missing_cols <- required_cols[!required_cols %in% colnames(maf_data)]
  if (length(missing_cols) > 0) {
    stop(paste0("MAF data is missing required columns: ", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Check if CLIN_SIG column exists, if not create it with default values
  if (!"CLIN_SIG" %in% colnames(maf_data)) {
    warning("CLIN_SIG column not found. Creating default column with all variants as 'unknown'.")
    maf_data$CLIN_SIG <- "unknown"
  }
  
  # Check if PolyPhen column exists, if not create it with default values
  if (!"PolyPhen" %in% colnames(maf_data)) {
    warning("PolyPhen column not found. Creating default column with all variants as NA.")
    maf_data$PolyPhen <- NA
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
  
  # Create pathogenic indicator (case-insensitive, check if "pathogenic" appears anywhere in CLIN_SIG)
  gene_maf$is_pathogenic <- grepl("pathogenic", tolower(gene_maf$CLIN_SIG))
  
  # Extract numeric PolyPhen scores from the formatted strings
  gene_maf$polyphen_score <- NA
  
  # Extract numeric values from PolyPhen strings like "benign(0.003)" or "probably_damaging(0.997)"
  polyphen_pattern <- "\\(([0-9.]+)\\)"
  matches <- regmatches(gene_maf$PolyPhen, regexec(polyphen_pattern, gene_maf$PolyPhen))
  
  for (i in seq_along(matches)) {
    if (length(matches[[i]]) > 1) {
      gene_maf$polyphen_score[i] <- as.numeric(matches[[i]][2])
    }
  }
  
  # Create PolyPhen border color based on score
  gene_maf$polyphen_border <- ifelse(is.na(gene_maf$polyphen_score), "black",
                                    ifelse(gene_maf$polyphen_score > 0.85, "red4", "olivedrab"))
  
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
  
  # Order samples and patients for vertical axis
  unique_patients <- unique(gene_maf$PatientID)
  unique_samples <- unique(gene_maf$Tumor_Sample_Barcode)
  
  # Create patient-sample mapping
  patient_sample_map <- gene_maf %>%
    dplyr::select(PatientID, Tumor_Sample_Barcode) %>%
    dplyr::distinct()
  
  # Handle patient ordering
  if (!is.null(patient_order)) {
    # Validate the patient order
    invalid_patients <- setdiff(patient_order, unique_patients)
    if (length(invalid_patients) > 0) {
      warning(paste0("The following patients in patient_order were not found in the data: ", 
                     paste(invalid_patients, collapse = ", ")))
    }
    
    missing_patients <- setdiff(unique_patients, patient_order)
    if (length(missing_patients) > 0) {
      warning(paste0("The following patients in the data were not specified in patient_order and will be placed at the bottom: ", 
                     paste(missing_patients, collapse = ", ")))
      # Add unspecified patients at the end of the order
      patient_order <- c(patient_order, missing_patients)
    }
  } else {
    # Default ordering by mutation count per patient if no patient_order provided
    message("No custom patient order provided. Ordering patients by mutation count.")
    patient_counts <- table(gene_maf$PatientID)
    patient_order <- names(sort(patient_counts, decreasing = TRUE))
  }
  
  # Handle sample ordering within each patient
  if (!is.null(sample_order)) {
    # For samples that are specified in sample_order, use that order
    # For samples not specified, order them by mutation count within each patient
    ordered_samples <- c()
    
    for (patient in patient_order) {
      patient_samples <- patient_sample_map$Tumor_Sample_Barcode[patient_sample_map$PatientID == patient]
      
      # Get samples for this patient that are in sample_order
      ordered_patient_samples <- intersect(sample_order, patient_samples)
      
      # Get samples for this patient that are NOT in sample_order
      unordered_patient_samples <- setdiff(patient_samples, sample_order)
      
      # For unordered samples, sort by mutation count
      if (length(unordered_patient_samples) > 0) {
        sample_counts <- table(gene_maf$Tumor_Sample_Barcode[gene_maf$PatientID == patient & 
                                                           gene_maf$Tumor_Sample_Barcode %in% unordered_patient_samples])
        unordered_patient_samples <- names(sort(sample_counts, decreasing = TRUE))
      }
      
      # Combine ordered and unordered samples for this patient
      ordered_samples <- c(ordered_samples, ordered_patient_samples, unordered_patient_samples)
    }
  } else {
    # Default ordering by mutation count within each patient
    message("No custom sample order provided. Ordering samples by mutation count within each patient.")
    ordered_samples <- c()
    
    for (patient in patient_order) {
      patient_samples <- patient_sample_map$Tumor_Sample_Barcode[patient_sample_map$PatientID == patient]
      patient_maf <- gene_maf[gene_maf$PatientID == patient, ]
      
      if (nrow(patient_maf) > 0) {
        sample_counts <- table(patient_maf$Tumor_Sample_Barcode)
        ordered_patient_samples <- names(sort(sample_counts, decreasing = TRUE))
        ordered_samples <- c(ordered_samples, ordered_patient_samples)
      }
    }
  }
  
  # Set factor levels for samples (reverse for plotting from top to bottom)
  gene_maf$Tumor_Sample_Barcode <- factor(gene_maf$Tumor_Sample_Barcode, 
                                         levels = rev(ordered_samples))
  
  # Set factor levels for patients (reverse for plotting from top to bottom)
  gene_maf$PatientID <- factor(gene_maf$PatientID, levels = rev(patient_order))
  
  # Make Hugo_Symbol a factor with specified order
  gene_maf$Hugo_Symbol <- factor(gene_maf$Hugo_Symbol, levels = genes)
  
  # Create the faceted plot with patient grouping
  p <- ggplot2::ggplot(gene_maf, 
                      ggplot2::aes(x = aa_pos, 
                                   y = Tumor_Sample_Barcode, 
                                   color = Variant_Classification,
                                   shape = variant_type)) +
    # Add regular points
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    # Add PolyPhen-based borders - separate for SNVs (circles) and INDELs (triangles)
    ggplot2::geom_point(data = gene_maf[gene_maf$variant_type == "SNV", ],
                       ggplot2::aes(x = aa_pos, y = Tumor_Sample_Barcode),
                       size = point_size,
                       shape = 21,  # Hollow circle
                       color = gene_maf$polyphen_border[gene_maf$variant_type == "SNV"],
                       fill = NA,
                       stroke = 1,
                       alpha = 1,
                       show.legend = FALSE,
                       inherit.aes = FALSE) +
    ggplot2::geom_point(data = gene_maf[gene_maf$variant_type == "INDEL", ],
                       ggplot2::aes(x = aa_pos, y = Tumor_Sample_Barcode),
                       size = point_size,
                       shape = 24,  # Hollow triangle
                       color = gene_maf$polyphen_border[gene_maf$variant_type == "INDEL"],
                       fill = NA,
                       stroke = 1,
                       alpha = 1,
                       show.legend = FALSE,
                       inherit.aes = FALSE) +
    # Add pathogenic variants with thick black border - circles for SNV, triangles for INDEL
    ggplot2::geom_point(data = gene_maf[gene_maf$is_pathogenic & gene_maf$variant_type == "SNV", ],
                       ggplot2::aes(x = aa_pos, y = Tumor_Sample_Barcode),
                       size = point_size,
                       stroke = 3,
                       color = "black",
                       fill = NA,
                       alpha = 1,
                       shape = 21,  # Hollow circle for SNV
                       inherit.aes = FALSE) +
    ggplot2::geom_point(data = gene_maf[gene_maf$is_pathogenic & gene_maf$variant_type == "INDEL", ],
                       ggplot2::aes(x = aa_pos, y = Tumor_Sample_Barcode),
                       size = point_size,
                       stroke = 3,
                       color = "black",
                       fill = NA,
                       alpha = 1,
                       shape = 24,  # Hollow triangle for INDEL
                       inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = variant_colors) +
    ggplot2::scale_shape_manual(values = c("SNV" = 16, "INDEL" = 17)) + # Circle for SNV, triangle for INDEL
    ggrepel::geom_label_repel(
      data = gene_maf[gene_maf$variant_type == "SNV", ],
      ggplot2::aes(label = hgvsp_label),
      size = text_size,
      color = "black",
      fill = scales::alpha(gene_maf$polyphen_border, 0.5),
      label.color = NA,  # No outline/border on the rectangle
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = 15
    ) +
    ggplot2::facet_grid(PatientID ~ Hugo_Symbol, scales = "free", space = "free_y") +
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
      strip.text = ggplot2::element_text(face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0)  # Keep patient labels horizontal
    )
  
  # Add manual legends using annotation instead of dummy geoms
  legend_guides <- list()
  
  # Add pathogenic legend if there are pathogenic variants
  if (any(gene_maf$is_pathogenic)) {
    legend_guides[["Clinical Significance"]] <- ggplot2::guide_legend(
      title = "Clinical Significance",
      override.aes = list(
        size = point_size + 1,
        stroke = 3,
        color = "black",
        fill = NA,
        shape = 21,
        alpha = 1
      ),
      order = 3
    )
    
    # Add a dummy aesthetic for pathogenic variants
    p <- p + 
      ggplot2::aes(linetype = ifelse(is_pathogenic, "Pathogenic", "Non-pathogenic")) +
      ggplot2::scale_linetype_manual(
        values = c("Pathogenic" = "solid", "Non-pathogenic" = "solid"),
        name = "Clinical Significance",
        guide = legend_guides[["Clinical Significance"]]
      )
  }
  
  # Create a simpler approach for the legends by using the guides() function
  p <- p + ggplot2::guides(
    color = ggplot2::guide_legend(title = "Variant Classification", order = 1),
    shape = ggplot2::guide_legend(title = "Variant Type", order = 2)
  )
  
  # Add text annotation for PolyPhen legend if there are PolyPhen scores
  if (any(!is.na(gene_maf$polyphen_score))) {
    # Add a text annotation explaining the border colors
    p <- p + 
      ggplot2::labs(caption = "Border colors: Red = PolyPhen > 0.85 (damaging), Green = PolyPhen ≤ 0.85 (benign), Black = No PolyPhen score\nThick black border = Pathogenic variant")
  } else if (any(gene_maf$is_pathogenic)) {
    p <- p + 
      ggplot2::labs(caption = "Thick black border = Pathogenic variant")
  }
  
  # Save the plot if an output directory is specified
  if (!is.null(output_dir)) {
    # Validate file type
    file_type <- tolower(file_type)
    if (!file_type %in% c("pdf", "png")) {
      warning("Invalid file_type. Must be 'pdf' or 'png'. Defaulting to 'pdf'.")
      file_type <- "pdf"
    }
    
    # Set file extension based on file type
    file_extension <- ifelse(file_type == "pdf", "pdf", "png")
    file_name <- file.path(output_dir, paste0(file_name,".",file_extension))
    message(paste0("Saving plot to: ", file_name))
    
    # Save as specified file type
    if (file_type == "pdf") {
      ggplot2::ggsave(
        file_name,
        plot = p,
        width = width,
        height = height
      )
    } else {
      ggplot2::ggsave(
        file_name,
        plot = p,
        width = width,
        height = height,
        dpi = dpi
      )
    }
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
# 
# # With default sample ordering (by mutation count) and PDF output:
# visualize_maf_aa_facet(maf, genes_of_interest, output_dir = "plots")
# 
# # With custom sample ordering:
# custom_order <- c("TCGA-A1-A0SD", "TCGA-A2-A0T2", "TCGA-AR-A1AR")
# visualize_maf_aa_facet(maf, genes_of_interest, sample_order = custom_order, output_dir = "plots")
# 
# # Saving as PNG with 300 DPI:
# visualize_maf_aa_facet(maf, genes_of_interest, output_dir = "plots", file_type = "png", dpi = 300)
# 
# # Or using a data frame directly:
# maf_df <- maf@data
# visualize_maf_aa_facet(maf_df, genes_of_interest, sample_order = custom_order, output_dir = "plots")
