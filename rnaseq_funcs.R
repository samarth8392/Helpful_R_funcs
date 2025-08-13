# ============================================================================
# Functions for analyzing differential expression
# ============================================================================



# Load required libraries
suppressPackageStartupMessages({
  library(DT)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
})

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Create Interactive DataTable
create_dt <- function(x, caption = "") {
  DT::datatable(x,
                extensions = 'Buttons',
                rownames = FALSE,
                options = list(
                  dom = 'Blfrtip',
                  buttons = c('copy', 'csv', 'excel'),
                  pageLength = 20,
                  lengthMenu = list(c(20, 35, 50, -1),
                                   c(20, 35, 50, "All")),
                  scrollX = TRUE, 
                  scrollY = 500
                ),
                caption = caption)
}

#' Remove Outliers Based on Percentiles
remove_outliers <- function(x) {
  lower_bound <- quantile(x, 0.05, na.rm = TRUE)
  upper_bound <- quantile(x, 0.95, na.rm = TRUE)
  return(dplyr::between(x, lower_bound, upper_bound))
}

#' Count Zero Values
count_zeros <- function(x) {
  sum(x == 0, na.rm = TRUE)
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

# Enhanced Volcano Plot Function for DESeq2 Results
# Adapted from methylation analysis volcano plot

# Enhanced create_enhanced_volcano_deseq2 function with EnhancedVolcano options
# Enhanced create_enhanced_volcano_deseq2 function with EnhancedVolcano options
create_enhanced_volcano_deseq2 <- function(deseq_results, 
                                           highlight_genes = NULL,
                                           comparison_name = "Comparison",
                                           # EnhancedVolcano-style parameters
                                           lab = NULL,
                                           x = 'log2FoldChange',
                                           y = 'padj',
                                           title = NULL,
                                           subtitle = NULL,
                                           caption = NULL,
                                           FCcutoff = 0.5,
                                           pCutoff = 0.05,
                                           xlim = NULL,
                                           ylim = NULL,
                                           # Color options
                                           col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                                           colAlpha = 0.7,
                                           # Label options
                                           label_top_genes = 10,
                                           labSize = 3.0,
                                           labCol = 'black',
                                           labFace = 'plain',
                                           boxedLabels = FALSE,
                                           # Point options
                                           pointSize = 1.2,
                                           # Grid options
                                           gridlines.major = TRUE,
                                           gridlines.minor = FALSE,
                                           # Legend options
                                           legendPosition = 'bottom',
                                           legendLabSize = 11,
                                           legendIconSize = 3.0,
                                           # Threshold line options
                                           hline = NULL,
                                           vline = NULL,
                                           drawConnectors = FALSE,
                                           save_plot = FALSE,
                                           filename = NULL) {
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # Convert to data frame and handle gene_id column properly
  if (is.data.frame(deseq_results)) {
    df <- deseq_results
    if (!"gene_id" %in% colnames(df)) {
      df <- df %>% tibble::rownames_to_column("gene_id")
    }
  } else {
    df <- as.data.frame(deseq_results) %>%
      tibble::rownames_to_column("gene_id")
  }
  
  # Handle lab parameter
  if (is.null(lab)) {
    if ("GeneName" %in% colnames(df)) {
      df$lab <- ifelse(is.na(df$GeneName) | df$GeneName == "", df$gene_id, df$GeneName)
    } else {
      df$lab <- df$gene_id
    }
  } else {
    if (is.character(lab) && length(lab) == 1 && lab %in% colnames(df)) {
      df$lab <- df[[lab]]
    } else if (length(lab) == nrow(df)) {
      df$lab <- lab
    } else {
      stop("lab parameter must be a column name in the data or a vector of same length as data")
    }
  }
  
  # Ensure required columns exist
  if (!x %in% colnames(df)) stop(paste("Column", x, "not found in data"))
  if (!y %in% colnames(df)) stop(paste("Column", y, "not found in data"))
  
  # Remove rows with NA values in essential columns
  df <- df %>%
    filter(!is.na(!!sym(x)) & !is.na(!!sym(y)))
  
  # Handle y-axis transformation (padj vs pvalue)
  if (y == 'padj' || y == 'pvalue') {
    df$y_transformed <- -log10(df[[y]])
    y_label <- paste0("-log10(", y, ")")
    if (is.null(hline)) hline <- -log10(pCutoff)
  } else {
    df$y_transformed <- df[[y]]
    y_label <- y
  }
  
  # Create significance categories using FCcutoff and pCutoff
  df <- df %>%
    mutate(
      significance = case_when(
        !!sym(y) < pCutoff & !!sym(x) >= FCcutoff ~ "Upregulated",
        !!sym(y) < pCutoff & !!sym(x) <= -FCcutoff ~ "Downregulated",
        !!sym(y) < pCutoff & abs(!!sym(x)) < FCcutoff ~ "Significant (small effect)",
        TRUE ~ "Not significant"
      )
    )
  
  # Set up colors (matching your original plot)
  colors <- c(
    "Not significant" = "grey70",        # Light grey for NS points
    "Significant (small effect)" = "#ff7f00",  # Orange for P-value only  
    "Upregulated" = "#e31a1c",          # Red for significant up
    "Downregulated" = "#e31a1c"         # Red for significant down (both directions same color)
  )
  
  # Get top genes for labeling (with yellow highlighting like your plot)
  if (label_top_genes > 0) {
    top_genes <- df %>%
      filter(significance %in% c("Upregulated", "Downregulated")) %>%
      arrange(!!sym(y)) %>%
      head(label_top_genes)
  } else {
    top_genes <- data.frame()
  }
  
  # Count significant genes
  sig_counts <- df %>% count(significance)
  up_count <- sig_counts$n[sig_counts$significance == "Upregulated"] %||% 0
  down_count <- sig_counts$n[sig_counts$significance == "Downregulated"] %||% 0
  
  # Set default title if not provided
  if (is.null(title)) title <- paste("Volcano Plot:", comparison_name)
  if (is.null(subtitle)) subtitle <- paste0("FC cutoff: ", FCcutoff, ", p-value cutoff: ", pCutoff)
  
  # Create the plot
  p <- ggplot(df, aes(x = !!sym(x), y = y_transformed)) +
    # Add all points first
    geom_point(aes(color = significance), alpha = colAlpha, size = pointSize) +
    
    # Add yellow circles around top genes (like your plot)
    {if (nrow(top_genes) > 0) {
      geom_point(data = top_genes, 
                 aes(x = !!sym(x), y = y_transformed),
                 color = "gold", fill = "yellow", 
                 shape = 21, size = pointSize + 2, stroke = 1.5, alpha = 0.8)
    }} +
    
    scale_color_manual(values = colors) +
    
    # Add threshold lines
    {if (!is.null(hline)) geom_hline(yintercept = hline, linetype = "dashed", color = "grey60", alpha = 0.8)} +
    {if (!is.null(vline)) geom_vline(xintercept = vline, linetype = "dashed", color = "grey60", alpha = 0.8)} +
    {if (is.null(vline)) geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed", color = "grey60", alpha = 0.8)} +
    
    # Add labels for top genes (black text like your plot)
    {if (nrow(top_genes) > 0) {
      geom_text_repel(
        data = top_genes,
        aes(label = lab),
        size = labSize,
        color = "black",
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 20,
        segment.color = 'black',
        segment.size = 0.3,
        segment.alpha = 0.6
      )
    }} +
    
    # Gene count annotation (like your plot style)
    annotate("text", x = Inf, y = Inf, 
             label = paste0("Significant DEGs: ", up_count + down_count, 
                           " (", up_count, " up, ", down_count, " down) | Total genes: ", nrow(df)),
             hjust = 1.02, vjust = 1.5, size = 3.5, color = "black") +
    
    # Set axis limits
    {if (!is.null(xlim)) scale_x_continuous(limits = xlim)} +
    {if (!is.null(ylim)) scale_y_continuous(limits = ylim)} +
    
    # Labels
    labs(
      title = title,
      subtitle = subtitle,
      caption = if(is.null(caption)) paste0("Significance thresholds: |log2FC| > ", FCcutoff, ", Adjusted P < ", pCutoff) else caption,
      x = paste0("log2(Fold Change)"),
      y = paste0("-log10(Adjusted P-value)"),
      color = "Significance"
    ) +
    
    # Theme (matching your plot)
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 1, size = 9, face = "italic", color = "grey50"),
      legend.position = "bottom",
      legend.title = element_text(size = legendLabSize, face = "bold"),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = legendIconSize, alpha = 1)))
  
  # Save plot if requested
  if (save_plot && !is.null(filename)) {
    ggsave(filename, p, width = 10, height = 8, dpi = 300)
    cat("Saved plot to:", filename, "\n")
  }
  
  return(p)
}

# Helper operator for default values
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Updated create_multiple_volcano_plots function
create_multiple_volcano_plots <- function(results_list, 
                                         comparison_names = NULL,
                                         highlight_genes = NULL,
                                         padj_threshold = 0.05,
                                         lfc_threshold = 1.0,
                                         save_plots = TRUE,
                                         plot_width = 10,
                                         plot_height = 8) {
  
  if (is.null(comparison_names)) {
    comparison_names <- names(results_list)
  }
  
  volcano_plots <- list()
  
  for (i in seq_along(results_list)) {
    cat("Creating volcano plot for:", comparison_names[i], "\n")
    
    volcano_plots[[i]] <- create_enhanced_volcano_deseq2(
      deseq_results = results_list[[i]],
      highlight_genes = highlight_genes,
      comparison_name = comparison_names[i],
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold,
      label_top_genes = 10
    )
    
    if (save_plots) {
      filename <- paste0("volcano_", gsub("[^A-Za-z0-9]", "_", comparison_names[i]), ".png")
      ggplot2::ggsave(filename, volcano_plots[[i]], 
                     width = plot_width, height = plot_height, dpi = 300)
      cat("Saved:", filename, "\n")
    }
  }
  
  names(volcano_plots) <- comparison_names
  return(volcano_plots)
}

# Updated create_multiple_volcano_plots function
create_multiple_volcano_plots <- function(results_list, 
                                         comparison_names = NULL,
                                         highlight_genes = NULL,
                                         padj_threshold = 0.05,
                                         lfc_threshold = 1.0,
                                         save_plots = TRUE,
                                         plot_width = 10,
                                         plot_height = 8) {
  
  if (is.null(comparison_names)) {
    comparison_names <- names(results_list)
  }
  
  volcano_plots <- list()
  
  for (i in seq_along(results_list)) {
    cat("Creating volcano plot for:", comparison_names[i], "\n")
    
    volcano_plots[[i]] <- create_enhanced_volcano_deseq2(
      deseq_results = results_list[[i]],
      highlight_genes = highlight_genes,
      comparison_name = comparison_names[i],
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold,
      label_top_genes = 10
    )
    
    if (save_plots) {
      filename <- paste0("volcano_", gsub("[^A-Za-z0-9]", "_", comparison_names[i]), ".png")
      ggplot2::ggsave(filename, volcano_plots[[i]], 
                     width = plot_width, height = plot_height, dpi = 300)
      cat("Saved:", filename, "\n")
    }
  }
  
  names(volcano_plots) <- comparison_names
  return(volcano_plots)
}

# =============================================================================
# BATCH VOLCANO PLOTS FOR MULTIPLE COMPARISONS
# =============================================================================

create_multiple_volcano_plots <- function(results_list, 
                                         comparison_names = NULL,
                                         highlight_genes = NULL,
                                         padj_threshold = 0.05,
                                         lfc_threshold = 1.0,
                                         save_plots = TRUE,
                                         plot_width = 10,
                                         plot_height = 8) {
  
  if (is.null(comparison_names)) {
    comparison_names <- names(results_list)
  }
  
  volcano_plots <- list()
  
  for (i in seq_along(results_list)) {
    cat("Creating volcano plot for:", comparison_names[i], "\n")
    
    volcano_plots[[i]] <- create_enhanced_volcano_deseq2(
      deseq_results = results_list[[i]],
      highlight_genes = highlight_genes,
      comparison_name = comparison_names[i],
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold,
      label_top_genes = 10
    )
    
    if (save_plots) {
      filename <- paste0("volcano_", gsub("[^A-Za-z0-9]", "_", comparison_names[i]), ".png")
      ggplot2::ggsave(filename, volcano_plots[[i]], 
                     width = plot_width, height = plot_height, dpi = 300)
      cat("Saved:", filename, "\n")
    }
  }
  
  names(volcano_plots) <- comparison_names
  return(volcano_plots)
}
