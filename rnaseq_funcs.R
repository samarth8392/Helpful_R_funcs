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

create_enhanced_volcano_deseq2 <- function(deseq_results, 
                                          highlight_genes = NULL,
                                          comparison_name = "DESeq2 Analysis",
                                          padj_threshold = 0.05,
                                          lfc_threshold = 1.0,
                                          label_top_genes = 10,
                                          point_size = 1.5,
                                          highlight_size = 3) {
  
  # Load required libraries
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  require(scales)
  
  # Input validation
  required_cols <- c("log2FoldChange", "padj", "baseMean")
  if (!all(required_cols %in% colnames(deseq_results))) {
    stop("deseq_results must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Prepare data
  plot_data <- deseq_results %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      # Create significance categories
      significance = dplyr::case_when(
        padj < padj_threshold & abs(log2FoldChange) > lfc_threshold ~ "Significant",
        padj < padj_threshold & abs(log2FoldChange) <= lfc_threshold ~ "P-value only",
        padj >= padj_threshold & abs(log2FoldChange) > lfc_threshold ~ "Effect size only",
        TRUE ~ "Not significant"
      ),
      # Direction of regulation
      direction = dplyr::case_when(
        log2FoldChange > 0 ~ "Upregulated",
        log2FoldChange < 0 ~ "Downregulated",
        TRUE ~ "No change"
      ),
      # Transform p-values
      neg_log10_padj = -log10(padj),
      # Check if gene should be highlighted
      is_highlight = if (!is.null(highlight_genes)) {
        gene_id %in% highlight_genes | 
        (exists("hugo_symbol", where = .) && hugo_symbol %in% highlight_genes)
      } else {
        FALSE
      },
      # Create labels for top genes
      label = case_when(
        is_highlight ~ if(exists("hugo_symbol", where = .)) hugo_symbol else gene_id,
        TRUE ~ ""
      )
    )
  
  # If no highlight genes provided, label top genes by significance
  if (is.null(highlight_genes) && label_top_genes > 0) {
    top_genes <- plot_data %>%
      dplyr::filter(significance == "Significant") %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = label_top_genes)
    
    plot_data <- plot_data %>%
      dplyr::mutate(
        label = case_when(
          gene_id %in% top_genes$gene_id ~ if(exists("hugo_symbol", where = .)) hugo_symbol else gene_id,
          TRUE ~ ""
        ),
        is_highlight = gene_id %in% top_genes$gene_id
      )
  }
  
  # Color scheme
  colors <- c(
    "Significant" = "#E31A1C",        # Red
    "P-value only" = "#FF7F00",       # Orange  
    "Effect size only" = "#1F78B4",   # Blue
    "Not significant" = "#CCCCCC"     # Gray
  )
  
  # Count significant genes
  sig_counts <- plot_data %>%
    dplyr::filter(significance == "Significant") %>%
    dplyr::summarise(
      total = dplyr::n(),
      up = sum(log2FoldChange > 0),
      down = sum(log2FoldChange < 0),
      .groups = 'drop'
    )
  
  # Create the volcano plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)) +
    # Main points
    ggplot2::geom_point(ggplot2::aes(color = significance, 
                                    size = is_highlight, 
                                    alpha = ifelse(is_highlight, 1, 0.7))) +
    
    # Highlighted points with special styling
    {if (any(plot_data$is_highlight)) {
      ggplot2::geom_point(data = plot_data %>% dplyr::filter(is_highlight), 
                         ggplot2::aes(x = log2FoldChange, y = neg_log10_padj),
                         color = "black", size = highlight_size + 1, shape = 21, 
                         fill = "yellow", stroke = 2, alpha = 1)
    }} +
    
    # Threshold lines
    ggplot2::geom_hline(yintercept = -log10(padj_threshold), 
                       linetype = "dashed", color = "gray50", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
                       linetype = "dashed", color = "gray50", alpha = 0.8) +
    
    # Gene labels
    {if (any(plot_data$label != "")) {
      ggrepel::geom_text_repel(data = plot_data %>% dplyr::filter(label != ""),
                              ggplot2::aes(label = label),
                              box.padding = 0.5,
                              point.padding = 0.5,
                              segment.color = "black",
                              segment.size = 0.5,
                              min.segment.length = 0,
                              max.overlaps = Inf,
                              size = 3.5,
                              fontface = "bold")
    }} +
    
    # Scales
    ggplot2::scale_color_manual(values = colors, name = "Significance") +
    ggplot2::scale_size_manual(values = c("FALSE" = point_size, "TRUE" = highlight_size), 
                              guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    
    # Labels
    ggplot2::labs(
      title = paste("Differential Gene Expression:", comparison_name),
      subtitle = paste0("Significant DEGs: ", sig_counts$total, 
                       " (", sig_counts$up, " up, ", sig_counts$down, " down) | ",
                       "Total genes: ", nrow(plot_data)),
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10]("Adjusted P-value")),
      caption = bquote("Significance thresholds: |log"[2]*"FC| > "*.(lfc_threshold)*
                      ", Adjusted P < "*.(padj_threshold))
    ) +
    
    # Theme
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, color = "gray30"),
      plot.caption = ggplot2::element_text(hjust = 1, size = 9, color = "gray50"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, size = 0.5)
    )
  
  return(p)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Basic usage with your DESeq2 results
volcano_plot <- create_enhanced_volcano_deseq2(
  deseq_results = res_gly1_HFCFD_shrunk,  # Your shrunk results
  comparison_name = "Glyphosate 1.0 vs Control (HFCFD)",
  padj_threshold = 0.05,
  lfc_threshold = 1.0,
  label_top_genes = 15
)

print(volcano_plot)

# With specific genes to highlight
genes_of_interest <- c("Adam11", "Dse", "Cd5l", "Orm3", "Ear2")

volcano_plot_highlighted <- create_enhanced_volcano_deseq2(
  deseq_results = res_gly1_HFCFD_shrunk,
  highlight_genes = genes_of_interest,
  comparison_name = "Glyphosate 1.0 vs Control (HFCFD)",
  padj_threshold = 0.05,
  lfc_threshold = 0.5,  # Lower threshold for more sensitive detection
  label_top_genes = 0   # Don't auto-label, only highlight specified genes
)

print(volcano_plot_highlighted)

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
