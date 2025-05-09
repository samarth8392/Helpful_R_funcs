# =====================================================================
# Enhanced Mutation Visualization Functions
# =====================================================================

#' Create detailed visualization for mutations in specific genes
#' 
#' @param maf MAF object from maftools
#' @param genes Character vector of gene names to visualize
#' @param plotType Type of plot ("lollipop", "custom", "oncostrip", "rainbowplot", "detailed")
#' @param showDomains Whether to show protein domains in lollipop plots
#' @param min_freq Minimum frequency of mutations to include in plot
#' @param customColors Color scheme for variant classifications
#' @param returnData Whether to return the plot data along with the plots
#' 
#' @return List of plots or plot data
#' 
visualize_gene_mutations <- function(maf, 
                                    genes, 
                                    plotType = "lollipop", 
                                    showDomains = TRUE,
                                    min_freq = 0,
                                    customColors = NULL,
                                    returnData = FALSE) {
  
  require(maftools)
  require(ggplot2)
  require(dplyr)
  
  result_list <- list()
  data_list <- list()
  
  # Process each gene
  for (gene in genes) {
    message(paste0("Processing gene: ", gene))
    
    # Check if gene exists in MAF
    if(!gene %in% unique(maf@data$Hugo_Symbol)) {
      message(paste0("Gene ", gene, " not found in the MAF object"))
      next
    }
    
    # Extract data for this gene
    gene_data <- maf@data[maf@data$Hugo_Symbol == gene,]
    
    if(nrow(gene_data) == 0) {
      message(paste0("No variants found for gene ", gene))
      next
    }
    
    # Store the data
    data_list[[gene]] <- gene_data
    
    # Process based on plot type
    if(plotType == "lollipop") {
      # Try standard lollipop plot
      tryCatch({
        plot <- lollipopPlot(maf = maf, 
                            gene = gene,
                            showMutationRate = TRUE,
                            labelPos = "all",
                            showDomainLabel = showDomains,
                            repel = TRUE,
                            printCount = TRUE)
        result_list[[gene]] <- plot
      }, error = function(e) {
        message(paste0("Error generating lollipop plot for ", gene, ": ", e$message))
        result_list[[gene]] <- NULL
      })
    } 
    else if(plotType == "oncostrip") {
      # Create oncostrip for just this gene
      tryCatch({
        plot <- oncostrip(maf = maf, 
                         genes = gene,
                         showTumorSampleBarcodes = TRUE,
                         annotationColor = customColors)
        result_list[[gene]] <- plot
      }, error = function(e) {
        message(paste0("Error generating oncostrip for ", gene, ": ", e$message))
        result_list[[gene]] <- NULL
      })
    }
    else if(plotType == "rainbowplot") {
      # Create rainfall plot which shows mutation distribution
      tryCatch({
        plot <- rainfallPlot(maf = maf, 
                            pointSize = 0.6,
                            gene = gene)
        result_list[[gene]] <- plot
      }, error = function(e) {
        message(paste0("Error generating rainfall plot for ", gene, ": ", e$message))
        result_list[[gene]] <- NULL
      })
    }
    else if(plotType == "custom") {
      # Create a custom visualization with flexible options
      # Extract mutation positions
      mutations <- try({
        data.frame(
          sample = gene_data$Tumor_Sample_Barcode,
          position = as.numeric(gsub(".*p\\.[A-Za-z]*(\\d+).*", "\\1", gene_data$HGVSp_Short)),
          type = gene_data$Variant_Classification,
          change = gene_data$HGVSp_Short,
          vaf = ifelse(!is.null(gene_data$t_vaf), gene_data$t_vaf, NA),
          impact = gene_data$IMPACT
        )
      }, silent = TRUE)
      
      if(inherits(mutations, "try-error") || all(is.na(mutations$position))) {
        message(paste0("Cannot extract mutation positions for gene ", gene))
        next
      }
      
      # Filter valid positions
      mutations <- mutations[!is.na(mutations$position),]
      
      if(nrow(mutations) > 0) {
        # Create basic plot
        p <- ggplot(mutations, aes(x = position)) +
          geom_histogram(aes(fill = type), binwidth = 5, color = "black", alpha = 0.8) +
          labs(title = paste0("Mutation distribution in ", gene),
               x = "Amino Acid Position",
               y = "Count") +
          theme_minimal() +
          theme(legend.title = element_text(size = 10),
                legend.text = element_text(size = 8),
                plot.title = element_text(hjust = 0.5, face = "bold"))
        
        result_list[[gene]] <- p
      } else {
        message("Cannot create plot: Unable to extract amino acid positions")
        result_list[[gene]] <- NULL
      }
    }
    else if(plotType == "detailed") {
      # Create a detailed mutation plot showing individual samples
      # Extract mutation positions
      mutations <- try({
        data.frame(
          sample = gene_data$Tumor_Sample_Barcode,
          position = as.numeric(gsub(".*p\\.[A-Za-z]*(\\d+).*", "\\1", gene_data$HGVSp_Short)),
          type = gene_data$Variant_Classification,
          change = gene_data$HGVSp_Short,
          vaf = ifelse(!is.null(gene_data$t_vaf), gene_data$t_vaf, NA),
          impact = gene_data$IMPACT
        )
      }, silent = TRUE)
      
      if(inherits(mutations, "try-error") || all(is.na(mutations$position))) {
        message(paste0("Cannot extract mutation positions for gene ", gene))
        next
      }
      
      # Filter valid positions
      mutations <- mutations[!is.na(mutations$position),]
      
      if(nrow(mutations) > 0) {
        # Summary of mutations
        pos_summary <- mutations %>%
          group_by(position, type) %>%
          summarise(count = n(), .groups = 'drop')
        
        # Create detailed chart
        p <- ggplot() +
          # Base protein line
          geom_segment(aes(x = min(mutations$position), xend = max(mutations$position),
                           y = 0, yend = 0), color = "gray50", size = 1) +
          # Mutation markers
          geom_point(data = pos_summary, 
                     aes(x = position, y = count, size = count, fill = type),
                     shape = 21, color = "black", alpha = 0.9) +
          # Mutation labels
          geom_text_repel(data = pos_summary %>% filter(count >= min_freq),
                          aes(x = position, y = count, label = paste0(position, " (", count, ")")),
                          size = 3, box.padding = 0.5, max.overlaps = 15) +
          # Styling
          labs(title = paste0("Detailed mutation map for ", gene),
               subtitle = paste0(nrow(mutations), " mutations in ", length(unique(mutations$sample)), " samples"),
               x = "Amino Acid Position",
               y = "Mutation Count",
               fill = "Variant Type") +
          theme_minimal() +
          theme(legend.position = "bottom",
                plot.title = element_text(face = "bold", hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5))
        
        result_list[[gene]] <- p
      } else {
        message("Cannot create plot: Unable to extract amino acid positions")
        result_list[[gene]] <- NULL
      }
    }
  }
  
  if(returnData) {
    return(list(plots = result_list, data = data_list))
  } else {
    return(result_list)
  }
}

#' Extract mutation details for specific genes
#'
#' @param maf MAF object from maftools
#' @param genes Character vector of gene names 
#' @param includeAll Whether to include all variants (TRUE) or only non-synonymous (FALSE)
#'
#' @return Data frame with mutation details
#'
extract_mutation_details <- function(maf, genes, includeAll = TRUE) {
  require(maftools)
  require(dplyr)
  
  all_mutations <- data.frame()
  
  for(gene in genes) {
    # Extract data for this gene
    if(includeAll) {
      gene_data <- maf@data[maf@data$Hugo_Symbol == gene,]
    } else {
      gene_data <- subsetMaf(maf, genes = gene)@data
    }
    
    if(nrow(gene_data) == 0) {
      message(paste0("No variants found for gene ", gene))
      next
    }
    
    # Process mutations to extract position info
    mutations <- try({
      mut_data <- data.frame(
        Gene = gene,
        Sample = gene_data$Tumor_Sample_Barcode,
        AA_Position = as.numeric(gsub(".*p\\.[A-Za-z]*(\\d+).*", "\\1", gene_data$HGVSp_Short)),
        HGVSp_Short = gene_data$HGVSp_Short,
        Variant_Class = gene_data$Variant_Classification,
        Variant_Type = gene_data$Variant_Type,
        Chromosome = gene_data$Chromosome,
        Start_Pos = gene_data$Start_Position,
        End_Pos = gene_data$End_Position
      )
      
      # Add VAF if available
      if("t_vaf" %in% colnames(gene_data)) {
        mut_data$VAF <- gene_data$t_vaf
      }
      
      # Add impact if available
      if("IMPACT" %in% colnames(gene_data)) {
        mut_data$Impact <- gene_data$IMPACT
      }
      
      mut_data
    }, silent = TRUE)
    
    if(!inherits(mutations, "try-error")) {
      all_mutations <- rbind(all_mutations, mutations)
    }
  }
  
  return(all_mutations)
}

#' Create a custom mutation lollipop plot with more control
#'
#' @param maf MAF object from maftools
#' @param gene Gene name to visualize
#' @param color_by Feature to color by (e.g., "Variant_Class")
#' @param label_top Top N mutations to label
#' @param repel Whether to use ggrepel for label positioning
#' @param show_domains Whether to show protein domains
#'
#' @return ggplot object
#'
create_custom_lollipop <- function(maf, gene, color_by = "Variant_Class", 
                                  label_top = 5, repel = TRUE, show_domains = TRUE) {
  require(maftools)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  # Extract mutations for the gene
  gene_data <- maf@data[maf@data$Hugo_Symbol == gene,]
  
  if(nrow(gene_data) == 0) {
    message(paste0("No mutations found for gene ", gene))
    return(NULL)
  }
  
  # Process mutations
  mutations <- try({
    data.frame(
      Sample = gene_data$Tumor_Sample_Barcode,
      Position = as.numeric(gsub(".*p\\.[A-Za-z]*(\\d+).*", "\\1", gene_data$HGVSp_Short)),
      Variant_Class = gene_data$Variant_Classification,
      Variant_Type = gene_data$Variant_Type,
      HGVSp_Short = gene_data$HGVSp_Short
    )
  }, silent = TRUE)
  
  if(inherits(mutations, "try-error") || all(is.na(mutations$Position))) {
    message("Cannot extract mutation positions")
    return(NULL)
  }
  
  # Filter valid positions
  mutations <- mutations[!is.na(mutations$Position),]
  
  if(nrow(mutations) == 0) {
    message("No valid positions found")
    return(NULL)
  }
  
  # Summarize mutations by position - avoid using n()
  pos_summary <- mutations %>%
    group_by(Position, Variant_Class) %>%
    summarise(
      Count = length(Sample),
      Samples = list(unique(Sample)), 
      Changes = list(unique(HGVSp_Short)), 
      .groups = 'drop'
    )
  
  # Get top mutations for labeling
  top_mutations <- pos_summary %>%
    arrange(desc(Count)) %>%
    head(label_top)
  
  # Create the plot
  p <- ggplot() +
    # Base protein line
    geom_segment(aes(x = min(mutations$Position), xend = max(mutations$Position),
                     y = 0, yend = 0), color = "gray50", size = 1) +
    # Lollipop stems
    geom_segment(data = pos_summary,
                 aes(x = Position, xend = Position, y = 0, yend = Count),
                 color = "gray30", size = 0.5) +
    # Lollipop heads
    geom_point(data = pos_summary,
               aes(x = Position, y = Count, fill = Variant_Class),
               shape = 21, color = "black", size = 3) +
    # Labels for top mutations
    labs(title = paste0("Mutation profile for ", gene),
         subtitle = paste0(nrow(mutations), " mutations in ", length(unique(mutations$Sample)), " samples"),
         x = "Amino Acid Position",
         y = "Mutation Count",
         fill = color_by) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Add labels with ggrepel if requested
  if(repel) {
    p <- p + geom_text_repel(data = top_mutations,
                            aes(x = Position, y = Count, 
                                label = paste0(Position, " (", Count, ")")),
                            size = 3, box.padding = 0.5, max.overlaps = 15)
  } else {
    p <- p + geom_text(data = top_mutations,
                      aes(x = Position, y = Count, 
                          label = paste0(Position, " (", Count, ")")),
                      size = 3, vjust = -0.5)
  }
  
  return(p)
}

#' Create a faceted mutation plot to compare multiple genes
#'
#' @param maf MAF object from maftools
#' @param genes Character vector of gene names to compare
#' @param color_by Feature to color by (e.g., "Variant_Class")
#' @param normalize Whether to normalize counts (percentage instead of absolute)
#'
#' @return ggplot object
#'
create_comparative_mutation_plot <- function(maf, genes, color_by = "Variant_Class", normalize = FALSE) {
  require(maftools)
  require(ggplot2)
  require(dplyr)
  
  all_mutations <- data.frame()
  
  for(gene in genes) {
    # Extract data for this gene
    gene_data <- maf@data[maf@data$Hugo_Symbol == gene,]
    
    if(nrow(gene_data) == 0) {
      message(paste0("No variants found for gene ", gene))
      next
    }
    
    # Process mutations to extract position info
    mutations <- try({
      mut_data <- data.frame(
        Gene = gene,
        Sample = gene_data$Tumor_Sample_Barcode,
        Position = as.numeric(gsub(".*p\\.[A-Za-z]*(\\d+).*", "\\1", gene_data$HGVSp_Short)),
        Variant_Class = gene_data$Variant_Classification,
        Variant_Type = gene_data$Variant_Type,
        HGVSp_Short = gene_data$HGVSp_Short
      )
      mut_data
    }, silent = TRUE)
    
    if(!inherits(mutations, "try-error")) {
      # Filter valid positions
      mutations <- mutations[!is.na(mutations$Position),]
      if(nrow(mutations) > 0) {
        all_mutations <- rbind(all_mutations, mutations)
      }
    }
  }
  
  if(nrow(all_mutations) == 0) {
    message("No valid mutations found for any gene")
    return(NULL)
  }
  
  # Summarize by gene and position
  if(normalize) {
    # Get total mutations per gene for normalization
    gene_totals <- all_mutations %>%
      group_by(Gene) %>%
      summarise(Total = n(), .groups = 'drop')
    
    pos_summary <- all_mutations %>%
      group_by(Gene, Position, Variant_Class) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      left_join(gene_totals, by = "Gene") %>%
      mutate(Percentage = Count / Total * 100)
    
    # Create plot with normalized values
    p <- ggplot(pos_summary, aes(x = Position, y = Percentage, fill = Variant_Class)) +
      geom_bar(stat = "identity", position = "stack", width = 5) +
      facet_wrap(~ Gene, scales = "free_x") +
      labs(title = "Comparative Mutation Profile",
           x = "Amino Acid Position",
           y = "Mutation Percentage (%)",
           fill = color_by) +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold"))
  } else {
    # Create plot with absolute counts
    pos_summary <- all_mutations %>%
      group_by(Gene, Position, Variant_Class) %>%
      summarise(Count = n(), .groups = 'drop')
    
    p <- ggplot(pos_summary, aes(x = Position, y = Count, fill = Variant_Class)) +
      geom_bar(stat = "identity", position = "stack", width = 5) +
      facet_wrap(~ Gene, scales = "free_x") +
      labs(title = "Comparative Mutation Profile",
           x = "Amino Acid Position",
           y = "Mutation Count",
           fill = color_by) +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold"))
  }
  
  return(p)
}

#' Create a detailed mutation map for a gene showing mutation hotspots
#'
#' @param maf MAF object from maftools
#' @param gene Gene name to visualize
#' @param color_by Feature to color by (e.g., "Variant_Classification")
#' @param point_size Base size for points
#' @param min_freq Minimum frequency to label a position
#'
#' @return ggplot object
#'
create_mutation_map <- function(maf, gene, color_by = "Variant_Classification", 
                              point_size = 3, min_freq = 2) {
  require(maftools)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  # Extract mutations for the gene
  gene_data <- maf@data[maf@data$Hugo_Symbol == gene,]
  
  if(nrow(gene_data) == 0) {
    message(paste0("No mutations found for gene ", gene))
    return(NULL)
  }
  
  # Check for protein length info
  prot_len <- NULL
  if(!is.null(maf@geneScheme) && gene %in% maf@geneScheme$gene) {
    prot_info <- maf@geneScheme[maf@geneScheme$gene == gene,]
    if(!is.null(prot_info$prot.length)) {
      prot_len <- prot_info$prot.length
    }
  }
  
  # Process mutations
  mutations <- try({
    data.frame(
      Sample = gene_data$Tumor_Sample_Barcode,
      Position = as.numeric(gsub(".*p\\.[A-Za-z]*(\\d+).*", "\\1", gene_data$HGVSp_Short)),
      Variant_Classification = gene_data$Variant_Classification,
      Variant_Type = gene_data$Variant_Type,
      HGVSp_Short = gene_data$HGVSp_Short
    )
  }, silent = TRUE)
  
  if(inherits(mutations, "try-error") || all(is.na(mutations$Position))) {
    message("Cannot extract mutation positions")
    return(NULL)
  }
  
  # Filter valid positions
  mutations <- mutations[!is.na(mutations$Position),]
  
  if(nrow(mutations) == 0) {
    message("No valid positions found")
    return(NULL)
  }
  
  # Summarize mutations by position
  pos_summary <- mutations %>%
    group_by(Position, Variant_Classification) %>%
    summarise(Count = n(), Changes = list(unique(HGVSp_Short)), .groups = 'drop')
  
  # Create the plot
  p <- ggplot() +
    # Base protein line
    geom_segment(aes(x = 0, 
                     xend = ifelse(is.null(prot_len), 
                                  max(mutations$Position) + 50, 
                                  prot_len),
                     y = 0, yend = 0), 
                color = "gray50", size = 1) +
    # Mutation markers
    geom_point(data = pos_summary,
               aes(x = Position, y = Count, size = Count, fill = Variant_Classification),
               shape = 21, color = "black", alpha = 0.9) +
    # Labels for frequent mutations
    geom_text_repel(data = pos_summary %>% filter(Count >= min_freq),
                   aes(x = Position, y = Count, 
                       label = paste0(Position, " (", Count, ")")),
                   size = 3, box.padding = 0.5, max.overlaps = 30) +
    # Styling
    scale_size_continuous(range = c(point_size, point_size*3)) +
    labs(title = paste0("Mutation map for ", gene),
         subtitle = paste0(nrow(mutations), " mutations in ", 
                          length(unique(mutations$Sample)), " samples"),
         x = "Amino Acid Position",
         y = "Mutation Count",
         fill = color_by,
         size = "Count") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}

#' Create a mutation heatmap showing mutation frequency across genes
#'
#' @param maf MAF object from maftools
#' @param genes Character vector of gene names
#' @param color_palette Color palette for the heatmap
#'
#' @return ggplot object
#'
create_mutation_heatmap <- function(maf, genes, color_palette = NULL) {
  require(maftools)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  # Initialize data frame to store mutation counts
  counts_df <- data.frame(Sample = unique(maf@data$Tumor_Sample_Barcode))
  
  # Extract mutation counts for each gene and sample
  for(gene in genes) {
    # Get mutations for this gene
    gene_data <- maf@data[maf@data$Hugo_Symbol == gene,]
    
    if(nrow(gene_data) == 0) {
      # No mutations found, add column of zeros
      counts_df[[gene]] <- 0
      next
    }
    
    # Count mutations per sample
    gene_counts <- table(gene_data$Tumor_Sample_Barcode)
    
    # Add to data frame
    counts_df[[gene]] <- 0
    for(sample in names(gene_counts)) {
      if(sample %in% counts_df$Sample) {
        counts_df[counts_df$Sample == sample, gene] <- gene_counts[sample]
      }
    }
  }
  
  # Convert to long format for plotting
  long_df <- counts_df %>%
    pivot_longer(cols = -Sample, names_to = "Gene", values_to = "Mutations")
  
  # Create the heatmap
  p <- ggplot(long_df, aes(x = Gene, y = Sample, fill = Mutations)) +
    geom_tile() +
    scale_fill_gradientn(colors = if(is.null(color_palette)) 
                                   c("white", "yellow", "orange", "red", "darkred") 
                                 else color_palette) +
    labs(title = "Mutation Frequency Heatmap",
         x = "Gene",
         y = "Sample",
         fill = "Mutations") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  return(p)
}

# Example usage:
# 
# # Load maftools and your MAF file
# library(maftools)
# maf <- read.maf("your_maf_file.maf")
# 
# # Create a custom lollipop plot for a specific gene
# create_custom_lollipop(maf, "TP53")
# 
# # Create detailed mutation maps for multiple genes
# genes_of_interest <- c("TP53", "KRAS", "EGFR")
# mutation_maps <- visualize_gene_mutations(maf, genes_of_interest, plotType = "detailed")
# 
# # Create a comparative plot showing mutation distributions across genes
# comp_plot <- create_comparative_mutation_plot(maf, genes_of_interest)
# 
# # Extract detailed mutation information
# mutation_details <- extract_mutation_details(maf, genes_of_interest)
# 
# # Create a mutation heatmap
# heatmap <- create_mutation_heatmap(maf, genes_of_interest)
# 
# # Visualize the plots
# print(mutation_maps$TP53)
# print(comp_plot)
# print(heatmap)
# 
# # Analyze hotspots across multiple genes
# for(gene in genes_of_interest) {
#   map <- create_mutation_map(maf, gene, min_freq = 3)
#   print(map)
# }
# 
# # Save results to PDF
# pdf("mutation_analysis_results.pdf", width = 10, height = 8)
# 
# # Save all visualizations
# for(gene in names(mutation_maps)) {
#   if(!is.null(mutation_maps[[gene]])) {
#     print(mutation_maps[[gene]])
#   }
# }
# 
# print(comp_plot)
# print(heatmap)
# 
# # Generate lollipop plots with domains for key genes
# for(gene in genes_of_interest) {
#   lollipop <- visualize_gene_mutations(maf, gene, plotType = "lollipop", showDomains = TRUE)
#   if(!is.null(lollipop[[gene]])) {
#     print(lollipop[[gene]])
#   }
# }
# 
# # Create oncostrips for the genes
# onco_plots <- visualize_gene_mutations(maf, genes_of_interest, plotType = "oncostrip")
# for(gene in names(onco_plots)) {
#   if(!is.null(onco_plots[[gene]])) {
#     print(onco_plots[[gene]])
#   }
# }
# 
# # Close PDF device
# dev.off()
# 
# # Export mutation details to CSV
# write.csv(mutation_details, "mutation_details.csv", row.names = FALSE)
# 
# # Create a summary report
# mutation_summary <- mutation_details %>%
#   group_by(Gene, Variant_Class) %>%
#   summarise(Count = n(), 
#             Samples = n_distinct(Sample),
#             .groups = 'drop')
# 
# # Print summary
# print(mutation_summary)
# 
# # Optional - Create additional custom visualizations for specific purposes
# # For example, a normalized comparison plot
# norm_comp <- create_comparative_mutation_plot(maf, genes_of_interest, normalize = TRUE)
# print(norm_comp)
#
# # End of example script
