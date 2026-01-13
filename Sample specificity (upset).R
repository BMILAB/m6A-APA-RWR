library(ggplot2)
library(UpSetR)
library(clusterProfiler)
library(org.At.tair.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(stringr)

## ----- Load m6APAreg genes across samples -----
gene_lists1 <- list(
  ALKBH10B = sig_alkbh10b,
  ALKBH10B_leaf_flower = sig_alkbh10b_heat,
  ALKBH9B = sig_alkbh9b,
  cold = sig_cold,
  fiona1 = sig_fiona1,
  MTA = sig_mta,
  vir = sig_vir
)

gene_lists2 <- list(
  ALKBH10B = SIG_alkbh10b,
  ALKBH10B_leaf_flower = SIG_alkbh10b_heat,
  ALKBH9B = SIG_alkbh9b,
  cold = SIG_cold,
  fiona1 = SIG_fiona1,
  MTA = SIG_mta,
  vir = SIG_vir
)

## ----- upset plot -----
p1 <- upset(fromList(gene_lists2),
            nsets = 7,
            nintersects = 25,
            order.by = "freq",
            decreasing = TRUE,
            mb.ratio = c(0.6, 0.4),
            number.angles = 30,
            point.size = 3.5,
            line.size = 1.5,
            mainbar.y.label = "Intersection Size",
            sets.x.label = "Set Size",
            text.scale = c(1.5, 1.5, 1.2, 1.2, 1.8, 1.2),
            sets.bar.color = "#53A6D9",
            main.bar.color = "#78aac8",
            matrix.color = "#704739"
)

p2 <- upset(fromList(gene_lists1),
            nsets = 7,
            nintersects = 25,
            order.by = "freq",
            decreasing = TRUE,
            mb.ratio = c(0.6, 0.4),
            number.angles = 30,
            point.size = 3.5,
            line.size = 1.5,
            mainbar.y.label = "Intersection Size",
            sets.x.label = "Set Size",
            text.scale = c(1.5, 1.5, 1.2, 1.2, 1.8, 1.2),
            sets.bar.color = "#7dba7f",
            main.bar.color = "#92C6A0",
            matrix.color = "#704739"
)


## ----- Extract sample-specific genes (unique to each sample) -----
specific_genes1 <- list()
for(sample_name in names(gene_lists1)) {
  current_genes <- gene_lists1[[sample_name]]
  other_genes <- unlist(gene_lists1[names(gene_lists1) != sample_name])
  specific_genes1[[sample_name]] <- setdiff(current_genes, other_genes)
}

sapply(specific_genes1, length)

specific_genes2 <- list()
for(sample_name in names(gene_lists2)) {
  current_genes <- gene_lists2[[sample_name]]
  other_genes <- unlist(gene_lists2[names(gene_lists2) != sample_name])
  specific_genes2[[sample_name]] <- setdiff(current_genes, other_genes)
}


# Add prefixes for distinction
names(specific_genes1) <- paste0("Set1_", names(specific_genes1))
names(specific_genes2) <- paste0("Set2_", names(specific_genes2))

# Define plot function
perform_GO_enrichment <- function(gene_list, sample_name) {
  if(length(gene_list) == 0) {
    return(NULL)
  }
  
  GO <- enrichGO(
    gene = gene_list,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    pAdjustMethod = "BH",
    minGSSize = 1
  )
  
  if(is.null(GO) || nrow(GO@result) == 0) {
    return(NULL)
  }
  
  result_df <- as.data.frame(GO@result)
  result_df$Sample <- sample_name
  result_df$Set <- ifelse(grepl("^Set1", sample_name), "m6APAreg(pcc>0)", "m6APAreg(pcc<0)")
  
  return(result_df)
}

# GO enrichment for Set1
enrich_results_set1 <- list()
for(sample_name in names(specific_genes1)) {
  result <- perform_GO_enrichment(specific_genes1[[sample_name]], sample_name)
  if(!is.null(result)) {
    enrich_results_set1[[sample_name]] <- result
  }
}

# Remove samples with <5 enriched pathways
filtered_results_set1 <- list()
for(sample_name in names(enrich_results_set1)) {
  sample_result <- enrich_results_set1[[sample_name]]
  sig_terms <- sum(sample_result$p.adjust < 0.05, na.rm = TRUE)
  if(sig_terms >= 5) {
    filtered_results_set1[[sample_name]] <- sample_result
  }
}

# GO enrichment for Set2
enrich_results_set2 <- list()
for(sample_name in names(specific_genes2)) {
  result <- perform_GO_enrichment(specific_genes2[[sample_name]], sample_name)
  if(!is.null(result)) {
    enrich_results_set2[[sample_name]] <- result
  }
}

# Remove samples with <5 enriched pathways
filtered_results_set2 <- list()
for(sample_name in names(enrich_results_set2)) {
  sample_result <- enrich_results_set2[[sample_name]]
  sig_terms <- sum(sample_result$p.adjust < 0.05, na.rm = TRUE)
  if(sig_terms >= 5) {
    filtered_results_set2[[sample_name]] <- sample_result
  }
}


## ----- Define plot function -----
create_GO_plot <- function(filtered_results, set_name, set_label) {
  
  if(length(filtered_results) == 0) {
    warning(paste(set_name, "No samples meet the criterion (number of enriched pathways ≥ 5); skip plotting."))
    return(NULL)
  }
  
  all_results <- do.call(rbind, filtered_results)
  all_results$log_pvalue <- -log10(all_results$pvalue)
  
  top_terms <- all_results %>%
    filter(p.adjust < 0.05) %>%
    group_by(Sample) %>%
    slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  if(nrow(top_terms) == 0) {
    warning(paste(set_name, "No significant GO terms; skip plotting."))
    return(NULL)
  }
  
  term_order <- top_terms %>%
    group_by(Description) %>%
    summarize(minP = min(p.adjust)) %>%
    arrange(minP) %>%
    pull(Description)
  
  top_terms$Description <- factor(top_terms$Description, levels = term_order)
  top_terms$Sample <- factor(top_terms$Sample, levels = names(filtered_results))
  
  if(set_name == "Set1") {
    color_low <- "#B5D3DC"  
    color_high <- "#4AA4B7"
    plot_title <- paste("GO Enrichment Analysis -", set_label)
  } else {
    color_low <- "#C7D5C3"   
    color_high <- "#85A982"  
    plot_title <- paste("GO Enrichment Analysis -", set_label)
  }
  
# plot
p <- ggplot(top_terms, aes(x = Sample, y = Description)) +
    geom_point(aes(size = Count, fill = log_pvalue), shape = 21, color = "black", stroke = 0.5) +
    scale_size_continuous(
      name = "Gene count",
      range = c(3, 8),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_fill_gradient(
      name = expression(-log[10]("P.adj")),
      low = color_low, 
      high = color_high,
      guide = guide_colorbar(reverse = FALSE)
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    labs(
      x = "Sample",
      y = "GO Biological Process",
      title = plot_title
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10, color = "black", lineheight = 0.8),
      axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 0),
      axis.title = element_text(size = 12, face = "bold"),
      axis.title.x.top = element_text(margin = margin(b = 10)),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.2, "cm"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, hjust = 0.5),
      panel.grid.major = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70"),
      plot.margin = margin(1, 1, 1, 1.5, "cm")
    )
  
  return(p)
}

p1 <- create_GO_plot(filtered_results_set1, "Set1", "m6APAreg (pcc > 0)")

p2 <- create_GO_plot(filtered_results_set2, "Set2", "m6APAreg (pcc < 0)")


