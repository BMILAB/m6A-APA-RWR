library(dplyr)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(patchwork)

## ----- Load m6A and APA data  -----
apa <- read.csv("apa_list_DE.csv")
m6a <- read.csv("m6a_list_DE.csv")

#Define sample groups
conditions <- data.frame(
  con = colnames(m6a),
  sample = c(
    rep("ALKBH10B", 4),
    rep("ALKBH10B-leaf-flower", 16),
    rep("ALKBH9B", 8),
    rep("fiona1", 6),
    rep("cold", 6),
    rep("MTA", 8),
    rep("vir", 8)
  )
)

corr_results <- matrix(NA, nrow = length(high_corr_genes), ncol = length(unique(conditions)),
                       dimnames = list(high_corr_genes, unique(conditions)))

condition_groups <- unique(conditions$sample)

corr_results <- matrix(NA, nrow = nrow(m6a), ncol = length(condition_groups),
                       dimnames = list(rownames(m6a), condition_groups))


## ----- Calculate m⁶A-APA correlation per gene per condition  ----- 
for(gene in rownames(m6a)){
  for(cond in condition_groups){
    
    # Get samples for current condition
    cond_samples <- conditions$con[conditions$sample == cond]
    
    cond_samples <- intersect(cond_samples, colnames(m6a))
    
    if(length(cond_samples) == 0) next
    
    # Get m⁶A & APA data for this gene and condition
    m6A_vals <- m6a[gene, cond_samples]
    apa_vals <- apa[gene, cond_samples]
    
    # Get non-NA indices
    valid_idx <- !is.na(m6A_vals) & !is.na(apa_vals)
    m6A_valid <- m6A_vals[valid_idx]
    apa_valid <- apa_vals[valid_idx]  
    
    # Ensure enough samples for correlation
    if(sum(valid_idx) >= 3) {
      # Calculate Spearman correlation
      corr_val <- cor(m6A_valid, apa_valid, method = "spearman")
    } else {
      corr_val <- NA
    }
    
    # 存储结果
    corr_results[gene, cond] <- corr_val
  }
}

## ----- Plot correlation heatmap with clustering  ----- 
pheatmap_obj <- pheatmap(corr_subset,
                         color = colorRampPalette(c("#AFC9CF","#D5E1E3","#EBBFC2","#E28187","#D93F49"))(100),
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_rownames = FALSE,
                         main = paste0("Condition-specific m6A-APA Correlation\n(Genes with data in ≥4 conditions; NA: ", na_percentage, "%)"),
                         border_color = NA,
                         fontsize_col = 10,
                         na_col = "gray90",
                         clustering_method = "ward.D2",
                         angle_col = 45,
                         silent = TRUE)

ggsave("heatmap(final).pdf", plot = pheatmap_obj, width = 8, height = 8, units = "in")

## ----- Group genes by clustering results  -----  
# Extract row dendrogram
row_dend <- as.hclust(pheatmap_obj$tree_row)

# Show dendrogram to choose k
plot(row_dend, cex = 0.6, main = "Gene Clustering Dendrogram")
rect.hclust(row_dend, h = 1.5, border = 2:5)
k <- 5
gene_clusters <- cutree(row_dend, k = k)

# Convert cluster results to dataframe
cluster_df <- data.frame(
  Gene = names(gene_clusters),
  Cluster = gene_clusters
)

# Extract genes by group
gene_groups <- split(cluster_df$Gene, cluster_df$Cluster)

plot_data <- corr_results_1 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Condition",
    values_to = "Correlation"
  ) %>%
  left_join(cluster_df, by = "Gene") %>%
  mutate(Cluster = factor(Cluster))

# Calculate mean correlation per group (NA removed)
group_means <- plot_data %>%
  group_by(Cluster, Condition) %>%
  summarise(
    Mean_Correlation = mean(Correlation, na.rm = TRUE),
    SE = sd(Correlation, na.rm = TRUE) / sqrt(sum(!is.na(Correlation))),
    .groups = "drop"
  )

cluster_colors <- c("#8386a8","#d15c6b","#f5cf36","#8fb943","#78b9d2")

# Define plot function
create_cluster_plot <- function(cluster_id) {
  cluster_data <- plot_data %>% filter(Cluster == cluster_id)
  cluster_means <- group_means %>% filter(Cluster == cluster_id)
  
  n_genes <- length(unique(cluster_data$Gene))
  
  ggplot() +
    geom_line(
      data = cluster_data,
      aes(x = Condition, y = Correlation, group = Gene),
      color = "gray70",
      alpha = 0.3,
      size = 0.5,
      na.rm = TRUE
    ) +

    geom_line(
      data = cluster_means,
      aes(x = Condition, y = Mean_Correlation, group = 1),
      color = cluster_colors[cluster_id],
      size = 2,
      na.rm = TRUE
    ) +

    geom_errorbar(
      data = cluster_means,
      aes(
        x = Condition,
        ymin = Mean_Correlation - SE,
        ymax = Mean_Correlation + SE
      ),
      width = 0.2,
      color = "black",
      size = 0.8
    ) +
    
    geom_point(
      data = cluster_means,
      aes(x = Condition, y = Mean_Correlation),
      color = "white",
      size = 4,
      shape = 21,
      fill = cluster_colors[cluster_id],
      stroke = 1.5
    ) +
    labs(
      title = paste("Cluster", cluster_id),
      subtitle = paste(n_genes, "genes"),
      y = "Correlation"
    ) +
    ylim(-1, 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(
        face = "bold",
        size = 14,
        color = cluster_colors[cluster_id]
      ),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
    )
}

# Generate group-specific figures
plot_list <- map(1:5, create_cluster_plot)

# Combine all plots
p <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title = "m6A-APA Correlation Patterns by Gene Cluster",
    subtitle = "Gray lines: individual genes | Colored line: cluster mean ± SE",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30")
    )
  )

# Save
pdf("cluster_correlation_profiles.pdf", width = 12, height = 16)
print(p)
dev.off()
