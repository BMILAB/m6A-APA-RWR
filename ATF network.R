library(dplyr)
library(tidyverse)

## ----- Load QAPA results  -----
data <- read.delim("QAPAresults.txt")

# Compute sample expression
TPM <- data %>%
  group_by(Gene, Gene_Name) %>%
  summarise(
    across(
      .cols = starts_with("sample"), 
      .fns = sum,                   
      .names = "Gene_{.col}"         
    ),
    .groups = "drop"    
  )

# Merge expression data
df_APA <- alkbh10b_APA %>%
  inner_join(alkbh10b_APA_flower_0h, by = c("gene")) %>%
  inner_join(alkbh10b_APA_flower_3h, by = c("gene")) %>%
  inner_join(alkbh10b_APA_leaf_0h, by = c("gene")) %>%
  inner_join(alkbh10b_APA_leaf_3h, by = c("gene")) %>%
  inner_join(alkbh9b_APA_ABA, by = c("gene")) %>%
  inner_join(alkbh9b_APA_MS0, by = c("gene")) %>%
  inner_join(cold_APA, by = c("gene")) %>%
  inner_join(fiona1_APA, by = c("gene")) %>%
  inner_join(mta_APA_chilling, by = c("gene")) %>%
  inner_join(mta_APA_control, by = c("gene")) %>%
  inner_join(vir_APA_0h, by = c("gene")) %>%
  inner_join(vir_APA_4h, by = c("gene"))

# Arabidopsis poly(A) factor
ath_factor <- c("AT1G30460", "AT5G23880", "AT5G58040", "AT1G61010", "AT3G06850", "AT5G60940", 
                "AT1G71820", "AT5G60950", "AT3G04600", "AT5G58040", "AT1G17960", "AT1G17980", 
                "AT2G36320", "AT3G06560", "AT5G57610", "AT1G30460", "AT1G49760", "AT5G04290")

factor_TPM <- tpm[tpm$gene %in% ath_factor, ]

# Calculate log₂(TPM + 1) for the TPM matrix
factor_tpm <- log2(as.matrix(factor_TPM[, -(1:2)]) + 1)

noreg_apa <- DE_apa[rownames(DE_apa) %in% no_reg$TAIR_ID, ]

## ----- Calculate the correlation between all factors and APA events  -----
# Initialize the result dataframe
ATF_APA_cor_results <- tibble(
  APA_Factor = character(),
  APA_Event = character(),
  Spearman_Rho = numeric(),
  P_value = numeric()
)

factors <- rownames(factor_tpm)
events <- rownames(df_APA)

# Traverse all combinations via nested for loops 
for (factor in factors) {
  for (event in events) {
    
    x <- as.numeric(factor_tpm[factor, ])  
    y <- as.numeric(df_APA[event, ])   
    
    valid_idx <- complete.cases(x, y)
    x_valid <- x[valid_idx]
    y_valid <- y[valid_idx]
    
    if (length(x_valid) < 2) {
      rho <- NA_real_
      p_val <- NA_real_
    } else {
      test <- cor.test(x_valid, y_valid, method = "spearman")
      rho <- test$estimate
      p_val <- test$p.value
    }
    
    ATF_APA_cor_results <- bind_rows(ATF_APA_cor_results, tibble(
      APA_Factor = factor,
      APA_Event = event,
      Spearman_Rho = rho,
      P_value = p_val
    ))
  }
}


## ----- Plot the ATF-m6APAreg regulatory network  -----
# The following data is available:
# 1. factor_apa: Correlation data between APA factors and APA events
#    Columns included: Source (APA factor), Target (APA event gene), Correlation, Pvalue
# 2. ppi_edges: PPI interaction data between APA event genes
#    Columns included: Source (Gene 1), Target (Gene 2)
# Load packages
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(tidygraph)
library(scales)
library(ggrepel)

# edges data
all_edges <- bind_rows(
  factor_APA %>% 
    mutate(Interaction = "Correlation",
           Correlation_Type = ifelse(Correlation > 0, "Positive", "Negative")),
  
  ppi_edges %>% 
    mutate(Interaction = "PPI")
)

# nodes data
all_nodes <- data.frame(
  ID = unique(c(all_edges$Source, all_edges$Target))
) %>%
  mutate(
    SubType = case_when(
      ID %in% factor_APA$Source ~ "APA factor",
      ID %in% reg_gene_big ~ "m6APAreg(pcc>0)",
      ID %in% reg_gene_small ~ "m6APAreg(pcc<0)",
      TRUE ~ "non-m6APAreg"
    ),
    DisplayName = ID
  )

# Create the complete network graph
g_full <- graph_from_data_frame(all_edges, vertices = all_nodes)

# Filter out the "Gene_Other" nodes
nodes_to_keep <- V(g_full)$name[V(g_full)$SubType != "NA"]
g <- induced_subgraph(g_full, vids = nodes_to_keep)

# Node centrality
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)

# Convert to a tidygraph object
g_tidy <- as_tbl_graph(g)

# Create a layout
set.seed(123)
layout <- create_layout(g_tidy, layout = "fr")

# Generate the network plot
p <- ggraph(layout) +
  geom_edge_link(
    aes(color = case_when(
      Interaction == "Correlation" & Correlation_Type == "Positive" ~ "Positive",
      Interaction == "Correlation" & Correlation_Type == "Negative" ~ "Negative",
      Interaction == "PPI" ~ "PPI"
    )),
    width = 0.7,
    alpha = 0.7,
    show.legend = TRUE
  ) +
  
  geom_node_point(
    aes(fill = SubType, 
        size = ifelse(SubType == "APA factor", degree * 4, degree * 4)),  
    shape = 21, 
    color = "grey20",
    stroke = 0.6,
    alpha = 0.9
  ) +
  
  scale_size_area(max_size = 15) +
  
  geom_text_repel(
    aes(x = x, y = y, 
        label = ifelse(SubType == "APA factor", DisplayName, "")),
    size = 4,
    color = "#f08f92",
    fontface = "bold",
    min.segment.length = 0,
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "grey50",
    family = "sans",
    force = 0.5,
    max.time = 1,
    max.iter = 10000,
    max.overlaps = Inf
  ) +
  
  geom_text_repel(
    aes(x = x, y = y, 
        label = ifelse(SubType == "m6APAreg(pcc>0)" & degree > quantile(degree,0.9), DisplayName, "")),
    size = 3.0,
    color = "grey30",
    min.segment.length = 0.1,
    box.padding = 0.2,
    point.padding = 0.2,
    max.overlaps = 1, 
    segment.color = "grey70",
    segment.alpha = 0.5,
    family = "sans",
    force = 0.3
  ) + 
  geom_text_repel(
    aes(x = x, y = y, 
        label = ifelse(SubType == "m6APAreg(pcc<0)" & degree > quantile(degree, 0.9), DisplayName, "")),
    size = 3.0,
    color = "grey30",
    min.segment.length = 0.1,
    box.padding = 0.2,
    point.padding = 0.2,
    max.overlaps = 1,
    segment.color = "grey70",
    segment.alpha = 0.5,
    family = "sans",
    force = 0.3
  ) +
  
  geom_text_repel(
    aes(x = x, y = y, 
        label = ifelse(SubType == "non-m6APAreg" & degree > quantile(degree, 0.95), DisplayName, "")),
    size = 3.0,
    color = "grey30",
    min.segment.length = 0.1,
    box.padding = 0.2,
    point.padding = 0.2,
    max.overlaps = 1, 
    segment.color = "grey70",
    segment.alpha = 0.5,
    family = "sans",
    force = 0.3
  ) +

  scale_edge_color_manual(
    name = "Edge Type",
    values = c("Positive" = "#f08f92", 
               "Negative" = "#9cbedb", 
               "PPI" = "#a9d5a5"),
    labels = c("Positive Correlation", 
               "Negative Correlation", 
               "Protein Interaction"),
    breaks = c("Positive", "Negative", "PPI")
  ) +
  
  scale_fill_manual(
    name = "Node Type",
    values = c(
      "APA factor" = "#f08f92", 
      "m6APAreg(pcc>0)" = "#9cbedb",  
      "m6APAreg(pcc<0)" = "#a9d5a5", 
      "non-m6APAreg" = "#cccccc"     
    ),
    labels = c("APA factor", "m6APAreg(pcc>0)", "m6APAreg(pcc<0)", "non-m6APAreg"),
    breaks = c("APA factor", "m6APAreg(pcc>0)", "m6APAreg(pcc<0)", "non-m6APAreg")
  ) +
  
  scale_size_continuous(guide = "none") +
  
  theme_void() +
  theme(
    legend.position = c(0.96, 0.01),
    legend.justification = c(1, 0),
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = "grey80", size = 0.3),
    legend.margin = margin(3, 5, 5, 5),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    plot.background = element_rect(fill = "#F8F8F8", color = NA),
    plot.margin = margin(0.6, 0.6, 0.6, 0.6, "cm")
  ) +
  
  guides(
    edge_color = guide_legend(
      title.position = "top",
      order = 1,
      override.aes = list(edge_width = 1.5)
    ),
    fill = guide_legend(
      title.position = "top",
      order = 2,
      override.aes = list(size = 3, shape = 21)
    )
  )

# Save as PDF
pdf("ATF_m6APAreg.pdf", width = 10, height = 8, useDingbats = FALSE)
print(p)
dev.off()




