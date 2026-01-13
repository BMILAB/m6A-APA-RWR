library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(tidygraph)
library(scales)
library(ggrepel)

## ----- Load PPI data  -----
ppi_edges <- reg_net_sig %>% 
  rename(Source = protein1, Target = protein2) %>%
  mutate(Interaction = "PPI")

## ----- Load miRNA regulatory data  -----
mirna_edges <- miRNA_gene %>% 
  rename(Source = miRNA_Acc., Target = UniProt_ID) %>%
  mutate(Interaction = "miRNA_regulation") %>%
  distinct(Source, Target, .keep_all = TRUE)

## ----- Select top 15 significant miRNAs (optional)  -----
miRNA_counts <- mirna_edges %>%
  group_by(Source) %>%
  summarise(target_count = n()) %>%
  arrange(desc(target_count)) %>%
  head(15)

mirna_edges_filtered <- mirna_edges %>%
  filter(Source %in% miRNA_counts$Source)

## ----- Combine data -----
all_edges <- bind_rows(ppi_edges, mirna_edges_filtered)

# Create node list
all_nodes <- data.frame(
  ID = unique(c(all_edges$Source, all_edges$Target))
) %>%
  dplyr::mutate(
    Type = dplyr::case_when(
      grepl("ath-miR", ID) ~ "miRNA",
      grepl("^[A-Z0-9]{6,10}$", ID) ~ "m6APAreg",
      TRUE ~ "Other"
    ),
    # Add readable names
    DisplayName = ifelse(Type == "miRNA", ID, ID)
  )

# Separate m6APAreg genes into PCC > 0 and PCC < 0 categories
all_nodes <- data.frame(
  ID = unique(c(all_edges$Source, all_edges$Target))
) %>%
  dplyr::mutate(
    Type = dplyr::case_when(
      grepl("ath-miR", ID) ~ "miRNA",
      grepl("^[A-Z0-9]{6,10}$", ID) ~ "m6APAreg",
      TRUE ~ "Other"
    ),
    # Annotate gene categories
    SubType = case_when(
      Type == "miRNA" ~ "miRNA",
      ID %in% reg_gene_big ~ "m6APAreg(pcc>0)",
      ID %in% reg_gene_small ~ "m6APAreg(pcc<0)",
      TRUE ~ "Gene_Other"
    ),
    DisplayName = ifelse(Type == "miRNA", ID, ID)
  )

## ----- Create graph -----
g <- graph_from_data_frame(all_edges, vertices = all_nodes)

g <- induced_subgraph(
  g,
  vids = subcomponent(g, v = V(g)[Type == "miRNA"], mode = "all")
)

# Calculate centrality
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)

# Convert to tidygraph
g_tidy <- as_tbl_graph(g)

# Generate layout
set.seed(123)
layout <- create_layout(g_tidy, layout = "fr")

## ----- Plot network -----
p <- ggraph(layout) +
  # Edge settings
  geom_edge_link(
    aes(color = Interaction, 
        linetype = Interaction),
    width = 0.7,
    alpha = 0.6,
    show.legend = TRUE
  ) +
  
  # Node settings
  geom_node_point(
    aes(fill = SubType, 
        size = betweenness * 100),
    shape = 21,
    color = "grey20",
    stroke = 0.5,
    alpha = 0.95
  ) +
  
  # Highlight miRNA nodes
  geom_node_point(
    aes(filter = Type == "miRNA"),
    shape = 21,
    color = "#f08f92",
    fill = "#f08f92",
    size = 5,
    stroke = 1.2
  ) +
  
  # Label settings
  geom_text_repel(
    aes(x = x, y = y, 
        label = ifelse(Type == "miRNA", DisplayName, "")),
    size = 3.5,
    color = "#f08f92",
    fontface = "bold",
    min.segment.length = 0.5, 
    box.padding = 0.7, 
    point.padding = 0.7,
    max.overlaps = Inf,
    segment.color = "grey70",
    segment.alpha = 0.3,
    segment.size = 0.3,
    nudge_y = 0.3, 
    direction = "both",
    family = "sans"
  ) +
  
  geom_text_repel(
    aes(x = x, y = y, 
        label = ifelse(Type == "m6APAreg" & degree > quantile(degree, 0.3), DisplayName, "")),
    size = 2.8,
    color = "grey30",
    min.segment.length = 0.1,
    box.padding = 0.25,
    point.padding = 0.3,
    max.overlaps = 30,     
    segment.color = "grey70",
    segment.alpha = 0.4,
    segment.size = 0.3,
    family = "sans"
  ) +
  
  # Legend settings
  scale_edge_color_manual(
    name = "Interaction Type",
    values = c("PPI" = "#c2e5cf", "miRNA_regulation" = "#f2b8ae"),
    labels = c("Protein Interaction", "miRNA Regulation"),
    breaks = c("PPI", "miRNA_regulation")
  ) +
  
  scale_edge_linetype_manual(
    name = "Interaction Type",
    values = c("PPI" = "solid", "miRNA_regulation" = "solid"),
    guide = "none"
  ) +
  
  scale_fill_manual(
    name = "Node Type",
    values = c(
      "miRNA" = "#f08f92", 
      "m6APAreg(pcc>0)" = "#9cbedb",  
      "m6APAreg(pcc<0)" = "#a9d5a5",
      "Gene_Other" = "#74a9c5"
    ),
    labels = c("miRNA","m6APAreg(pcc>0)","m6APAreg(pcc<0)", "Other Gene"),
    breaks = c("miRNA","m6APAreg(pcc>0)","m6APAreg(pcc<0)", "Other Gene")
  ) +
  
  scale_size_continuous(
    name = "Node Centrality",
    range = c(3, 12),
    guide = "none"
  ) +
  
  # Theme settings
  theme_void() +
  theme(
    legend.position = c(0.98, 0.03),
    legend.justification = c(1, 0),
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = "grey80", size = 0.3),
    legend.margin = margin(3, 5, 5, 5),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    plot.background = element_rect(fill = "#F8F8F8", color = NA),
    plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm")
  ) +
  
  # Guides settings
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
pdf("miRNA-m6APAreg.pdf", width = 14, height = 12, useDingbats = FALSE)
print(p)
dev.off()
