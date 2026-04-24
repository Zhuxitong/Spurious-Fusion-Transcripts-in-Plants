############################################################
# 05. Fusion Network Analysis
#
# Purpose:
# Construct gene-level fusion networks from putative fusion
# transcript pairs and identify network modules.
#
# Related figures:
# Fig. 4a-e
#
# Input:
#   mh63_flagleaf.pc
#
# Required R package:
#   igraph
#
# Output:
#   isoseq_gene_modules.txt
#   isoseq_gene_modules_stats.txt
############################################################


library(igraph)
# isoseq network ----------------------------------------------------------

isoseq_gene_pairs <- read.table('mh63_flagleaf.pc', header = F) %>% select(V1,V2)
isoseq_g <- graph_from_data_frame(isoseq_gene_pairs, directed = F)
isoseq_communities <- cluster_louvain(isoseq_g)
print(isoseq_communities)
isoseq_module_ids <- membership(isoseq_communities)
isoseq_module_df <- data.frame(
  gene = names(isoseq_module_ids),
  module = as.numeric(isoseq_module_ids)
)
write_tsv(isoseq_module_df, "isoseq_gene_modules.txt")
isoseq_module_size <- sizes(isoseq_communities)
isoseq_module_size_df <- data.frame(
  module_id = as.integer(names(isoseq_module_size)),  
  gene_count = as.vector(isoseq_module_size)          
)
isoseq_module_avg_degree <- data.frame(
  gene = V(isoseq_g)$name,
  degree = degree(isoseq_g),
  module_id = isoseq_module_ids
) %>%
  group_by(module_id) %>%
  summarise(avg_degree = mean(degree))

isoseq_module_stats <- merge(isoseq_module_size_df, isoseq_module_avg_degree, by = 'module_id')
isoseq_module_stats <- isoseq_module_stats[order(-isoseq_module_stats$gene_count), ]
write_tsv(isoseq_module_stats, "isoseq_gene_modules_stats.txt")

isoseq_module_id <- 1
isoseq_module_genes <- names(which(membership(isoseq_communities) == isoseq_module_id))
isoseq_module_edges <- isoseq_gene_pairs[
  isoseq_gene_pairs$V1 %in% isoseq_module_genes & 
    isoseq_gene_pairs$V2 %in% isoseq_module_genes, 
]
isoseq_module_graph <- graph_from_data_frame(isoseq_module_edges, directed = FALSE)

V(isoseq_module_graph)$color <- "lightblue"
V(isoseq_module_graph)$size <- 8
V(isoseq_module_graph)$label.cex <- 0.8

set.seed(123)
plot(isoseq_module_graph,
     layout = layout_with_fr(isoseq_module_graph),
     main = paste("Module", isoseq_module_id, "Gene-Gene Interactions")
)
write_tsv(isoseq_module_edges, 
          file = paste0("isoseq_module_", isoseq_module_id, "_edges.txt"))

isoseq_target_genes <- c("OsMH_10G0001300")
isoseq_target_modules <- unique(isoseq_module_ids[names(isoseq_module_ids) %in% isoseq_target_genes])
isoseq_target_modules
isoseq_module_genes <- names(isoseq_module_ids[isoseq_module_ids %in% isoseq_target_modules])
isoseq_neighbors_all <- lapply(isoseq_target_genes, function(gene) neighbors(isoseq_g, gene)) %>% unlist() %>% unique()
isoseq_neighbor_names <- V(g)$name[neighbors_all]


# module ID and module size -----------------------------------------------
p <- bind_rows(isoseq_module_stats %>% add_column(tech='IsoSeq') %>% slice(1:10),
          drs_module_stats %>% add_column(tech='dRNA') %>% slice(1:10),
          cDNA_module_stats %>% add_column(tech='cDNA') %>% slice(1:10)) %>% 
  mutate_at(.vars = 'module_id', ~factor(.x, levels=module_id %>% unique())) %>% 
  mutate_at(.vars = 'tech', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>% 
  mutate(category_sorted=reorder_within(module_id, -gene_count, tech)) %>%
  ggplot(aes(x=category_sorted, y=gene_count, fill=tech)) + geom_bar(stat = 'identity') + facet_grid(.~tech, scales = 'free') + labs(x = "Module ID (Ranked by Size)", y = "Number of Genes") + theme(axis.text.x = element_text(angle = 90, size=4)) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none', axis.text.x = element_text(angle = 90, size=5)) + scale_fill_manual(values = c('#81B3D6','#6D60BB','#F57F73'))
p+ canvas(5,4)
ggsave(plot = p, filename = 'top_10_module_vs_size.pdf', width = 5, height = 4)


bind_rows(isoseq_module_stats %>% add_column(tech='IsoSeq') %>% add_column(x=sprintf('%03d', 1:nrow(isoseq_module_stats))),
          drs_module_stats %>% add_column(tech='dRNA') %>% add_column(x=sprintf('%03d', 1:nrow(drs_module_stats))),
          cDNA_module_stats %>% add_column(tech='cDNA') %>% add_column(x=sprintf('%03d', 1:nrow(cDNA_module_stats)))) %>% 
  group_by(tech) %>%
  mutate(cum_gene=cumsum(gene_count),
         cum_perc = cumsum(gene_count)/sum(gene_count)) %>% 
  write_tsv('module_cum_gene_count.txt')

p <- read.table('module_cum_gene_count_top10.txt', header = T) %>% 
  mutate_at(.vars = 'tech', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>%
  mutate_at(.vars = 'x', ~factor(.x, levels=x %>% unique())) %>% 
  mutate(category_sorted=reorder_within(module_id, cum_gene, tech)) %>% 
  ggplot(aes(x=category_sorted, y = cum_perc, fill=tech)) + geom_point(aes(color=tech)) + geom_line(aes(group=tech, color=tech)) + facet_grid(.~tech, scales = 'free') +  theme(axis.text.x = element_text(angle = 90, size=4)) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none', axis.text.x = element_text(angle = 90, size=5)) + scale_color_manual(values = c('#81B3D6','#6D60BB','#F57F73'))
ggsave(plot = p, filename = 'top_10_module_vs_size_cumfreq.pdf', width = 5, height = 4)