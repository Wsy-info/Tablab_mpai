### Figure 4
### EFG

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)

### load data
load("brain.harmony_singlet_round2.RData")

### load pyScenic results
load("source_target.RData")
source_target_mef2c <- source_target %>% filter(from == "Mef2c(+)")
gene_target <- source_target_mef2c$to

### get the number of down-regulated target genes between Old vs Young
Idents(brain.harmony) <- "origin"
expression_target <- as.data.frame(AverageExpression(brain.harmony, features = gene_target)$RNA)
length(which(expression_target$Young > expression_target$Old))

### Figure 4G
gene_target <- rownames(expression_target)[which(expression_target$Young > expression_target$Old)]
DotPlot(brain.harmony, features = gene_target, dot.scale = 15, idents = c("Young", "Old", "MET", "NR", "D+Q", "SPD")) +
  scale_size(range = c(2, 14)) +
  RotatedAxis() +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  labs(x = "", y = "")
ggsave("Mef2c_target_dotplot.pdf", width = 8, height = 3)


### Figure 4F
brain.harmony <- AddModuleScore(brain.harmony, features = list(gene_target), name = "Target")
### ridge plot
color_nmn <- c("#d8c91c", "#47a3bc", "#a3c8dc", "#a91f24", "#f09594",
               "#8e4b98", "#F2B71F", "#0160A6")
names(color_nmn) <- c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")

ggplot(data = brain.harmony@meta.data,
             aes(x = Target1, y = origin, fill = origin)) +
  scale_fill_manual(values = color_nmn) +
  geom_density_ridges(alpha = 0.8,
                      color = 'white',
                      rel_min_height = 0.01,
                      scale = 1.8,
                      quantile_lines = TRUE,
                      quantiles = 2
  ) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  ) +
  labs(x = "Regulon gene set score", y = "")
ggsave("Mef2c_target_gene_score_ridge_plot.pdf", width = 6, height = 4)


### GO term enrichment
entrez_ids <- bitr(gene_target, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # 选择生物过程（BP）、细胞组分（CC）或分子功能（MF）
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results <- simplify(go_results)
simplified_go_results <- simplified_go_results@result %>% arrange(pvalue)

simplified_go_results_use <- simplified_go_results %>% filter(Count >= 3)

### prepare data for network (Cytoscape_v3.10.2)
net_node <- data.frame(Node = c(gene_target, simplified_go_results_use$Description, "Mef2c(+)"),
                       Type = c(rep("gene", 13), rep("pathway", 8), "TF"))
net_edge <- c()
for(i in 1:(dim(simplified_go_results_use)[1])){
  gene_use <- strsplit(simplified_go_results_use[i, 8], split = "\\/")[[1]]
  net_edge <- rbind(net_edge, data.frame(from = gene_use, to = simplified_go_results_use[i, 2]))
}
net_edge <- rbind(net_edge, data.frame(from = "Mef2c(+)", to = gene_target))
write.table(net_edge, file = "net_edge.txt", sep = "\t", row.names = F, quote = F)
write.table(net_node, file = "net_node.txt", sep = "\t", row.names = F, quote = F)