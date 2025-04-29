### Figure 7
### E

library(dplyr)
library(Seurat)

### load data
load("liver.harmony_singlet_round2.RData")

### DEGs for celltype compare in different groups
all_markers <- c()
for(i in c("Young", "MET", "NR", "D+Q", "SPD")){
  liver.harmony_i <- subset(liver.harmony, subset = origin %in% c(i, "Old"))
  Idents(liver.harmony_i) <- "origin"
  for(j in levels(liver.harmony_i@meta.data$celltype)){
    liver.harmony_i_j <- subset(liver.harmony_i, subset = celltype == j)
    Idents(liver.harmony_i_j) <- "origin"
    markers_i_j <- FindAllMarkers(liver.harmony_i_j, only.pos = T)
    markers_i_j$celltype <- paste0(i, "_", j)

    all_markers <- rbind(all_markers, markers_i_j)
  }
}

save(all_markers, file = "liver_DEGs_celltype.RData")

### KEGG_gene
load("KEGG_mmu_gene_list.RData")
metabolic_related_gene <- unique(as.character(unlist(kegg_mmu_list)))

### calculate the proportion of DEGs in metabolic pathways
### downregulated DEGs
down_ratio <- c()
for(i in c("Young", "MET", "NR", "D+Q", "SPD")){
  all_markers_i <- all_markers %>% filter(cluster == i)
  ratio_i <- all_markers_i %>% group_by(celltype) %>% reframe(Ratio = (length(intersect(gene, metabolic_related_gene)) / length(gene))) %>% as.data.frame()
  ratio_i$group <- i
  down_ratio <- rbind(down_ratio, ratio_i)
}


### upregulated DEGs
all_markers_old <- all_markers %>% filter(cluster == "Old")
all_markers_old$origin <- apply(all_markers_old, 1, function(x){
  strsplit(x[8], split = "_")[[1]][1]
})

up_ratio <- c()
for(i in c("Young", "MET", "NR", "D+Q", "SPD")){
  all_markers_i <- all_markers_old %>% filter(origin == i)
  ratio_i <- all_markers_i %>% group_by(celltype) %>% reframe(Ratio = (length(intersect(gene, metabolic_related_gene)) / length(gene))) %>% as.data.frame()
  ratio_i$group <- i
  up_ratio <- rbind(up_ratio, ratio_i)
}

### prepare data for plot
all_ratio <- rbind(data.frame(up_ratio %>% filter(group != "Young"), type = "down"),
                   data.frame(up_ratio %>% filter(group == "Young"), type = "up"),
                   data.frame(down_ratio %>% filter(group != "Young"), type = "up"),
                   data.frame(down_ratio %>% filter(group == "Young"), type = "down")
)
all_ratio$group <- factor(all_ratio$group, levels = c("Young", "MET", "NR", "D+Q", "SPD"),
                          labels = c("Old/Young", "MET/Old", "NR/Old", "D+Q/Old", "SPD/Old"))
all_ratio$celltype <- apply(all_ratio, 1, function (x){
  strsplit(x[1], split = "_")[[1]][2]
})
all_ratio$celltype <- factor(all_ratio$celltype, levels = levels(liver.harmony@meta.data$celltype))

color_nmn <- c("#47a3bc", "#a3c8dc", "#a91f24", "#f09594",
               "#8e4b98")

all_ratio$group_type <- paste(all_ratio$group, all_ratio$type)
all_ratio$group_type <- factor(all_ratio$group_type, levels = c("Old/Young up", "Old/Young down",
                                                                "MET/Old up", "MET/Old down",
                                                                "NR/Old up", "NR/Old down",
                                                                "D+Q/Old up", "D+Q/Old down",
                                                                "SPD/Old up", "SPD/Old down"
))

### set celltype and color
color_liver <- c("#3CB371", "#CBDF9AFF", "#4682B4", "#483D8B", "#008B8B",
                 "#FF8C00", "#FD6467")
cluster_order <- c("HSC", "EC", "Cho", "Kup", "PC-Hep",
                   "MZ-Hep", "PP-Hep")
col_cluster <- setNames(color_liver, cluster_order)

### plot
ggplot(all_ratio, aes(x = group_type, y = Ratio, fill = type)) +
  geom_boxplot(aes(color = group),
               position = position_dodge(0.9),
               outlier.color = NA,
               width = 0.75, alpha = 0.3
  ) +
  geom_jitter(size = 2, position = position_jitter(0.3),
              aes(color = celltype)) +
  scale_color_manual(values = col_cluster) +
  theme_classic() +
  scale_fill_manual(values = c("lightblue", "pink")) +
  labs(x = "", y = "Ratio") +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "none"
  )

ggsave("metabolic_gene_ratio.pdf", width = 6, height = 6)