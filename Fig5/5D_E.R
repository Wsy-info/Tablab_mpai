### Figure 5
### DE

library(Seurat)
library(sceasy)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

### Brain
load("heart.harmony_singlet_round2.RData")

### use aging-related DEGs
### only retain Young and Old groups
heart_O_Y <- subset(heart.harmony, subset = origin %in% c("Young", "Old"))
heart_O_Y@meta.data$origin <- as.character(heart_O_Y@meta.data$origin)
heart_O_Y@meta.data$origin <- factor(heart_O_Y@meta.data$origin, levels = c("Young", "Old"))
Idents(heart_O_Y) <- "origin"

### get DEGs for each celltype
all_markers_O_Y <- c()
for(i in levels(heart_O_Y@meta.data$celltype)){
  heart_O_Y_i <- subset(heart_O_Y, subset = celltype == i)
  Idents(heart_O_Y_i) <- "origin"
  markers_O_Y <- FindAllMarkers(heart_O_Y_i)
  markers_O_Y$celltype <- i
  all_markers_O_Y <- rbind(all_markers_O_Y, markers_O_Y)
}

heart_O_Y_aging <- subset(heart_O_Y, features = unique(all_markers_O_Y$gene))
### write for h5ad
convertFormat(heart_O_Y_aging, from="seurat", to="anndata", outFile='heart_O_Y_aging.h5ad')

heart_aging <- subset(heart.harmony, features = unique(all_markers_O_Y$gene))
### write for h5ad
convertFormat(heart_aging, from="seurat", to="anndata", outFile='heart_aging.h5ad')


### Fig5/5D_E.py get predict results
### set color for origin
color_nmn <- c("#47a3bc", "#d8c91c", "#a3c8dc", "#a91f24", "#f09594",
               "#8e4b98", "#F1B412", "#0066B2")
names(color_nmn) <- c("Old", "Young", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")

### prepare data for plot
heart_n_o_y <- subset(heart.harmony, subset = origin %in% c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))
result_heart <- read.csv("prediction_heart_aging_bys.csv")
heart_n_o_y@meta.data$log_label <- result_heart$y_pred
x <- sort((as.matrix(table(result_heart$y_pred, heart_n_o_y@meta.data$origin))[1, 3:8]) / colSums(as.matrix(table(result_heart$y_pred, heart_n_o_y@meta.data$origin))[, 3:8]))
x <- data.frame(origin = names(x), Ratio = x)
x$origin <- factor(x$origin, levels = c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

### plot
celltype_pre_age <- heart_n_o_y@meta.data %>% group_by(celltype) %>% reframe(Ratio = (as.matrix(table(log_label, origin))[1, 3:8]) / colSums(as.matrix(table(log_label, origin))[, 3:8])) %>% as.data.frame()
celltype_pre_age$origin <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
celltype_pre_age <- celltype_pre_age %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_pre_age) <- celltype_pre_age[, 1]
celltype_pre_age <- celltype_pre_age[, -1]

### for Figure 5E
pheatmap(celltype_pre_age,
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         scale = "none",
         breaks = seq(from = 0, to = 1, length.out = 100),
         filename = "heart_ratio_by_cell_type_heatmap_all.pdf",
         width = 5, height = 3
)

### for Figure S13A
pheatmap(celltype_pre_age[c(1, 2, 5), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "heart_ratio_by_cell_type_heatmap_MN.pdf",
         width = 5, height = 2
)

### for Figure S13D
pheatmap(celltype_pre_age[c(1, 4, 6), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "heart_ratio_by_cell_type_heatmap_MS.pdf",
         width = 5, height = 2
)

### barplot
### for Figure 5D
ggplot(data = x %>% filter(origin %in% c("MET", "NR", "D+Q", "SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Heart")
ggsave("heart_ratio.pdf", width = 4, height = 4)

### For Figure S13A
ggplot(data = x %>% filter(origin %in% c("MET", "NR", "MET+NR")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Heart")
ggsave("heart_ratio_MN.pdf", width = 4, height = 4)

### For Figure S13D
ggplot(data = x %>% filter(origin %in% c("MET", "SPD", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Heart")
ggsave("heart_ratio_MS.pdf", width = 4, height = 4)