### Figure 1
### C

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(cowplot)
library(viridis)
library(ggthemes)

### 1.Brain
### load data
load("brain.harmony_singlet_round2.RData")
### set color for celltype
color_brain <- c("#0073C2FF", "#EEC000FF", "#CD534CFF", "#00CC33FF", "#984EA3",
                 "#8F7700FF", "#FCCDE5", "#C46DA0FF", "#F0027F", "#374E55FF",
                 "#FD6467", "#0B775E", "#FF59FF", "#FDC086", "#5B1A18",
                 "#69C8ECFF", "#3da8e0", "#CBDF9AFF", "#802268FF")
cluster_order <- c("OLG", "MG", "ASC", "OPC", "PC",
                   "InN_1", "InN_2", "InN_3", "InN_4", "InN_5",
                   "ExN_1", "ExN_2", "DoN", "D1-MSN", "D2-MSN",
                   "VLMC", "ABC", "EC", "CPC")
col_cluster <- setNames(color_brain, cluster_order)

### UMAP plot
### get the umap data
df <- DimPlot(brain.harmony, reduction = "umap")$data
df$origin <- brain.harmony@meta.data$origin

pdf("brain_origin_harmony.pdf", width = 8.5)
ggplot() +
  geom_point(data = df,
             mapping = aes(UMAP_1, UMAP_2, color = ident),
             size = 0.3, stroke = 0.5) +
  facet_wrap(~origin, ncol = 3) +
  scale_color_manual("", values = col_cluster) +
  theme_classic() +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
dev.off()





### 2.Heart
### load data
load("heart.harmony_singlet_round2.RData")
### set color for celltype
color_heart <- c('#cb4936', '#ff0000', '#984EA3', '#1616c8', '#0D63A5',
                 '#A12568', '#F7C394', '#F499C1', "#20b2aa", '#ade87c',
                 '#EE4590', '#0d6c0d', '#6A5ACD')
cluster_order <-  c("Fib", "nmSC", "PC", "SMC", "Vas_EC",
                    "Endo_EC", "LEC", "EPI", "VEN", "ATR",
                    "Mac1", "Mac2", "TC")
col_cluster <- setNames(color_heart, cluster_order)

### UMAP plot
### get the umap data
df <- DimPlot(heart.harmony, reduction = "umap")$data
df$origin <- heart.harmony@meta.data$origin

pdf("heart_origin_harmony.pdf", width = 8.5)
ggplot() +
  geom_point(data = df,
             mapping = aes(UMAP_1, UMAP_2, color = ident),
             size = 0.3, stroke = 0.5) +
  facet_wrap(~origin, ncol = 3) +
  scale_color_manual("", values = col_cluster) +
  theme_classic() +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
dev.off()





### 3.Kidney
### load data
load("kidney.harmony_singlet_round2.RData")
color_kidney <- c('#CBDF9AFF', '#cb4936', '#1616c8', '#007ABA', '#8B58A4',
                  '#DF75AE', '#00B7CA', '#A8C7E9', '#E91E25', '#F37121',
                  '#FBB36E', '#F58D93', '#00A163', '#8CCA7C', '#EE4590')
cluster_order <- c("EC", "Fib", "SMC", "Podo", "PT-S1",
                   "PT-S2", "PT-S3", "LOH", "DCT1", "DCT2",
                   "CNT", "CD-IC", "CD-PC", "CD-trans", "Mac1")
col_cluster <- setNames(color_kidney, cluster_order)

### UMAP plot
### get the umap data
df <- DimPlot(kidney.harmony, reduction = "umap")$data
df$origin <- kidney.harmony@meta.data$origin

pdf("kidney_origin_harmony.pdf", width = 8.5)
ggplot() +
  geom_point(data = df,
             mapping = aes(UMAP_1, UMAP_2, color = ident),
             size = 0.3, stroke = 0.5) +
  facet_wrap(~origin, ncol = 3) +
  scale_color_manual("", values = col_cluster) +
  theme_classic() +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
dev.off()





### 4.Liver
### load data
load("liver.harmony_singlet_round2.RData")
color_liver <- c("#3CB371", "#CBDF9AFF", "#4682B4", "#483D8B", "#008B8B",
                 "#FF8C00", "#FD6467")
cluster_order <- c("HSC", "EC", "Cho", "Kup", "PC-Hep",
                   "MZ-Hep", "PP-Hep")
col_cluster <- setNames(color_liver, cluster_order)

### UMAP plot
### get the umap data
df <- DimPlot(liver.harmony, reduction = "umap")$data
df$origin <- liver.harmony@meta.data$origin

pdf("liver_origin_harmony.pdf", width = 8.5)
ggplot() +
  geom_point(data = df,
             mapping = aes(UMAP_1, UMAP_2, color = ident),
             size = 0.3, stroke = 0.5) +
  facet_wrap(~origin, ncol = 3) +
  scale_color_manual("", values = col_cluster) +
  theme_classic() +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
dev.off()