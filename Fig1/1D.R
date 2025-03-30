### Figure 1
### D

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
### set color for origin
color_origin <- c("#d8c91c", "#47a3bc", "#a3c8dc", "#a91f24", "#f09594",
                  "#8e4b98", "#fdbd10", "#0066b2")
origin_order <- c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
col_origin <- setNames(color_origin, origin_order)

### get the proportion of origin for each celltype
brain.harmony.main <- subset(brain.harmony, subset = origin %in% levels(brain.harmony@meta.data$origin)[1:6])
data_df <- brain.harmony.main@meta.data

cellratio_df <- as.data.frame(prop.table(table(data_df$origin, data_df$celltype), margin = 2))
colnames(cellratio_df) <- c("Origin", "celltype", "Ratio")
cellratio_df$Origin <- as.character(cellratio_df$Origin)
cellratio_df$Origin <- factor(cellratio_df$Origin, levels = levels(brain.harmony@meta.data$origin)[1:6])


### save
pdf("brain_proportion.pdf")
ggplot(cellratio_df, aes(x = celltype, y = Ratio, fill = Origin)) +
  geom_col(width = 0.5, color = NA) +
  scale_fill_manual(values = col_origin) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL, y = 'Cell proportion', x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "none"
  )
dev.off()





### 2.Heart
### load data
load("heart.harmony_singlet_round2.RData")

### get the proportion of origin for each celltype
heart.harmony.main <- subset(heart.harmony, subset = origin %in% levels(heart.harmony@meta.data$origin)[1:6])
data_df <- heart.harmony.main@meta.data

cellratio_df <- as.data.frame(prop.table(table(data_df$origin, data_df$celltype), margin = 2))
colnames(cellratio_df) <- c("Origin", "celltype", "Ratio")
cellratio_df$Origin <- as.character(cellratio_df$Origin)
cellratio_df$Origin <- factor(cellratio_df$Origin, levels = levels(heart.harmony@meta.data$origin)[1:6])


### save
pdf("heart_proportion.pdf")
ggplot(cellratio_df, aes(x = celltype, y = Ratio, fill = Origin)) +
  geom_col(width = 0.5, color = NA) +
  scale_fill_manual(values = col_origin) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL, y = 'Cell proportion', x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "none"
  )
dev.off()





### 3.Kidney
### load data
load("kidney.harmony_singlet_round2.RData")

### get the proportion of origin for each celltype
kidney.harmony.main <- subset(kidney.harmony, subset = origin %in% levels(kidney.harmony@meta.data$origin)[1:6])
data_df <- kidney.harmony.main@meta.data

cellratio_df <- as.data.frame(prop.table(table(data_df$origin, data_df$celltype), margin = 2))
colnames(cellratio_df) <- c("Origin", "celltype", "Ratio")
cellratio_df$Origin <- as.character(cellratio_df$Origin)
cellratio_df$Origin <- factor(cellratio_df$Origin, levels = levels(kidney.harmony@meta.data$origin)[1:6])


### save
pdf("kidney_proportion.pdf")
ggplot(cellratio_df, aes(x = celltype, y = Ratio, fill = Origin)) +
  geom_col(width = 0.5, color = NA) +
  scale_fill_manual(values = col_origin) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL, y = 'Cell proportion', x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "none"
  )
dev.off()





### 4.Liver
### load data
load("liver.harmony_singlet_round2.RData")

### get the proportion of origin for each celltype
liver.harmony.main <- subset(liver.harmony, subset = origin %in% levels(liver.harmony@meta.data$origin)[1:6])
data_df <- liver.harmony.main@meta.data

cellratio_df <- as.data.frame(prop.table(table(data_df$origin, data_df$celltype), margin = 2))
colnames(cellratio_df) <- c("Origin", "celltype", "Ratio")
cellratio_df$Origin <- as.character(cellratio_df$Origin)
cellratio_df$Origin <- factor(cellratio_df$Origin, levels = levels(liver.harmony@meta.data$origin)[1:6])


### save
pdf("liver_proportion.pdf")
ggplot(cellratio_df, aes(x = celltype, y = Ratio, fill = Origin)) +
  geom_col(width = 0.5, color = NA) +
  scale_fill_manual(values = col_origin) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL, y = 'Cell proportion', x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "none"
  )
dev.off()