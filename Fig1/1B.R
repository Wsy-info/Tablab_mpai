### Figure 1
### B

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
cell <- df %>% group_by(ident) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

rownames(cell) <- cell$ident
A <- cell[cluster_order, ]

### 1) UAMP plot
p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = ident)) +
      geom_point(size = 0.3, shape = 16) +
      scale_color_manual("", values = col_cluster) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text_repel(data = A, aes(label = ident),color = "black", size = 5, point.padding = 0.3)


### 2) barplot
count_df <- data.frame(Cell_Type = names(table(brain.harmony@meta.data$celltype)),
                       Cell_Count = unname(table(brain.harmony@meta.data$celltype)))[, c(1, 3)]
count_df$Cell_Type <- factor(count_df$Cell_Type, levels =  rev(cluster_order))
colnames(count_df) <- c("Cell_Type", "Cell_Count")

p2 <- ggplot(count_df, aes(x = Cell_Type, y = Cell_Count, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text(aes(label = Cell_Count), hjust = -0.05) +
      coord_flip() + ylim(0, 132000) +
      scale_fill_manual(values = rev(color_brain))

### 3) legend
df <- data.frame(x = 0, y = rev(levels(brain.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = df$y)

p3 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = rev(color_brain)) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

### save plot
pdf("brain_cell_harmony.pdf", width = 11.2)
plot_grid(plotlist = list(p1, p2, p3), ncol = 3, rel_widths = c(3, 1.3, 0.6))
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
cell <- df %>% group_by(ident) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

rownames(cell) <- cell$ident
A <- cell[cluster_order, ]

### 1) UAMP plot
p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = ident)) +
      geom_point(size = 0.3, shape = 16) +
      scale_color_manual("", values = col_cluster) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text_repel(data = A, aes(label = ident),color = "black", size = 5, point.padding = 0.3)


### 2) barplot
count_df <- data.frame(Cell_Type = names(table(heart.harmony@meta.data$celltype)),
                       Cell_Count = unname(table(heart.harmony@meta.data$celltype)))[, c(1, 3)]
count_df$Cell_Type <- factor(count_df$Cell_Type, levels =  rev(cluster_order))
colnames(count_df) <- c("Cell_Type", "Cell_Count")

p2 <- ggplot(count_df, aes(x = Cell_Type, y = Cell_Count, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text(aes(label = Cell_Count), hjust = -0.05) +
      coord_flip() + ylim(0, 132000) +
      scale_fill_manual(values = rev(color_heart))

### 3) legend
df <- data.frame(x = 0, y = rev(levels(heart.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = df$y)

p3 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = rev(color_heart)) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

### save plot
pdf("heart_cell_harmony.pdf", width = 11.2)
plot_grid(plotlist = list(p1, p2, p3), ncol = 3, rel_widths = c(3, 1.3, 0.6))
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
cell <- df %>% group_by(ident) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

rownames(cell) <- cell$ident
A <- cell[cluster_order, ]

### 1) UAMP plot
p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = ident)) +
      geom_point(size = 0.3, shape = 16) +
      scale_color_manual("", values = col_cluster) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text_repel(data = A, aes(label = ident),color = "black", size = 5, point.padding = 0.3)


### 2) barplot
count_df <- data.frame(Cell_Type = names(table(kidney.harmony@meta.data$celltype)),
                       Cell_Count = unname(table(kidney.harmony@meta.data$celltype)))[, c(1, 3)]
count_df$Cell_Type <- factor(count_df$Cell_Type, levels =  rev(cluster_order))
colnames(count_df) <- c("Cell_Type", "Cell_Count")

p2 <- ggplot(count_df, aes(x = Cell_Type, y = Cell_Count, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text(aes(label = Cell_Count), hjust = -0.05) +
      coord_flip() + ylim(0, 132000) +
      scale_fill_manual(values = rev(color_kidney))

### 3) legend
df <- data.frame(x = 0, y = rev(levels(kidney.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = df$y)

p3 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = rev(color_kidney)) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

### save plot
pdf("kidney_cell_harmony.pdf", width = 11.2)
plot_grid(plotlist = list(p1, p2, p3), ncol = 3, rel_widths = c(3, 1.3, 0.6))
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
cell <- df %>% group_by(ident) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

rownames(cell) <- cell$ident
A <- cell[cluster_order, ]

### 1) UAMP plot
p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = ident)) +
      geom_point(size = 0.3, shape = 16) +
      scale_color_manual("", values = col_cluster) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text_repel(data = A, aes(label = ident),color = "black", size = 5, point.padding = 0.3)


### 2) barplot
count_df <- data.frame(Cell_Type = names(table(liver.harmony@meta.data$celltype)),
                       Cell_Count = unname(table(liver.harmony@meta.data$celltype)))[, c(1, 3)]
count_df$Cell_Type <- factor(count_df$Cell_Type, levels =  rev(cluster_order))
colnames(count_df) <- c("Cell_Type", "Cell_Count")

p2 <- ggplot(count_df, aes(x = Cell_Type, y = Cell_Count, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text(aes(label = Cell_Count), hjust = -0.05) +
      coord_flip() + ylim(0, 132000) +
      scale_fill_manual(values = rev(color_liver))

### 3) legend
df <- data.frame(x = 0, y = rev(levels(liver.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = df$y)

p3 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = rev(color_liver)) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

### save plot
pdf("liver_cell_harmony.pdf", width = 11.2)
plot_grid(plotlist = list(p1, p2, p3), ncol = 3, rel_widths = c(3, 1.3, 0.6))
dev.off()