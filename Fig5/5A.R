### Figure 5
### A
### top

library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggthemes)

### Heart
### load data
load("heart.harmony_singlet_round2.RData")

### get the DEGs for Old group vs non-Old groups
get_degs <- function(seurat_harmony, group_name){
  ### set the droup order
  group_list <- levels(seurat_harmony@meta.data$origin)[which(levels(seurat_harmony@meta.data$origin) %in% c("Old", group_name))]

  ### get label for each celltype and origin
  seurat_harmony_Oi <- subset(seurat_harmony, subset = origin %in% group_list)
  seurat_harmony_Oi@meta.data$celltype_condition <- paste(seurat_harmony_Oi@meta.data$celltype, seurat_harmony_Oi@meta.data$origin, sep = "_")
  marker_condition <- data.frame()
  Idents(seurat_harmony_Oi) <- "celltype_condition"

  for (ci in as.character(unique(seurat_harmony_Oi@meta.data$celltype))){
    tmp.marker <- FindMarkers(seurat_harmony_Oi, logfc.threshold = 0, min.pct = 0.1, only.pos = F,
                              ident.1 = paste0(ci, "_", group_list[2]), ident.2 = paste0(ci, "_", group_list[1]), test.use = "wilcox")

    tmp.marker$gene <- rownames(tmp.marker)
    tmp.marker$condition <- ifelse(tmp.marker$avg_log2FC > 0, paste0(ci, "_", group_list[2]), paste0(ci, "_", group_list[1]))
    tmp.marker$cluster <- ci

    tmp.marker <- as.data.frame(tmp.marker)
    tmp.marker <- tmp.marker %>% arrange(desc(avg_log2FC))

    marker_condition <- marker_condition %>% rbind(tmp.marker)
  }

  marker_condition$sig <- ""
  marker_condition$sig[abs(marker_condition$avg_log2FC) > 0.25 & marker_condition$p_val_adj < 0.05] <- "sig"
  marker_condition$sig2 <- paste(marker_condition$cluster, marker_condition$sig, sep = "_")
  marker_condition$sig2[str_detect(marker_condition$sig2, "_$")] <- "not_sig"
  marker_condition$sig2 <- str_replace(marker_condition$sig2, "_sig", "")
  deg <- Reduce(rbind, by(marker_condition, marker_condition$cluster, function(x){x[x$sig == "sig", ]}))

  write.csv(deg, file = paste0("degs_sig_", group_list[2], "_", group_list[1], ".csv"), quote = F)
  return(deg)
}


### Aging and Drug DEGs(MET, NR, D+Q, SPD)
old_young <- get_degs(heart.harmony, "Young")
met_old <- get_degs(heart.harmony, "MET")
nr_old <- get_degs(heart.harmony, "NR")
dq_old <- get_degs(heart.harmony, "D+Q")
spd_old <- get_degs(heart.harmony, "SPD")

### prepare data for barplot
### Old vs. Young
cell <-  paste(rep(levels(heart.harmony), each = 2), rep(c("Old", "Young"), times = 2), sep = "_")
count <- ifelse(cell %in% names(table(old_young$condition)), table(old_young$condition)[match(cell, names(table(old_young$condition)))], 0)
df1 <- data.frame(cell = cell, count = count)
df1$class <- ifelse(grepl("Old", df1$cell), "Upregulated", "Downregulated")
df1$cell <- paste("Old", df1$cell, sep = "!")
df1$cell <- gsub("_[A-Za-z].*", "", df1$cell)
cell_order <- names(table(df1$cell))[order(aggregate(df1$count, by = list(df1$cell), sum)$x, decreasing = T)]
df1$cell <- factor(df1$cell, levels = cell_order)

### MET vs. Old
cell <-  paste(rep(levels(heart.harmony), each = 2), rep(c("MET", "Old"), times = 2), sep = "_")
count <- ifelse(cell %in% names(table(met_old$condition)), table(met_old$condition)[match(cell, names(table(met_old$condition)))], 0)
df2 <- data.frame(cell = cell, count = count)
df2$class <- ifelse(grepl("MET", df2$cell), "Upregulated", "Downregulated")
df2$cell <- paste("MET", df2$cell, sep = "!")
df2$cell <- gsub("_[A-Za-z].*", "", df2$cell)
cell_order <- names(table(df2$cell))[order(aggregate(df2$count, by = list(df2$cell), sum)$x, decreasing = T)]
df2$cell <- factor(df2$cell, levels = cell_order)

### NR
cell <-  paste(rep(levels(heart.harmony), each = 2), rep(c("NR", "Old"), times = 2), sep = "_")
count <- ifelse(cell %in% names(table(nr_old$condition)), table(nr_old$condition)[match(cell, names(table(nr_old$condition)))], 0)
df3 <- data.frame(cell = cell, count = count)
df3$class <- ifelse(grepl("NR", df3$cell), "Upregulated", "Downregulated")
df3$cell <- paste("NR", df3$cell, sep = "!")
df3$cell <- gsub("_[A-Za-z].*", "", df3$cell)
cell_order <- names(table(df3$cell))[order(aggregate(df3$count, by = list(df3$cell), sum)$x, decreasing = T)]
df3$cell <- factor(df3$cell, levels = cell_order)

### D+Q
cell <-  paste(rep(levels(heart.harmony), each = 2), rep(c("D+Q", "Old"), times = 2), sep = "_")
count <- ifelse(cell %in% names(table(dq_old$condition)), table(dq_old$condition)[match(cell, names(table(dq_old$condition)))], 0)
df4 <- data.frame(cell = cell, count = count)
df4$class <- ifelse(grepl("D\\+Q", df4$cell), "Upregulated", "Downregulated")
df4$cell <- paste("D+Q", df4$cell, sep = "!")
df4$cell <- gsub("_[A-Za-z].*", "", df4$cell)
cell_order <- names(table(df4$cell))[order(aggregate(df4$count, by = list(df4$cell), sum)$x, decreasing = T)]
df4$cell <- factor(df4$cell, levels = cell_order)

### SPD
cell <-  paste(rep(levels(heart.harmony), each = 2), rep(c("SPD", "Old"), times = 2), sep = "_")
count <- ifelse(cell %in% names(table(spd_old$condition)), table(spd_old$condition)[match(cell, names(table(spd_old$condition)))], 0)
df5 <- data.frame(cell = cell, count = count)
df5$class <- ifelse(grepl("SPD", df5$cell), "Upregulated", "Downregulated")
df5$cell <- paste("SPD", df5$cell, sep = "!")
df5$cell <- gsub("_[A-Za-z].*", "", df5$cell)
cell_order <- names(table(df5$cell))[order(aggregate(df5$count, by = list(df5$cell), sum)$x, decreasing = T)]
df5$cell <- factor(df5$cell, levels = cell_order)

### all
df <- rbind(df1, df2, df3, df4, df5)
p1 <- ggplot(df, aes(x = cell, y = count, fill = class)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell type", y = "DEGs number") +
  scale_fill_manual(values = c("#709ec4", "#eca3c6")) +
  theme_few() +
  theme(text = element_text(colour = "black", size = 15),
        axis.text.x = element_text(angle = 90, hjust  = 1, vjust = 0.5),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.85),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  scale_x_discrete(labels = setNames(gsub(".*!", "", df$cell), df$cell))

#设置颜色
color_heart <- c('#cb4936', '#ff0000', '#984EA3', '#1616c8', '#0D63A5',
                 '#A12568', '#F7C394', '#F499C1', '#20b2aa', '#ade87c',
                 '#EE4590', '#0d6c0d', '#6A5ACD')
cluster_order <- c("Fib", "nmSC", "PC", "SMC", "Vas_EC",
                   "Endo_EC", "LEC", "EPI", "VEN", "ATR",
                   "Mac1", "Mac2", "TC")
col_cluster <- setNames(color_heart, cluster_order)

col <- data.frame(x = 0, y = c(levels(df1$cell), levels(df2$cell), levels(df3$cell), levels(df4$cell), levels(df5$cell)), stringsAsFactors = F )
col$y <- factor(col$y, levels = col$y)
col$z <- col_cluster[match(gsub(".*!", "", col$y), names(col_cluster))]
p2 <- ggplot(col, aes(x, y, color = factor(y))) +
      geom_point(size = 5, show.legend = F) +
      scale_color_manual(values = col$z) +
      theme_classic() + coord_flip() +
      scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position = "right") +
      theme(
         plot.margin = margin(r=0),
         axis.title = element_blank(),
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank()
    )

pdf("Heart_solo_DEGs.pdf", width = 18)
p1 + p2 + plot_layout(ncol = 1, heights = c(1.3, 1))
dev.off()


### Figure 5
### A
### bottom

### get pathway name and type in heart heatmap
pathway_select <- read.xlsx("heatmap.xlsx", sheet = 1)
pathway_select_heart <- pathway_select %>% filter(origin == "Heart")
rownames(pathway_select_heart) <- pathway_select_heart$Pathway_name

### loa heatmap data
load("heart_heatmap.RData")
df_heart <- melt(df_final)
df_heart$Pathway_name <- rownames(df_final)
df_heart$pathway_type <- pathway_select_heart[df_heart$Pathway_name, "pathway_type"]
df_heart$direction <- pathway_select_heart[df_heart$Pathway_name, "direction"]
colnames(df_heart)[1:2] <- c("drug_celltype", "value")

### calculate sum of heatmap score
df_heart_mean <- df_heart %>% group_by(drug_celltype, pathway_type, direction) %>% reframe(sum_score = sum(value))
df_heart_mean$origin <- apply(df_heart_mean, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
df_heart_mean <- df_heart_mean %>% group_by(origin, pathway_type, direction) %>% reframe(Score = sum(sum_score))

### save
write.csv(df_heart_mean, file = "sankey_data_heart.csv")
