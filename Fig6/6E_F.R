### Figure 6
### EF

library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggpubr)


### Kidney
### load data
load("kidney.harmony_singlet_round2.RData")

### get the DEGs for Old group vs non-Old groups
get_degs <- function(seurat_harmony, group_name){
  ### set the droup order
  group_list <- levels(seurat_harmony@meta.data$origin)[which(levels(seurat_harmony@meta.data$origin) %in% c("Old", group_name))]

  ### get label for each celltype and origin
  seurat_harmony_Oi <- subset(seurat_harmony, subset = origin %in% group_list)
  seurat_harmony_Oi@meta.data$celltype_condition <- paste(seurat_harmony_Oi@meta.data$celltype, seurat_harmony_Oi@meta.data$origin, sep = "_")
  marker_condition <- data.frame()
  Idents(seurat_harmony_Oi) <- "celltype_condition"

  for ( ci in as.character(unique(seurat_harmony_Oi@meta.data$celltype))){
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


### get all degs
deg1 <- get_degs(kidney.harmony, "Young")
deg2 <- get_degs(kidney.harmony, "MET")
deg3 <- get_degs(kidney.harmony, "NR")
deg4 <- get_degs(kidney.harmony, "D+Q")
deg5 <- get_degs(kidney.harmony, "SPD")
deg6 <- get_degs(kidney.harmony, "MET+NR")
deg7 <- get_degs(kidney.harmony, "MET+SPD")
cluster_order <- c("EC", "Fib", "SMC", "Podo", "PT-S1",
                   "PT-S2", "PT-S3", "LOH", "DCT1", "DCT2",
                   "CNT", "CD-IC", "CD-PC", "CD-trans", "Mac1")


### get the Rev-aging genes
list1 <- list()
for(i in 1:length(cluster_order)){
  gene1 <- deg1[deg1$condition == paste(cluster_order[i], "Young", sep = "_"), ]$gene
  gene2 <- deg2[deg2$condition == paste(cluster_order[i], "MET", sep = "_"), ]$gene
  gene3 <- deg3[deg3$condition == paste(cluster_order[i], "NR", sep = "_"), ]$gene
  gene4 <- deg4[deg4$condition == paste(cluster_order[i], "D+Q", sep = "_"), ]$gene
  gene5 <- deg5[deg5$condition == paste(cluster_order[i], "SPD", sep = "_"), ]$gene
  gene6 <- deg6[deg6$condition == paste(cluster_order[i], "MET+NR", sep = "_"), ]$gene
  gene7 <- deg7[deg7$condition == paste(cluster_order[i], "MRT+SPD", sep = "_"), ]$gene
  list1[[i]] <- list("Young_MET" = intersect(gene1, gene2),
                     "Young_NR" = intersect(gene1, gene3),
                     "Young_D+Q" = intersect(gene1, gene4),
                     "Young_SPD" = intersect(gene1, gene5),
                     "Young_MET+NR" = intersect(gene1, gene6),
                     "Young_MET+SPD" = intersect(gene1, gene7))
}
names(list1) <- cluster_order

list2 <- list()
for(i in 1:length(cluster_order)){
  gene1 <- deg1[deg1$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene2 <- deg2[deg2$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene3 <- deg3[deg3$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene4 <- deg4[deg4$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene5 <- deg5[deg5$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene6 <- deg6[deg6$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene7 <- deg7[deg7$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  list2[[i]] <- list("Old_Old(MET)" = intersect(gene1, gene2),
                     "Old_Old(NR)" = intersect(gene1, gene3),
                     "Old_Old(D+Q)" = intersect(gene1, gene4),
                     "Old_Old(SPD)" = intersect(gene1, gene5),
                     "Old_Old(MET+NR)" = intersect(gene1, gene6),
                     "Old_Old(MET+SPD)" = intersect(gene1, gene7))
}
names(list2) <- cluster_order

### prepare data for barplot
met <- list()
for(i in 1:length(list1)){
 df <- length(list1[[i]][[1]]) + length(list2[[i]][[1]])
 met[[i]] <- df
}

nr <- list()
for(i in 1:length(list1)){
 df <- length(list1[[i]][[2]]) + length(list2[[i]][[2]])
 nr[[i]] <- df
}

dq <- list()
for(i in 1:length(list1)){
 df <- length(list1[[i]][[3]]) + length(list2[[i]][[3]])
 dq[[i]] <- df
}

spd <- list()
for(i in 1:length(list1)){
 df <- length(list1[[i]][[4]]) + length(list2[[i]][[4]])
 spd[[i]] <- df
}

mn <- list()
for(i in 1:length(list1)){
 df <- length(list1[[i]][[5]]) + length(list2[[i]][[5]])
 mn[[i]] <- df
}

ms <- list()
for(i in 1:length(list1)){
 df <- length(list1[[i]][[6]]) + length(list2[[i]][[6]])
 ms[[i]] <- df
}

### for Figure 6E: Rev-aging genes
### get the up-reglated genes
met <- list()
for(i in 1:length(list1)){
 df <- list1[[i]][[1]]
 met[[i]] <- df
}

nr <- list()
for(i in 1:length(list1)){
 df <- list1[[i]][[2]]
 nr[[i]] <- df
}

dq <- list()
for(i in 1:length(list1)){
 df <- list1[[i]][[3]]
 dq[[i]] <- df
}

spd <- list()
for(i in 1:length(list1)){
 df <- list1[[i]][[4]]
 spd[[i]] <- df
}

mn <- list()
for(i in 1:length(list1)){
 df <- list1[[i]][[5]]
 mn[[i]] <- df
}

ms <- list()
for(i in 1:length(list1)){
 df <- list1[[i]][[6]]
 ms[[i]] <- df
}

df_up <- cbind(length(unique(unlist(met))), length(unique(unlist(nr))), length(unique(unlist(dq))), length(unique(unlist(spd))))
colnames(df_up) <-  c("MET", "NR", "D+Q", "SPD")
rownames(df_up) <- "Up-regulated"

### get the down-reglated genes
met <- list()
for(i in 1:length(list1)){
 df <- list2[[i]][[1]]
 met[[i]] <- df
}

nr <- list()
for(i in 1:length(list1)){
 df <- list2[[i]][[2]]
 nr[[i]] <- df
}

dq <- list()
for(i in 1:length(list1)){
 df <- list2[[i]][[3]]
 dq[[i]] <- df
}

spd <- list()
for(i in 1:length(list1)){
 df <- list2[[i]][[4]]
 spd[[i]] <- df
}

mn <- list()
for(i in 1:length(list1)){
 df <- list2[[i]][[5]]
 mn[[i]] <- df
}

ms <- list()
for(i in 1:length(list1)){
 df <- list2[[i]][[6]]
 ms[[i]] <- df
}

df_down <- cbind(length(unique(unlist(met))), length(unique(unlist(nr))), length(unique(unlist(dq))), length(unique(unlist(spd))))
colnames(df_down) <-  c("MET", "NR", "D+Q", "SPD")
rownames(df_down) <- "Down-regulated"

count <- rbind(df_up, df_down)
rownames(count) <- c("Up-regulated", "Down-regulated")
colnames(count) <- c("MET", "NR", "D+Q", "SPD")

### prepare data for Rose diagram
df <- data.frame(class = rep(colnames(count), times = colSums(count)),
                 count = c(rep("Up-regulated", 33), rep("Down-regulated", 181), rep("Up-regulated", 39), rep("Down-regulated", 168),
                           rep("Up-regulated", 8), rep("Down-regulated", 258), rep("Up-regulated", 20), rep("Down-regulated", 191)))



df_down <- cbind(length(unique(unlist(met))), length(unique(unlist(nr))), length(unique(unlist(dq))), length(unique(unlist(spd))))
colnames(df_down) <-  c("MET", "NR", "D+Q", "SPD")
rownames(df_down) <- "Up-regulated"

count <- rbind(df_up, df_down)
rownames(count) <- c("Up-regulated", "Down-regulated")
colnames(count) <- c("MET", "NR", "D+Q", "SPD")

df <- data.frame(class = rep(colnames(count), times = colSums(count)),
                 count = c(rep("Up-regulated", 33), rep("Down-regulated", 181), rep("Up-regulated", 39), rep("Down-regulated", 168),
                           rep("Up-regulated", 8), rep("Down-regulated", 258), rep("Up-regulated", 20), rep("Down-regulated", 191)))

pdf("Kidney_Rose_plot_revaging.pdf")
ggplot(df, aes(class, fill = count)) +
  geom_bar(alpha = 0.9) +
  coord_polar() +
  xlab("") + scale_fill_manual(values = c("#709ec4", "#eca3c6"), labels = c("Down-regulated", "Up-regulated")) +
  ylab("") + labs(title = "Rev-aging DEGs (Kidney)") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        text = element_text(size = 15),
        panel.grid = element_line(size = 0.5, color = "gray60"),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()





### for Figure 6F: Pro-aging genes
### get the Pro-aging genes
list1 <- list()
for(i in 1:length(cluster_order)){
  gene1 <- deg1[deg1$condition == paste(cluster_order[i], "Young", sep = "_"), ]$gene
  gene2 <- deg2[deg2$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene3 <- deg3[deg3$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene4 <- deg4[deg4$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene5 <- deg5[deg5$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene6 <- deg6[deg6$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene7 <- deg7[deg7$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  list1[[i]] <- list("Young_Old(MET)" = intersect(gene1, gene2),
                     "Young_Old(NR)" = intersect(gene1, gene3),
                     "Young_Old(D+Q)" = intersect(gene1, gene4),
                     "Young_Old(SPD)" = intersect(gene1, gene5),
                     "Young_Old(MET+NR)" = intersect(gene1, gene6),
                     "Young_Old(MET+SPD)" = intersect(gene1, gene7))
}

names(list1) <- cluster_order

list2 <- list()
for(i in 1:length(cluster_order)){
  gene1 <- deg1[deg1$condition == paste(cluster_order[i], "Old", sep = "_"), ]$gene
  gene2 <- deg2[deg2$condition == paste(cluster_order[i], "MET", sep = "_"), ]$gene
  gene3 <- deg3[deg3$condition == paste(cluster_order[i], "NR", sep = "_"), ]$gene
  gene4 <- deg4[deg4$condition == paste(cluster_order[i], "D+Q", sep = "_"), ]$gene
  gene5 <- deg5[deg5$condition == paste(cluster_order[i], "SPD", sep = "_"), ]$gene
  gene6 <- deg6[deg6$condition == paste(cluster_order[i], "MET+NR", sep = "_"), ]$gene
  gene7 <- deg7[deg7$condition == paste(cluster_order[i], "MET+SPD", sep = "_"), ]$gene
  list2[[i]] <- list("Old_MET" = intersect(gene1, gene2),
                     "Old_NR" = intersect(gene1, gene3),
                     "Old_D+Q" = intersect(gene1, gene4),
                     "Old_SPD" = intersect(gene1, gene5),
                     "Old_MET+NR" = intersect(gene1, gene6),
                     "Old_MET+SPD" = intersect(gene1, gene7))
}
names(list2) <- cluster_order

### get the up-reglated genes
met <- list()
for(i in 1: length(list1)){
 df <- list1[[i]][[1]]
 met[[i]] <- df
}

nr <- list()
for(i in 1: length(list1)){
 df <- list1[[i]][[2]]
 nr[[i]] <- df
}

dq <- list()
for(i in 1: length(list1)){
 df <- list1[[i]][[3]]
 dq[[i]] <- df
}

spd <- list()
for(i in 1: length(list1)){
 df <- list1[[i]][[4]]
 spd[[i]] <- df
}

mn <- list()
for(i in 1: length(list1)){
 df <- list1[[i]][[5]]
 mn[[i]] <- df
}

ms <- list()
for(i in 1: length(list1)){
 df <- list1[[i]][[6]]
 ms[[i]] <- df
}

df_up <- cbind(length(unique(unlist(met))), length(unique(unlist(nr))), length(unique(unlist(dq))), length(unique(unlist(spd))))
colnames(df_up) <-  c("MET", "NR", "D+Q", "SPD")
rownames(df_up) <- "Up-regulated"


### get the down-reglated genes
met <- list()
for(i in 1: length(list1)){
 df <- list2[[i]][[1]]
 met[[i]] <- df
}

nr <- list()
for(i in 1: length(list1)){
 df <- list2[[i]][[2]]
 nr[[i]] <- df
}

dq <- list()
for(i in 1: length(list1)){
 df <- list2[[i]][[3]]
 dq[[i]] <- df
}

spd <- list()
for(i in 1: length(list1)){
 df <- list2[[i]][[4]]
 spd[[i]] <- df
}

mn <- list()
for(i in 1: length(list1)){
 df <- list2[[i]][[5]]
 mn[[i]] <- df
}

ms <- list()
for(i in 1: length(list1)){
 df <- list2[[i]][[6]]
 ms[[i]] <- df
}

df_down <- cbind(length(unique(unlist(met))), length(unique(unlist(nr))), length(unique(unlist(dq))), length(unique(unlist(spd))))
colnames(df_down) <-  c("MET", "NR", "D+Q", "SPD")
rownames(df_down) <- "Down-regulated"

count <- rbind(df_up, df_down)
rownames(count) <- c("Up-regulated", "Down-regulated")
colnames(count) <- c("MET", "NR", "D+Q", "SPD")

### prepare data for Rose diagram
df <- data.frame(class = rep(colnames(count), times = colSums(count)),
                 count = c(rep("Up-regulated", 3), rep("Down-regulated", 35), rep("Up-regulated", 48), rep("Down-regulated", 105),
                           rep("Up-regulated", 6), rep("Down-regulated", 11), rep("Up-regulated", 24), rep("Down-regulated", 7)))

pdf("Kidney_Rose_plot_proaging.pdf")
ggplot(df, aes(class, fill = count)) +
  geom_bar(alpha = 0.9) +
  coord_polar() +
  xlab("") + scale_fill_manual(values = c("#709ec4", "#eca3c6"), labels = c("Down-regulated", "Up-regulated")) +
  ylab("") + labs(title = "Pro-aging DEGs (Brain)") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        text = element_text(size = 15),
        panel.grid = element_line(size = 0.5, color = "gray60"),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()