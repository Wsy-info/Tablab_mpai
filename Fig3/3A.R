### Figure 3
### A

library(Seurat)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(stringr)

### Brain
### load data
load("brain.harmony_singlet_round2.RData")

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


### get all degs
deg1 <- get_degs(brain.harmony, "Young")
deg2 <- get_degs(brain.harmony, "MET")
deg3 <- get_degs(brain.harmony, "NR")
deg4 <- get_degs(brain.harmony, "D+Q")
deg5 <- get_degs(brain.harmony, "SPD")
deg6 <- get_degs(brain.harmony, "MET+NR")
deg7 <- get_degs(brain.harmony, "MET+SPD")
cluster_order <- c("OLG", "MG", "ASC", "OPC", "PC",
                   "InN_1", "InN_2", "InN_3", "InN_4", "InN_5",
                   "ExN_1", "ExN_2", "DoN", "D1-MSN", "D2-MSN",
                   "VLMC", "ABC", "EC", "CPC")


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

df_final <- cbind(as.numeric(met), as.numeric(nr), as.numeric(dq), as.numeric(spd), as.numeric(mn), as.numeric(ms))
colnames(df_final) <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
rownames(df_final) <- cluster_order
data1 <- melt(df_final)
data1$group <- "Rev-aging"


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
capture.output(list1, file = "Brain_Pro-aging_genes_1.csv")

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
capture.output(list2, file = "Brain_Pro-aging_genes_2.csv")

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

df_final <- cbind(as.numeric(met), as.numeric(nr), as.numeric(dq), as.numeric(spd), as.numeric(mn), as.numeric(ms))
colnames(df_final) <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
rownames(df_final) <- cluster_order
data2 <- melt(df_final)
data2$group <- "Pro-aging"
data <- rbind(data1, data2)


### for MET, NR, D+Q, SPD
data1 <- data[data$Var2 %in% c("MET", "NR", "D+Q", "SPD"), ]

### sort the number of Rev-aging genes
cell <- data1[1:19, ]$Var1
cell_order <- cell[order(aggregate(data1[data1$group == "Rev-aging", ]$value, by = list(data1[data1$group == "Rev-aging", ]$Var1), sum)$x, decreasing = T)]
data1$Var1 <- factor(data1$Var1, levels = cell_order)

pdf("brain_solo_Rev_Pro.pdf")
ggplot(data1, aes(x = Var1, y = ifelse(group == "Rev-aging", value, -value), fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Change", values = c("#a3c8dc", "#a91f24", "#f09594", "#8e4b98")) +
  scale_y_continuous(
    labels = abs,
    expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "", y = "Gene number", title = "The number of Rev-aging and Pro-aging DEGs in Brain") +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(text = element_text(colour = "black", size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        panel.background = element_rect(size = 1, color = "black"))
dev.off()



### for Figure 8B
### for MET+NR
data1 <- data[data$Var2 %in% c("MET", "NR", "MET+NR"), ]
### sort the number of Rev-aging genes
cell <- data1[1:19, ]$Var1
cell_order <- cell[order(aggregate(data1[data1$group == "Rev-aging", ]$value, by = list(data1[data1$group == "Rev-aging", ]$Var1), sum)$x, decreasing = T)]
data1$Var1 <- factor(data1$Var1, levels = cell_order)

pdf("brain_MN_Rev_Pro.pdf")
ggplot(data1, aes(x = Var1, y = ifelse(group == "Rev-aging", value, -value), fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Change", values = c("#a3c8dc", "#a91f24", "#fdbd10")) +
  scale_y_continuous(
    labels = abs,
    expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "", y = "Gene number", title = "The number of Rev-aging and Pro-aging DEGs in Brain") +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(text = element_text(colour = "black", size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        panel.background = element_rect(size = 1, color = "black"))
dev.off()


### for Figure 8C
### for MET+SPD
data1 <- data[data$Var2 %in% c("MET", "SPD", "MET+SPD"), ]
### sort the number of Rev-aging genes
cell <- data1[1:19, ]$Var1
cell_order <- cell[order(aggregate(data1[data1$group == "Rev-aging", ]$value, by = list(data1[data1$group == "Rev-aging", ]$Var1), sum)$x, decreasing = T)]
data1$Var1 <- factor(data1$Var1, levels = cell_order)

pdf("brain_MS_Rev_Pro.pdf")
ggplot(data1, aes(x = Var1, y = ifelse(group == "Rev-aging", value, -value), fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Change", values = c("#a3c8dc", "#8e4b98", "#0066b2")) +
  scale_y_continuous(
    labels = abs,
    expand=expansion(mult=c(0.1, 0.1))) +
  labs(x = "", y = "Gene number", title = "The number of Rev-aging and Pro-aging DEGs in Brain") +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(text = element_text(colour = "black", size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        panel.background = element_rect(size = 1, color = "black"))
dev.off()