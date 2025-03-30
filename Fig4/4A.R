### Figure 4
### A

library(dplyr)
library(Seurat)
library(reshape2)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

### get information of regulon
sce_SCENIC <- open_loom("brain.harmony.group_SCENIC_new.loom")
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name = 'RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(sce_SCENIC)
close_loom(sce_SCENIC)

### load information of snRNA-seq
load("brain.RData")
sub_regulonAUC <- regulonAUC[, match(col_ba, colnames(regulonAUC))]
dim(sub_regulonAUC)
identical(colnames(sub_regulonAUC), col_ba)

cellTypes <- data.frame(row.names = col_ba,
                        celltype = min_origin_ba)
### save
save(sub_regulonAUC, cellTypes, regulons,
     file = 'for_rss_and_visual.Rdata')

### average of regulons activity
### heatmap
selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes),
                       cellTypes[, selectedResolution])
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ]
dim(sub_regulonAUC)
regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells)
  rowMeans(getAUC(sub_regulonAUC)[, cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale = T))
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)
save(regulonActivity_byGroup_Scaled, file = "regulonActivity_byGroup_Scaled.RData")

### get regulons-regulatory network
source_target <- NULL
regulons_sub <- regulons
for (i in names(regulons_sub)) {
  gene <- regulons_sub[[i]]
  info <- data.frame("from" = rep(i,length(gene)),"to" = gene)
  source_target <- rbind(info, source_target)
}
source_target$TFs <- source_target$from
name <- unique(c(source_target$from,source_target$to))

### save net data
save(source_target, file = "source_target.RData")

### load data
load("brain.harmony_singlet_round2.RData")

### TF
tf <- data.frame(TF = unique(source_target$from))
all_tf <- apply(tf, 1, function(x){
  gsub("\\(\\+\\)", "", x)
})
source_target$from <- apply(source_target, 1, function(x){
  gsub("\\(\\+\\)", "", x[1])
})

### average ecpression of TF
tf_expression <- AverageExpression(brain.harmony, features = all_tf, add.ident = 'origin')$RNA
tf_expression <- melt(tf_expression)

### prepare dataframe for tf information
tf_expression$celltype <- apply(tf_expression, 1, function(x){
  y <- strsplit(x[2], split = "_")[[1]]
  if(length(y) == 2){
    return(y[1])
  } else {
    return(paste(y[1], y[2], sep = '_'))
  }
})

tf_expression$origin <- apply(tf_expression, 1, function(x){
  y <- strsplit(x[2], split = "_")[[1]]
  if(length(y) == 2){
    return(y[2])
  } else {
    return(y[3])
  }
})

all_cell_type <- levels(Idents(brain.harmony))
all_origin <- levels(brain.harmony@meta.data$origin)

all_tf_freq <- data.frame(Index = c("up_freq", "down_freq"), row.names = c("up_freq", "down_freq"))
for(i in all_tf){
  ### Old group
  tf_expression_i_old <- tf_expression %>% filter(Var1 == i & origin == 'Old')
  rownames(tf_expression_i_old) <- tf_expression_i_old$celltype
  tf_expression_i_old <- tf_expression_i_old[all_cell_type, ]

  for(j in all_origin[-2]){
    tf_expression_j <- tf_expression %>% filter(Var1 == i & origin == j)
    rownames(tf_expression_j) <- tf_expression_j$celltype
    tf_expression_j <- tf_expression_j[all_cell_type, ]

    tf_expression_i_j <- tf_expression_j
    tf_expression_i_j$value <- (tf_expression_j$value - tf_expression_i_old$value)
    tf_expression_i_j$fc <- log2(tf_expression_j$value / tf_expression_i_old$value)

    tf_freq_j <- tf_expression_i_j %>%
      filter(abs(fc) > 0.25) %>%
      reframe(up_freq = sum(value > 0), down_freq = sum(value < 0)) %>%
      t()
    colnames(tf_freq_j) <- paste(i, j)

    all_tf_freq <- cbind(all_tf_freq, tf_freq_j)
  }
}


Idents(brain.harmony) <- "origin"
tf_label_exp <- AverageExpression(brain.harmony, features = all_tf)$RNA

### TF color
tf_label <- apply(tf_label_exp, 1, function(x){
  (log2(x[1] / x[2]) + log2(x[3] / x[2]) + log2(x[4] / x[2]) + log2(x[5] / x[2]) + log2(x[6] / x[2]) + log2(x[7] / x[2]) + log2(x[8] / x[2])) / 7
})

### Up TF & Down TF
tf_label_binary <- apply(tf_label_exp, 1, function(x){
  ifelse(x[1] > x[2], "Down", "Up")
})
up_tf <- names(tf_label_binary)[which(tf_label_binary == "Up")]
down_tf <-  names(tf_label_binary)[which(tf_label_binary == "Down")]


### cell type
up_tf_celltype <- c()
up_number_all <- c()
up_number_all_sum <- c()
up_number_all_union <- c()
down_tf_celltype <- c()
down_number_all <- c()
down_number_all_sum <- c()
down_number_all_union <- c()

for(i in all_cell_type){
  brain.harmony_i <- subset(brain.harmony, subset = celltype == i)
  Idents(brain.harmony_i) <- "origin"
  ### each origin group
  for(j in levels(brain.harmony@meta.data$origin)[-2]){
    degs_i <- FindMarkers(brain.harmony_i, ident.1 = "Old", ident.2 = j)
    degs_i$gene <- rownames(degs_i)
    degs_i <- degs_i %>% filter(p_val < 0.05)

    ### up-regulated
    degs_i_up <- degs_i %>% filter(avg_log2FC > 0)
    up_tf_i <- degs_i_up$gene[which(degs_i_up$gene %in% up_tf)]
    if(length(up_tf_i) > 0){
      up_tf_celltype_i <- data.frame(TF = up_tf_i, celltype = i, origin = j)
      up_tf_celltype <- rbind(up_tf_celltype, up_tf_celltype_i)
      target_gene <- unique((source_target %>% filter(from %in% up_tf_i))$to)
      target_number <- intersect(target_gene, degs_i_up$gene)
      up_number_all <- c(up_number_all, target_number)
    } else {
      up_number_all <- c(up_number_all, 0)
    }

    ### down-regulated
    degs_i_down <- degs_i %>% filter(avg_log2FC < 0)
    down_tf_i <- degs_i_down$gene[which(degs_i_down$gene %in% down_tf)]
    if(length(down_tf_i) > 0){
      down_tf_celltype_i <- data.frame(TF = down_tf_i, celltype = i)
      down_tf_celltype <- rbind(down_tf_celltype, down_tf_celltype_i)
      target_gene <- unique((source_target %>% filter(from %in% down_tf_i))$to)
      target_number <- intersect(target_gene, degs_i_down$gene)
      down_number_all <- c(down_number_all, target_number)
    } else {
      down_number_all <- c(down_number_all, 0)
    }
  }
  up_number_all_sum <- c(up_number_all_sum, length(up_number_all))
  down_number_all_sum <- c(down_number_all_sum, length(down_number_all))
  up_number_all_union <- c(up_number_all_union, length(unique(up_number_all)))
  down_number_all_union <- c(down_number_all_union, length(unique(down_number_all)))
}

### network (for Cytoscape_v3.10.2)
up_tf_celltype_out <- unique(up_tf_celltype[, 1:2])
down_tf_celltype_out <- unique(down_tf_celltype[, 1:2])
write.table(up_tf_celltype_out, file = "up_tf_celltype.txt", sep = "\t", row.names = F, quote = F)
write.table(down_tf_celltype_out, file = "down_tf_celltype.txt", sep = "\t", row.names = F, quote = F)

### node information
### Up
### TF
node <- data.frame(node = names(table(up_tf_celltype$TF)),
                   size = as.numeric(table(up_tf_celltype$TF)),
                   type = "TF",
                   color = as.numeric(tf_label[names(table(up_tf_celltype$TF))])
)

### Celltype
names(up_number_all_union) <- all_cell_type
node_2 <- data.frame(node = names(up_number_all_union),
                   size = as.numeric(up_number_all_union),
                   type = "Cell type",
                   color = 0
)

node <- rbind(node, node_2)
write.table(node, file = "node_up.txt", sep = "\t", row.names = F, quote = F)


### Down
### TF
node <- data.frame(node = names(table(down_tf_celltype$TF)),
                   size = as.numeric(table(down_tf_celltype$TF)),
                   type = "TF",
                   color = as.numeric(tf_label[names(table(down_tf_celltype$TF))])
)

### Celltype
names(down_number_all_union) <- all_cell_type
node_2 <- data.frame(node = names(down_number_all_union),
                   size = as.numeric(down_number_all_union),
                   type = "Cell type",
                   color = 0
)

node <- rbind(node, node_2)
write.table(node, file = "node_down.txt", sep = "\t", row.names = F, quote = F)
