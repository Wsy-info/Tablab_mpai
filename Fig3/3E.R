### Figure 3
### E

library(Seurat)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)

### load data
load("brain.harmony_singlet_round2.RData")

### cellchat for each origin group
### Old
### Create a CellChat object
brain.harmony_old <- subset(brain.harmony, subset = origin == "Old")
cellchat_old <- brain.harmony_old@assays$RNA@data
meta <- brain.harmony_old@meta.data
cellchat_old <- createCellChat(object = cellchat_old, meta = meta, group.by = "celltype")

### set default labels
cellchat_old <- setIdent(cellchat_old,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_old@idents))

### set the ligand-receptor interaction database
cellchatDB <- CellChatDB.mouse
cellchat_old@DB <- cellchatDB

### processing the expression data for cell-cell communication
cellchat_old <- subsetData(cellchat_old)
cellchat_old <- identifyOverExpressedGenes(cellchat_old)
cellchat_old <- identifyOverExpressedInteractions(cellchat_old)
cellchat_old <- projectData(cellchat_old, PPI.mouse)

### inference of cell-cell communication network
cellchat_old <- computeCommunProb(cellchat_old)
cellchat_old <- computeCommunProbPathway(cellchat_old)
cellchat_old <- aggregateNet(cellchat_old)
cellchat_old <- netAnalysis_computeCentrality(cellchat_old, slot.name = "netP")
save(cellchat_old, file = "cellchat_old.RData")

### Young
### Create a CellChat object
brain.harmony_young <- subset(brain.harmony, subset = origin == "Young")
cellchat_young <- brain.harmony_young@assays$RNA@data
meta <- brain.harmony_young@meta.data
cellchat_young <- createCellChat(object = cellchat_young, meta = meta, group.by = "celltype")

### set default labels
cellchat_young <- setIdent(cellchat_young,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_young@idents))

### processing the expression data for cell-cell communication
cellchat_young <- subsetData(cellchat_young)
cellchat_young <- identifyOverExpressedGenes(cellchat_young)
cellchat_young <- identifyOverExpressedInteractions(cellchat_young)
cellchat_young <- projectData(cellchat_young, PPI.mouse)

### inference of cell-cell communication network
cellchat_young <- computeCommunProb(cellchat_young)
cellchat_young <- computeCommunProbPathway(cellchat_young)
cellchat_young <- aggregateNet(cellchat_young)
cellchat_young <- netAnalysis_computeCentrality(cellchat_young, slot.name = "netP")
save(cellchat_young, file = "cellchat_young.RData")

### MET
### Create a CellChat object
brain.harmony_met <- subset(brain.harmony, subset = origin == "MET")
cellchat_met <- brain.harmony_met@assays$RNA@data
meta <- brain.harmony_met@meta.data
cellchat_met <- createCellChat(object = cellchat_met, meta = meta, group.by = "celltype")

### set default labels
cellchat_met <- setIdent(cellchat_met,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_met@idents))

### set the ligand-receptor interaction database
cellchatDB <- CellChatDB.mouse
cellchat_met@DB <- cellchatDB

### processing the expression data for cell-cell communication
cellchat_met <- subsetData(cellchat_met)
cellchat_met <- identifyOverExpressedGenes(cellchat_met)
cellchat_met <- identifyOverExpressedInteractions(cellchat_met)
cellchat_met <- projectData(cellchat_met, PPI.mouse)

### inference of cell-cell communication network
cellchat_met <- computeCommunProb(cellchat_met)
cellchat_met <- computeCommunProbPathway(cellchat_met)
cellchat_met <- aggregateNet(cellchat_met)
cellchat_met <- netAnalysis_computeCentrality(cellchat_met, slot.name = "netP")
save(cellchat_met, file = "cellchat_met.RData")


### NR
### Create a CellChat object
brain.harmony_nr <- subset(brain.harmony, subset = origin == "NR")
cellchat_nr <- brain.harmony_nr@assays$RNA@data
meta <- brain.harmony_nr@meta.data
cellchat_nr <- createCellChat(object = cellchat_nr, meta = meta, group.by = "celltype")

### set default labels
cellchat_nr <- setIdent(cellchat_nr,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_nr@idents))

### set the ligand-receptor interaction database
cellchatDB <- CellChatDB.mouse
cellchat_nr@DB <- cellchatDB

### processing the expression data for cell-cell communication
cellchat_nr <- subsetData(cellchat_nr)
cellchat_nr <- identifyOverExpressedGenes(cellchat_nr)
cellchat_nr <- identifyOverExpressedInteractions(cellchat_nr)
cellchat_nr <- projectData(cellchat_nr, PPI.mouse)

### inference of cell-cell communication network
cellchat_nr <- computeCommunProb(cellchat_nr)
cellchat_nr <- computeCommunProbPathway(cellchat_nr)
cellchat_nr <- aggregateNet(cellchat_nr)
cellchat_nr <- netAnalysis_computeCentrality(cellchat_nr, slot.name = "netP")
save(cellchat_nr, file = "cellchat_nr.RData")


### D+Q
### Create a CellChat object
brain.harmony_dq <- subset(brain.harmony, subset = origin == "D+Q")
cellchat_dq <- brain.harmony_dq@assays$RNA@data
meta <- brain.harmony_dq@meta.data
cellchat_dq <- createCellChat(object = cellchat_dq, meta = meta, group.by = "celltype")


### set default labels
cellchat_dq <- setIdent(cellchat_dq,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_dq@idents))


### set the ligand-receptor interaction database
cellchatDB <- CellChatDB.mouse
cellchat_dq@DB <- cellchatDB

### processing the expression data for cell-cell communication
cellchat_dq <- subsetData(cellchat_dq)
cellchat_dq <- identifyOverExpressedGenes(cellchat_dq)
cellchat_dq <- identifyOverExpressedInteractions(cellchat_dq)
cellchat_dq <- projectData(cellchat_dq, PPI.mouse)

### inference of cell-cell communication network
cellchat_dq <- computeCommunProb(cellchat_dq)
cellchat_dq <- computeCommunProbPathway(cellchat_dq)
cellchat_dq <- aggregateNet(cellchat_dq)
cellchat_dq <- netAnalysis_computeCentrality(cellchat_dq, slot.name = "netP")
save(cellchat_dq, file = "cellchat_dq.RData")


### SPD
### Create a CellChat object
brain.harmony_spd <- subset(brain.harmony, subset = origin == "SPD")
cellchat_spd <- brain.harmony_spd@assays$RNA@data
meta <- brain.harmony_spd@meta.data
cellchat_spd <- createCellChat(object = cellchat_spd, meta = meta, group.by = "celltype")

### set default labels
cellchat_spd <- setIdent(cellchat_spd,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_spd@idents))

### set the ligand-receptor interaction database
cellchatDB <- CellChatDB.mouse
cellchat_spd@DB <- cellchatDB

### processing the expression data for cell-cell communication
cellchat_spd <- subsetData(cellchat_spd)
cellchat_spd <- identifyOverExpressedGenes(cellchat_spd)
cellchat_spd <- identifyOverExpressedInteractions(cellchat_spd)
cellchat_spd <- projectData(cellchat_spd, PPI.mouse)

### inference of cell-cell communication network
cellchat_spd <- computeCommunProb(cellchat_spd)
cellchat_spd <- computeCommunProbPathway(cellchat_spd)
cellchat_spd <- aggregateNet(cellchat_spd)
cellchat_spd <- netAnalysis_computeCentrality(cellchat_spd, slot.name = "netP")
save(cellchat_spd, file = "cellchat_spd.RData")


### MET+NR
### Create a CellChat object
brain.harmony_mn <- subset(brain.harmony, subset = origin == "MET+NR")
cellchat_mn <- brain.harmony_mn@assays$RNA@data
meta <- brain.harmony_mn@meta.data
cellchat_mn <- createCellChat(object = cellchat_mn, meta = meta, group.by = "celltype")

### set default labels
cellchat_mn <- setIdent(cellchat_mn,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_mn@idents))

### set the ligand-receptor interaction database
cellchatDB <- CellChatDB.mouse
cellchat_mn@DB <- cellchatDB

### processing the expression data for cell-cell communication
cellchat_mn <- subsetData(cellchat_mn)
cellchat_mn <- identifyOverExpressedGenes(cellchat_mn)
cellchat_mn <- identifyOverExpressedInteractions(cellchat_mn)
cellchat_mn <- projectData(cellchat_mn, PPI.mouse)

### inference of cell-cell communication network
cellchat_mn <- computeCommunProb(cellchat_mn)
cellchat_mn <- computeCommunProbPathway(cellchat_mn)
cellchat_mn <- aggregateNet(cellchat_mn)
cellchat_mn <- netAnalysis_computeCentrality(cellchat_mn, slot.name = "netP")
save(cellchat_mn, file = "cellchat_mn.RData")

### MET+SPD
### Create a CellChat object
brain.harmony_ms <- subset(brain.harmony, subset = origin == "MET+SPD")
cellchat_ms <- brain.harmony_ms@assays$RNA@data
meta <- brain.harmony_ms@meta.data
cellchat_ms <- createCellChat(object = cellchat_ms, meta = meta, group.by = "celltype")


### set default labels
cellchat_ms <- setIdent(cellchat_ms,  ident.use = "celltype")
### number of cells in each cell group
groupSize <- as.numeric(table(cellchat_ms@idents))


### set the ligand-receptor interaction database
cellchatDB <- CellChatDB.mouse
cellchat_ms@DB <- cellchatDB


### processing the expression data for cell-cell communication
cellchat_ms <- subsetData(cellchat_ms)
cellchat_ms <- identifyOverExpressedGenes(cellchat_ms)
cellchat_ms <- identifyOverExpressedInteractions(cellchat_ms)
cellchat_ms <- projectData(cellchat_ms, PPI.mouse)

### inference of cell-cell communication network
cellchat_ms <- computeCommunProb(cellchat_ms)
cellchat_ms <- computeCommunProbPathway(cellchat_ms)
cellchat_ms <- aggregateNet(cellchat_ms)
cellchat_ms <- netAnalysis_computeCentrality(cellchat_ms, slot.name = "netP")
save(cellchat_ms, file = "cellchat_ms.RData")



### merge all origin groups
object.list <- list("Old" = cellchat_old,
                    "Young" = cellchat_young,
                    "MET" = cellchat_met,
                    "NR" = cellchat_nr,
                    "D+Q" = cellchat_dq,
                    "SPD" = cellchat_spd,
                    "MN" = cellchat_mn,
                    "MS" = cellchat_ms)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
save(cellchat, file = "cellchat_merge_Brain.RData")





### analysis
### load cellchat data
load("cellchat_merge_Brain.RData")
load("cellchat_old.RData")
load("cellchat_young.RData")
load("cellchat_met.RData")
load("cellchat_nr.RData")
load("cellchat_dq.RData")
load("cellchat_spd.RData")
load("cellchat_mn.RData")
load("cellchat_ms.RData")

object.list <- list(Old = cellchat_old,
                    Young = cellchat_young,
                    MET = cellchat_met,
                    NR = cellchat_nr,
                    `D+Q` = cellchat_dq,
                    SPD = cellchat_spd,
                    MN = cellchat_mn,
                    MS = cellchat_ms
)

### set origin color
color_nmn <- c("#47a3bc", "#d8c91c", "#a3c8dc", "#a91f24", "#f09594",
               "#8e4b98", "#F2B71F", "#0160A6")
names(color_nmn) <- c("Old", "Young", "MET", "NR", "D+Q",
                      "SPD", "MN", "MS")


### visulization for prob
### get comparison for Young/MET/NR/D+Q/SPD/MN/MS group vs. Old group
### 1) get all differential ligand-receptor
all_diff_lr <- c()
for(comparison_2 in 2:8){
  o_y_1 <- netVisual_bubble(cellchat, comparison = c(1, comparison_2), angle.x = 45, remove.isolate = F)
  o_y_1 <- o_y_1$data
  diff_lr <- unique(o_y_1$interaction_name_2)
  all_diff_lr <- union(all_diff_lr, diff_lr)
}

### get prob matrix, row for all differential ligand-receptor, col for origin groups
all_diff_lr_matrix <- matrix(NA, length(all_diff_lr), 8)
rownames(all_diff_lr_matrix) <- all_diff_lr
colnames(all_diff_lr_matrix) <- names(color_nmn)

for(comparison_2 in 2:8){
  comparison_2_name <- names(color_nmn)[comparison_2]
  o_y_1 <- netVisual_bubble(cellchat, comparison = c(1, comparison_2), angle.x = 45, remove.isolate = F)
  o_y_1 <- o_y_1$data
  for(i in all_diff_lr){
    if(i %in% as.character(o_y_1$interaction_name_2)){
      o_y_1_i <- o_y_1 %>% filter(interaction_name_2 == i) %>% group_by(dataset) %>% reframe(sum(prob))
         if(!"Old" %in% o_y_1_i$dataset){
           all_diff_lr_matrix[i, comparison_2_name] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == comparison_2_name), 2])
         } else {
           if(!comparison_2_name %in% o_y_1_i$dataset){
             all_diff_lr_matrix[i, 1] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == "Old"), 2])
           } else {
             all_diff_lr_matrix[i, comparison_2_name] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == comparison_2_name), 2])
             all_diff_lr_matrix[i, 1] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == "Old"), 2])
           }
         }
    }

  }
}

### change NA to 0
all_diff_lr_matrix[which(is.na(all_diff_lr_matrix), arr.ind = T)] <- 0

### only retain Young, Old, MET, NR, D+Q and SPD groups
all_diff_lr_df_all_dan <- as.data.frame(all_diff_lr_matrix)
all_diff_lr_df_all_dan <- all_diff_lr_df_all_dan[, 1:6]
all_diff_lr_df_all_dan$freq_all <- apply(all_diff_lr_df_all_dan, 1, function(x){
  length(which(x != 0))
})


### get pathway type information
pathway_type <- read.xlsx("brain_pathway_type.xlsx")

### get GO information to find ligand-receptors in each pathway
GO_data <- clusterProfiler::get_GO_data("org.Mm.eg.db", "BP", "SYMBOL")
names(GO_data$PATHID2EXTID) <- GO_data$PATHID2NAME[match(names(GO_data$PATHID2EXTID), names(GO_data$PATHID2NAME))]

pathway_lr <- apply(pathway_type, 1, function(x){
  as.character(unlist(GO_data$PATHID2EXTID[names(GO_data$PATHID2EXTID) == x[1]]))
})
names(pathway_lr) <- pathway_type$Description

### comparsion <- c(1, 2), 1 > 2
pathway_circos <- function(comparsion){
  young_old_lr <- rownames(all_diff_lr_df_all_dan %>% filter(!!sym(comparsion[1]) > !!sym(comparsion[2])))
  young_old_score_df <- data.frame(pathway_type, origin = paste0(comparsion, collapse = '/'), score = 0)
  for(i in 1:length(pathway_type$Description)){
    score_all <- 0
    for(j in young_old_lr){
      l_r_i <- unlist(strsplit(j, split = "  - "))
      ### If one of the ligand-receptors is in the pathway
      if(l_r_i[1] %in% pathway_lr[[i]]){
        score_all <- score_all + all_diff_lr_df_all_dan[j, comparsion[1]]
      }
      if(l_r_i[2] %in% pathway_lr[[i]]){
        score_all <- score_all + all_diff_lr_df_all_dan[j, comparsion[1]]
      }
    }
    young_old_score_df[i, 4] <- score_all
  }
  return(young_old_score_df)
}


young_old_score_df <- pathway_circos(c("Young", "Old"))
young_old_score_df <- young_old_score_df[-c(38, 39), ]

met_old_score_df <- pathway_circos(c("MET", "Old"))
met_old_score_df <- met_old_score_df[-c(38, 39), ]

dq_old_score_df <- pathway_circos(c("D+Q", "Old"))
dq_old_score_df <- dq_old_score_df[-c(38, 39), ]

nr_old_score_df <- pathway_circos(c("NR", "Old"))
nr_old_score_df <- nr_old_score_df[-c(38, 39), ]

spd_old_score_df <- pathway_circos(c("SPD", "Old"))
spd_old_score_df <- spd_old_score_df[-c(38, 39), ]

dq_old_score_df_down <- pathway_circos(c("Old", "D+Q"))
dq_old_score_df_down <- dq_old_score_df_down[c(38, 39), ]

spd_old_score_df_down <- pathway_circos(c("Old", "SPD"))
spd_old_score_df_down <- spd_old_score_df_down[c(38, 39), ]

### combined
pathway_circos_combine <- rbind(young_old_score_df, rbind(met_old_score_df, rbind(dq_old_score_df, rbind(nr_old_score_df, rbind(spd_old_score_df, rbind(dq_old_score_df_down, spd_old_score_df_down))))))

node_info_1 <- pathway_circos_combine %>% group_by(Description) %>% reframe(Size = sum(score), type = type) %>% as.data.frame() %>% unique()
node_info <- pathway_circos_combine %>% group_by(origin) %>% reframe(Size = sum(score)) %>% as.data.frame()
node_info$type <- "Group"
colnames(node_info)[1] <- 'Node'
colnames(node_info_1)[1] <- 'Node'
node_info <- rbind(node_info, node_info_1)

### prepare for the network (cytoscapev3.10.2)
write.table(node_info, file = "pathway_circos_combine_node.txt", sep = "\t", quote = F, row.names = F)
write.table(pathway_circos_combine, file = "pathway_circos_combine.txt", sep = "\t", quote = F, row.names = F)


### barplot for ligand-receptors of each comparsion-group
bar_info <- pathway_circos_combine %>% group_by(origin, type) %>% reframe(Size = sum(score)) %>% as.data.frame()

ggplot(data = bar_info %>% filter(type == "development"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(high = '#85B94B', low = "white") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_development_bar_plot.pdf", width = 5, height = 2)

ggplot(data = bar_info %>% filter(type == "synapse-related"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(brewer.pal(8, "Purples")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_synapse_bar_plot.pdf", width = 5, height = 2)

ggplot(data = bar_info %>% filter(type == "Immune related"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(high = '#C568CC', low = "white") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_Immune_bar_plot.pdf", width = 5, height = 2)

ggplot(data = bar_info %>% filter(type == "brain function"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(high = '#E77CAF', low = "white") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_function_bar_plot.pdf", width = 5, height = 2)