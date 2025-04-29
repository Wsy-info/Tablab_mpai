### Figure 7
### A

library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(reshape2)
library(circlize)
library(stringr)
library(Seurat)
library(clusterProfiler)
library(openxlsx)

### prepare pathway information
pathway_select <- read.csv("heatmap.csv")
GO_data <- clusterProfiler:::get_GO_data("org.Mm.eg.db", "BP", "SYMBOL")
names(GO_data$PATHID2EXTID) <- GO_data$PATHID2NAME[match(names(GO_data$PATHID2EXTID), names(GO_data$PATHID2NAME))]

### Liver
### load data
load("liver.harmony_singlet_round2.RData")
liver.harmony@meta.data$condition <- paste(liver.harmony@meta.data$origin, liver.harmony@meta.data$celltype, sep = "_")

### Change the order of columns in a plot
combinations <- expand.grid(origin = levels(liver.harmony@meta.data$origin), celltype = levels(liver.harmony@meta.data$celltype))
combinations$str <- apply(combinations, 1, function(row) {
  paste(row['origin'], row['celltype'], sep = "_")
})
sorted_combinations <- combinations[order(match(combinations$origin, levels(liver.harmony@meta.data$origin)), match(combinations$celltype, levels(liver.harmony@meta.data$celltype))), ]$str
liver.harmony@meta.data$condition <- factor(liver.harmony@meta.data$condition, levels = sorted_combinations)


### Old_up
select_pathway1 <- pathway_select %>% filter(origin == "Liver" & direction == "up")

score1 <- lapply(select_pathway1$Pathway_name, function(x){
  liver.harmony <- AddModuleScore(liver.harmony, features = list(as.character(unlist(GO_data$PATHID2EXTID[names(GO_data$PATHID2EXTID) == x]))))
  p1 <- VlnPlot(liver.harmony, features= "Cluster1", pt.size = 0, group.by = "condition")$data
  p1 %>% group_by(ident) %>% summarize(mean_value = mean(Cluster1))})
names(score1) <- select_pathway1$Pathway_name
df1 <- Reduce(rbind, score1)
df1$pathway <- factor(rep(select_pathway1$Pathway_name, each = 7 * 8), levels = select_pathway1$Pathway_name)


### old-young_up
mean_value <- by(df1, df1$pathway, function(x){a1 <- x[grepl("Young", x$ident),]
                                              x$mean_value - a1$mean_value})
df1$mean_value <- unlist(mean_value)
df1 <- df1[!grepl("Young", df1$ident), ]
df_up <- df1[!grepl("MET\\+NR|MET\\+SPD", df1$ident),]
df_old_up <- df_up[grepl("Old", df_up$ident), ]
df_old_up <- dcast(df_old_up, pathway ~ ident, value.var = "mean_value")


### drug-old_down(MET, NR, D+Q, SPD)
mean_value <- by(df1, df1$pathway, function(x){a1 <- x[grepl("Old", x$ident),]
                                              x$mean_value - a1$mean_value})
df1$mean_value <- unlist(mean_value)
df1 <- df1[!grepl("Old|Young", df1$ident), ]
df_drug_down <- df1[!grepl("MET\\+NR|MET\\+SPD", df1$ident),]
df_drug_down <- dcast(df_drug_down, pathway ~ ident, value.var = "mean_value")


### Old_down
select_pathway2 <- pathway_select %>% filter(origin == "Liver" & direction == "down")

score2 <- lapply(select_pathway2$Pathway_name, function(x){
  liver.harmony <- AddModuleScore(liver.harmony, features = list(as.character(unlist(GO_data$PATHID2EXTID[names(GO_data$PATHID2EXTID) == x]))))
  p1 <- VlnPlot(liver.harmony, features = "Cluster1", pt.size = 0, group.by = "condition")$data
  p1 %>% group_by(ident) %>% summarize(mean_value = mean(Cluster1))})
names(score2) <- select_pathway2
df2 <- Reduce(rbind, score2)
df2$pathway <- factor(rep(select_pathway2$Pathway_name, each = 7 * 8), levels = select_pathway2$Pathway_name)


### old-young_down
mean_value <- by(df2, df2$pathway, function(x){a1 <- x[grepl("Young", x$ident),]
                                              x$mean_value - a1$mean_value})
df2$mean_value <- unlist(mean_value)
df2 <- df2[!grepl("Young", df2$ident), ]
df_down <- df2[!grepl("MET\\+NR|MET\\+SPD", df2$ident),]
df_old_down <- df_down[grepl("Old", df_down$ident), ]
df_old_down <- dcast(df_old_down, pathway ~ ident, value.var = "mean_value")


### drug-old_up(MET, NR, D+Q, SPD)
mean_value <- by(df2, df2$pathway, function(x){a1 <- x[grepl("Old", x$ident),]
                                              x$mean_value - a1$mean_value})
df2$mean_value <- unlist(mean_value)
df2 <- df2[!grepl("Old|Young", df2$ident), ]
df_drug_up <- df2[!grepl("MET\\+NR|MET\\+SPD", df2$ident),]
df_drug_up <- dcast(df_drug_up, pathway ~ ident, value.var = "mean_value")


### old_up_down heatmap
df_final <- rbind(df_old_up, df_old_down)
df_final <- df_final[, -1]
rownames(df_final) <- c(select_pathway1$Pathway_name, select_pathway2$Pathway_name)


### annotation for column
color_liver <- c("#3CB371", "#CBDF9AFF", "#4682B4", "#483D8B", "#008B8B",
                 "#FF8C00", "#FD6467")
cluster_order <- c("HSC", "EC", "Cho", "Kup", "PC-Hep",
                   "MZ-Hep", "PP-Hep")
col_cluster <- setNames(color_liver, cluster_order)
col_annots <- data.frame(Group = str_match(colnames(df_final), "(^[^_]*)_.*$")[, 2],
                         Cell = str_match(colnames(df_final), "^[^_]*_(.*)$")[, 2])
col_annots$Group <- factor(col_annots$Group, levels = c("Old"))
col_annots$Cell <- factor(col_annots$Cell, levels =  c("HSC", "EC", "Cho", "Kup", "PC-Hep",
                                                       "MZ-Hep", "PP-Hep"))

group_colors <- c("Old" = "#47a3bc")
stage_colors <- col_cluster[match(unique(col_annots$Cell), names(col_cluster))]
top_annotation <- HeatmapAnnotation(df = col_annots,
                                    col = list(Group = group_colors, Cell = stage_colors))

### annotation for row
row_annots <- data.frame(Pathway = c(select_pathway1$pathway_type, select_pathway2$pathway_type))
row_annots$Pathway <- factor(row_annots$Pathway, levels = c("up_Signaling related", "Immune related", "down_Signaling related"))
pathway_colors <- c("up_Signaling related" = "#FF7F00",
                    "Immune related" = "#CAB2D6",
                    "down_Signaling related" = "#FF7F00"
)
left_annotation <- rowAnnotation(df = row_annots,
                                 col = list(Pathway = pathway_colors))

p1 <- Heatmap(as.matrix(df_final), border = T ,rect_gp = gpar(col = "white", lwd = 0.5),
        width = ncol(df_final) * unit(3, "mm"),
        height = nrow(df_final) * unit(5, "mm"),
        cluster_rows = F,
        cluster_columns = F,
        column_labels = col_annots$Cell,
        column_split = col_annots$Group,
        row_split = row_annots$Pathway,
        top_annotation = top_annotation,
        col = colorRamp2(c(-0.05, 0, 0.05), c("#0051b0", "white", "#ff6a6f")),
        left_annotation = left_annotation,
        heatmap_legend_param = list(
        title = "Difference score",
        direction = "horizontal"
     ),
     column_names_side = 'top',
     row_names_side = 'left',
     column_title_gp = gpar(fontsize = 15, fontface = "bold"),
     row_title_gp = gpar(fontsize = 0))





### drug_up_down(MET, NR, D+Q, SPD)
df_final <- rbind(df_drug_down, df_drug_up)
df_final <- df_final[, -1]
rownames(df_final) <- c(select_pathway1$Pathway_name, select_pathway2$Pathway_name)


### annotation for column
color_liver <- c("#3CB371", "#CBDF9AFF", "#4682B4", "#483D8B", "#008B8B",
                 "#FF8C00", "#FD6467")
cluster_order <- c("HSC", "EC", "Cho", "Kup", "PC-Hep",
                   "MZ-Hep", "PP-Hep")
col_cluster <- setNames(color_liver, cluster_order)
col_annots <- data.frame(Group = str_match(colnames(df_final), "(^[^_]*)_.*$")[, 2],
                         Cell = str_match(colnames(df_final), "^[^_]*_(.*)$")[, 2])
col_annots$Group <- factor(col_annots$Group, levels = c("MET", "NR", "D+Q", "SPD"))
col_annots$Cell <- factor(col_annots$Cell, levels =  c("HSC", "EC", "Cho", "Kup", "PC-Hep",
                                                       "MZ-Hep", "PP-Hep"))

group_colors <- c("MET" = "#a3c8dc",
                  "NR" = "#a91f24",
                  "D+Q" = "#f09594",
                  "SPD" = "#8e4b98")
stage_colors <- col_cluster[match(unique(col_annots$Cell), names(col_cluster))]
top_annotation <- HeatmapAnnotation(df = col_annots,
                                    col = list(Group = group_colors, Cell = stage_colors))


p2 <- Heatmap(as.matrix(df_final), border = T ,rect_gp = gpar(col = "white", lwd = 0.5),
        width = ncol(df_final) * unit(3, "mm"),
        height = nrow(df_final) * unit(5, "mm"),
        cluster_rows = F,
        cluster_columns = F,
        column_labels = col_annots$Cell,
        column_split = col_annots$Group,
        row_split = row_annots$Pathway,
        top_annotation = top_annotation,
        col = colorRamp2(c(-0.05, 0, 0.05), c("#0051b0", "white", "#ff6a6f")),
        left_annotation = left_annotation,
        heatmap_legend_param = list(
        title = "Difference score",
        direction = "horizontal"
     ),
     column_names_side = 'top',
     row_names_side = 'left',
     column_title_gp = gpar(fontsize = 15, fontface = "bold"),
     row_title_gp = gpar(fontsize = 0))


pdf("Liver_DEGs_heatmap.pdf", height = 20, width = 10)
p1 + p2
dev.off()
