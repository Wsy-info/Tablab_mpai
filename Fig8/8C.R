### Figure 8
### C

library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)
library(reshape2)

### load data
load("brain.harmony_singlet_round2.RData")
load("heart.harmony_singlet_round2.RData")
load("kidney.harmony_singlet_round2.RData")
load("liver.harmony_singlet_round2.RData")

### function for pathway enrichment
get_go_pathway <- function(gene_list){
  entrez_ids <- bitr(gene_list,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
  go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                         OrgDb = org.Mm.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)
  filter_go_results <- go_results@result %>% arrange(pvalue) %>% filter(pvalue < 0.05 & Count >= 3)
  return(filter_go_results)
}



common_cell <- c("EC", "PC", "Fib", "Mac1", "SMC")
### get upregulated and downregulated pathways in EC
i <- "EC"
brain.harmony.common <- subset(brain.harmony, subset = celltype == i)
heart.harmony.common <- subset(heart.harmony, subset = celltype %in% c("Vas_EC", "Endo_EC"))
kidney.harmony.common <- subset(kidney.harmony, subset = celltype == i)
liver.harmony.common <- subset(liver.harmony, subset = celltype == i)

### Figure S15
### get upregulated and downregulated pathways in Fib
# i <- "Fib"
# heart.harmony.common <- subset(heart.harmony, subset = celltype == i)
# kidney.harmony.common <- subset(kidney.harmony, subset = celltype == i)

Idents(brain.harmony.common) <- "origin"
Idents(heart.harmony.common) <- "origin"
Idents(kidney.harmony.common) <- "origin"
Idents(liver.harmony.common) <- "origin"


### 1. MET
### 1.1 Brain
met_brain <- FindMarkers(brain.harmony.common, ident.1 = "MET", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
met_brain_up <- met_brain %>% filter(avg_log2FC > 0)
met_brain_up_pathway <- get_go_pathway(rownames(met_brain_up))
### down
met_brain_down <- met_brain %>% filter(avg_log2FC < 0)
met_brain_down_pathway <- get_go_pathway(rownames(met_brain_down))

### 1.2 Heart
met_heart <- FindMarkers(heart.harmony.common, ident.1 = "MET", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
met_heart_up <- met_heart %>% filter(avg_log2FC > 0)
met_heart_up_pathway <- get_go_pathway(rownames(met_heart_up))
### down
met_heart_down <- met_heart %>% filter(avg_log2FC < 0)
met_heart_down_pathway <- get_go_pathway(rownames(met_heart_down))

### 1.3 Kidney
met_kidney <- FindMarkers(kidney.harmony.common, ident.1 = "MET", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
met_kidney_up <- met_kidney %>% filter(avg_log2FC > 0)
met_kidney_up_pathway <- get_go_pathway(rownames(met_kidney_up))
### down
met_kidney_down <- met_kidney %>% filter(avg_log2FC < 0)
met_kidney_down_pathway <- get_go_pathway(rownames(met_kidney_down))

### 1.4 Liver
met_liver <- FindMarkers(liver.harmony.common, ident.1 = "MET", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
met_liver_up <- met_liver %>% filter(avg_log2FC > 0)
met_liver_up_pathway <- get_go_pathway(rownames(met_liver_up))
### down
met_liver_down <- met_liver %>% filter(avg_log2FC < 0)
met_liver_down_pathway <- get_go_pathway(rownames(met_liver_down))



### 2. NR
### 2.1 Brain
nr_brain <- FindMarkers(brain.harmony.common, ident.1 = "NR", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
nr_brain_up <- nr_brain %>% filter(avg_log2FC > 0)
nr_brain_up_pathway <- get_go_pathway(rownames(nr_brain_up))
### down
nr_brain_down <- nr_brain %>% filter(avg_log2FC < 0)
nr_brain_down_pathway <- get_go_pathway(rownames(nr_brain_down))

### 2.2 Heart
nr_heart <- FindMarkers(heart.harmony.common, ident.1 = "NR", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
nr_heart_up <- nr_heart %>% filter(avg_log2FC > 0)
nr_heart_up_pathway <- get_go_pathway(rownames(nr_heart_up))
### down
nr_heart_down <- nr_heart %>% filter(avg_log2FC < 0)
nr_heart_down_pathway <- get_go_pathway(rownames(nr_heart_down))

### 2.3 Kidney
nr_kidney <- FindMarkers(kidney.harmony.common, ident.1 = "NR", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
nr_kidney_up <- nr_kidney %>% filter(avg_log2FC > 0)
nr_kidney_up_pathway <- get_go_pathway(rownames(nr_kidney_up))
### down
nr_kidney_down <- nr_kidney %>% filter(avg_log2FC < 0)
nr_kidney_down_pathway <- get_go_pathway(rownames(nr_kidney_down))

### 2.4 Liver
nr_liver <- FindMarkers(liver.harmony.common, ident.1 = "NR", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
nr_liver_up <- nr_liver %>% filter(avg_log2FC > 0)
nr_liver_up_pathway <- get_go_pathway(rownames(nr_liver_up))
### down
nr_liver_down <- nr_liver %>% filter(avg_log2FC < 0)
nr_liver_down_pathway <- get_go_pathway(rownames(nr_liver_down))


### 3. D+Q
### 3.1 Brain
dq_brain <- FindMarkers(brain.harmony.common, ident.1 = "D+Q", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
dq_brain_up <- dq_brain %>% filter(avg_log2FC > 0)
dq_brain_up_pathway <- get_go_pathway(rownames(dq_brain_up))
### down
dq_brain_down <- dq_brain %>% filter(avg_log2FC < 0)
dq_brain_down_pathway <- get_go_pathway(rownames(dq_brain_down))

### 3.2 Heart
dq_heart <- FindMarkers(heart.harmony.common, ident.1 = "D+Q", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
dq_heart_up <- dq_heart %>% filter(avg_log2FC > 0)
dq_heart_up_pathway <- get_go_pathway(rownames(dq_heart_up))
### down
dq_heart_down <- dq_heart %>% filter(avg_log2FC < 0)
dq_heart_down_pathway <- get_go_pathway(rownames(dq_heart_down))

### 3.3 Kidney
dq_kidney <- FindMarkers(kidney.harmony.common, ident.1 = "D+Q", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
dq_kidney_up <- dq_kidney %>% filter(avg_log2FC > 0)
dq_kidney_up_pathway <- get_go_pathway(rownames(dq_kidney_up))
### down
dq_kidney_down <- dq_kidney %>% filter(avg_log2FC < 0)
dq_kidney_down_pathway <- get_go_pathway(rownames(dq_kidney_down))

### 3.4 Liver
dq_liver <- FindMarkers(liver.harmony.common, ident.1 = "D+Q", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
dq_liver_up <- dq_liver %>% filter(avg_log2FC > 0)
dq_liver_up_pathway <- get_go_pathway(rownames(dq_liver_up))
### down
dq_liver_down <- dq_liver %>% filter(avg_log2FC < 0)
dq_liver_down_pathway <- get_go_pathway(rownames(dq_liver_down))


### 4. SPD
### 4.1 Brain
spd_brain <- FindMarkers(brain.harmony.common, ident.1 = "SPD", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
spd_brain_up <- spd_brain %>% filter(avg_log2FC > 0)
spd_brain_up_pathway <- get_go_pathway(rownames(spd_brain_up))
### down
spd_brain_down <- spd_brain %>% filter(avg_log2FC < 0)
spd_brain_down_pathway <- get_go_pathway(rownames(spd_brain_down))

### 4.2 Heart
spd_heart <- FindMarkers(heart.harmony.common, ident.1 = "SPD", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
spd_heart_up <- spd_heart %>% filter(avg_log2FC > 0)
spd_heart_up_pathway <- get_go_pathway(rownames(spd_heart_up))
### down
spd_heart_down <- spd_heart %>% filter(avg_log2FC < 0)
spd_heart_down_pathway <- get_go_pathway(rownames(spd_heart_down))

### 4.3 Kidney
spd_kidney <- FindMarkers(kidney.harmony.common, ident.1 = "SPD", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
spd_kidney_up <- spd_kidney %>% filter(avg_log2FC > 0)
spd_kidney_up_pathway <- get_go_pathway(rownames(spd_kidney_up))
### down
spd_kidney_down <- spd_kidney %>% filter(avg_log2FC < 0)
spd_kidney_down_pathway <- get_go_pathway(rownames(spd_kidney_down))

### 4.4 Liver
spd_liver <- FindMarkers(liver.harmony.common, ident.1 = "SPD", ident.2 = "Old", logfc.threshold = 0.25) %>% filter(p_val < 0.05)
### up
spd_liver_up <- spd_liver %>% filter(avg_log2FC > 0)
spd_liver_up_pathway <- get_go_pathway(rownames(spd_liver_up))
### down
spd_liver_down <- spd_liver %>% filter(avg_log2FC < 0)
spd_liver_down_pathway <- get_go_pathway(rownames(spd_liver_down))


### plot
pathway_result <- read.xlsx("common_celltype_pathway.xlsx", sheet = 1)
pathway_result_up <- pathway_result %>% filter(direction == "up")
up_pathway <- unique(pathway_result_up$Pathway_name)
pathway_result_down <- pathway_result %>% filter(direction == "down")
down_pathway <- unique(pathway_result_down$Pathway_name)


### upregulated pathway
up_pathway_plot <- rbind(met_brain_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "MET_Brain"),
                         nr_brain_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "NR_Brain"),
                         dq_brain_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "D+Q_Brain"),
                         spd_brain_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "SPD_Brain"),

                         met_heart_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "MET_Heart"),
                         nr_heart_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "NR_Heart"),
                         dq_heart_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "D+Q_Heart"),
                         spd_heart_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "SPD_Heart"),

                         met_kidney_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "MET_Kidney"),
                         nr_kidney_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "NR_Kidney"),
                         dq_kidney_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "D+Q_Kidney"),
                         spd_kidney_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "SPD_Kidney"),

                         met_liver_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "MET_Liver"),
                         nr_liver_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "NR_Liver"),
                         dq_liver_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "D+Q_Liver"),
                         spd_liver_up_pathway %>% filter(Description %in% up_pathway) %>% mutate(group = "SPD_Liver")
)
up_pathway_plot$Description <- factor(up_pathway_plot$Description,
                                      levels = rev(up_pathway))
up_pathway_plot$group <- factor(up_pathway_plot$group, levels = c("MET_Brain", "NR_Brain", "D+Q_Brain", "SPD_Brain",
                                                                  "MET_Heart", "NR_Heart", "D+Q_Heart", "SPD_Heart",
                                                                  "MET_Kidney", "NR_Kidney", "D+Q_Kidney", "SPD_Kidney",
                                                                  "MET_Liver", "NR_Liver", "D+Q_Liver", "SPD_Liver"
))
up_pathway_plot$pvalue <- -1 * log10(up_pathway_plot$pvalue)


p1 <- ggplot(data = up_pathway_plot, aes(x = group, y = Description)) +
  geom_point(aes(color = pvalue, size = Count)) +
  scale_color_distiller(name = "-log10(pvalue)", palette = "Reds", direction = 1) +
  scale_size_continuous(range = c(2, 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 12, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
        ) +
  labs(x = "", y = "", title = "Upregulated Pathway")


### downregulated pathway
down_pathway_plot <- rbind(met_brain_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "MET_Brain"),
                           nr_brain_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "NR_Brain"),
                           dq_brain_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "D+Q_Brain"),
                           spd_brain_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "SPD_Brain"),

                           met_heart_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "MET_Heart"),
                           nr_heart_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "NR_Heart"),
                           dq_heart_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "D+Q_Heart"),
                           spd_heart_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "SPD_Heart"),

                           met_kidney_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "MET_Kidney"),
                           nr_kidney_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "NR_Kidney"),
                           dq_kidney_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "D+Q_Kidney"),
                           spd_kidney_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "SPD_Kidney"),

                           met_liver_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "MET_Liver"),
                           nr_liver_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "NR_Liver"),
                           dq_liver_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "D+Q_Liver"),
                           spd_liver_down_pathway %>% filter(Description %in% down_pathway) %>% mutate(group = "SPD_Liver")
)
down_pathway_plot$Description <- factor(down_pathway_plot$Description, levels = rev(down_pathway))
down_pathway_plot$group <- factor(down_pathway_plot$group, levels = c("MET_Brain", "NR_Brain", "D+Q_Brain", "SPD_Brain",
                                                                      "MET_Heart", "NR_Heart", "D+Q_Heart", "SPD_Heart",
                                                                      "MET_Kidney", "NR_Kidney", "D+Q_Kidney", "SPD_Kidney",
                                                                      "MET_Liver", "NR_Liver", "D+Q_Liver", "SPD_Liver"
))
down_pathway_plot$pvalue <- -1 * log10(down_pathway_plot$pvalue)

p2 <- ggplot(data = down_pathway_plot, aes(x = group, y = Description)) +
  geom_point(aes(color = pvalue, size = Count)) +
  scale_color_distiller(name = "-log10(pvalue)", palette = "Blues", direction = 1) +
  scale_size_continuous(range = c(2, 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 12, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
        ) +
  labs(x = "", y = "", title = "Downregulated Pathway")
