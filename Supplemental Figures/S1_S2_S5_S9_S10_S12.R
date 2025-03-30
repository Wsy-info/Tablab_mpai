### Figure
### S1, S2, S5, S9, S10, S12

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(cowplot)
library(viridis)
library(ggthemes)
library(ggalluvial)
library(Scillus)
library(scRNAtoolVis)

### 1.Brain
load("brain.harmony_singlet_round2.RData")

### Figure S1A
### Vlnplot of qc
Idents(brain.harmony) <- "orig.ident"

color_sample <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F",
                  "#E41A1C","#E7298A", "#C3BC3F", "#FF7F00", "#F1788D",
                  "#223D6C", "#D20A13", "#FFD121", "#088247", "#B037C4",
                  "#58CDD9", "#7A142C", "#5D90BA", "#7CC767", "#91612D",
                  "#6E568C", "#E0367A", "#D8D155", "#64495D")

pdf("brain_QC_nFeature.pdf", width = 8)
VlnPlot(brain.harmony, features = c("nFeature_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Genes") + xlab("") + labs("nFeature_RNA")
dev.off()

pdf("brain_QC_nCount.pdf", width = 8)
VlnPlot(brain.harmony, features = c("nCount_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Counts") + xlab("") + labs("nCount_RNA")
dev.off()

pdf("brain_QC_mt.pdf", width = 8)
VlnPlot(brain.harmony, features = c("percent.mt"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black",size = 11,face = "bold"),
        legend.title = element_text(colour = "black",size = 15,face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("MT-genes(%)") + xlab("") + labs("percent.mt")
dev.off()


### Figure S2A
### DotPlot
gene <- c("Cldn11", "Mag", "Plp1", "C1qa", "Csf1r",
          "Cx3cr1", "Aqp4", "Gja1", "Slc1a3", "Pdgfra",
          "Cspg4", "Tnr", "Vtn", "Atp13a5", "Rgs5",
          "Gad1", "Gad2", "Slc32a1", "Pvalb", "Meis2",
          "Vip", "Sst", "Npy", "Slc17a6", "Slc17a7",
          "Slc18a2", "Slc6a3", "Th", "Drd1", "Tac1",
          "Drd2", "Penk", "Slc6a13", "Slc47a1", "Cldn5",
          "Flt1", "Prlr", "Ttr")
marker_df <- data.frame(gene = gene, cluster =  c(rep("OLG", 3), rep("MG", 3), rep("ASC", 3), rep("OPC", 3), rep("PC", 3), rep("InN_1", 4), rep("InN_2", 1), rep("InN_3", 1), rep("InN_4", 1), rep("InN_5", 1),
                                                  rep("ExN_1", 1), rep("ExN_2", 1),  rep("DoN", 3), rep("D1-MSN", 2), rep("D2-MSN", 2), rep("VLMC", 1), rep("ABC", 1), rep("EC", 2), rep("CPC", 2)))
marker_df$cluster <- factor(marker_df$cluster, levels = cluster_order)

p5 <- DotPlot(brain.harmony, features = marker_df$gene) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    RotatedAxis() +
    theme(
         panel.border = element_rect(color = "black"),
         panel.spacing = unit(1, "mm"),
         axis.title = element_blank(),
         axis.text.y = element_blank(),
         axis.text.x = element_text(face = "italic"),
         strip.text.x = element_text(angle = 90)
     ) +
  scale_size(range=c(2.1, 10))


df <- data.frame(x = 0, y =rev(levels(brain.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = rev(df$y))

p6 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = color_brain) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
        plot.margin = margin(r = 0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

pdf("brain_marker_dotplot.pdf", width = 15.8, height = 6.3)
plot_grid(p6, p5, align = "h", axis = "bt", rel_widths = c(1.5, 9))
dev.off()


### Figure S5
### Heatmap_marker
markers <- read.table("DEGs_celltype_brain.txt", sep = '\t', row.names = 1, header = T)
marker_top30 <- data.frame(markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))

pdf("brain_celltype_heatmap.pdf", width = 5, height = 9)
AverageHeatmap(object = brain.harmony,
               markerGene = marker_top30$gene,
               annoCol = TRUE, showRowNames = F,row_title = "Top30 marker genes",
               myanCol =  c("#0073C2FF", "#EEC000FF", "#CD534CFF", "#00CC33FF", "#984EA3",
                            "#8F7700FF", "#FCCDE5", "#C46DA0FF", "#F0027F", "#374E55FF",
                            "#FD6467", "#0B775E", "#FF59FF", "#FDC086", "#5B1A18",
                            "#69C8ECFF", "#3da8e0", "#CBDF9AFF", "#802268FF"))
dev.off()


### 2.Heart
### load data
load("heart.harmony_singlet_round2.RData")

### Figure S1B
### Vlnplot of qc
Idents(heart.harmony) <- "orig.ident"
color_sample <- c("#4E79A7", "#A0CBE8", "#8CD17D", "#499894", "#F28E2B",
                  "#FFBE7D", "#B6992D", "#E15759", "#FF9D9A", "#79706E",
                  "#D37295", "#FABFD2", "#B07AA1", "#B07AA1", "#D4A6C8",
                  "#9D7660", "#E58606", "#5D69B1", "#24796C", "#499894",
                  "#DAA51B", "#000000", "#99C945", "#ED645A")

pdf("heart_QC_nFeature.pdf", width = 8)
VlnPlot(heart.harmony, features = c("nFeature_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Genes") + xlab("") + labs("nFeature_RNA")
dev.off()

pdf("heart_QC_nCount.pdf", width = 8)
VlnPlot(heart.harmony, features = c("nCount_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Counts") + xlab("") + labs("nCount_RNA")
dev.off()

pdf("heart_QC_mt.pdf", width = 8)
VlnPlot(heart.harmony, features = c("percent.mt"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("MT-genes(%)") + xlab("") + labs("percent.mt")
dev.off()


### Figure S2B
### DotPlot
gene <- c("Abca8a", "Arhgap24", "Dcn", "Scn7a", "Csmd1",
          "Cdh19", "Trpc3", "Notch3", "Abcc9", "Myh11",
          "Mylk", "Kcnab1", "Egfl7", "Flt1" ,"Fabp4",
          "Npr3", "Flt4", "Mmrn1", "Reln", "Pkhd1l1",
          "Muc16", "Upk1b", "Myh6", "Nexn", "Tnnt2",
          "Myl4", "Myl7", "Nppa", "Cd74", "H2-Ab1",
          "Ccr2", "F13a1", "Mrc1", "Rbpj", "Skap1",
          "Il7r", "Grap2")
marker_df <- data.frame(gene = gene, cluster =  c(rep("Fib", 3), rep("nmSC", 3), rep("PC", 3), rep("SMC", 3), rep("Vas_EC", 3), rep("Endo_EC", 1), rep("LEC", 3), rep("EPI", 3), rep("VEN", 3), rep("ATR", 3), rep("Mac1", 3), rep("Mac2", 3), rep("TC", 3)))

p5 <- DotPlot(heart.harmony, features = marker_df$gene) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  RotatedAxis() +
  theme(
       panel.border = element_rect(color = "black"),
       panel.spacing = unit(1, "mm"),
       axis.title = element_blank(),
       axis.text.y = element_blank(),
       axis.text.x = element_text(face = "italic"),
       strip.text.x = element_text(angle = 90)
   ) +
  scale_size(range=c(2.1,10))


df <- data.frame(x = 0, y =rev(levels(heart.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = rev(df$y))

p6 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = color_heart) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
        plot.margin = margin(r = 0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

pdf("heart_marker_dotplot.pdf", width = 15.8, height = 6.3)
plot_grid(p6, p5, align = "h", axis = "bt", rel_widths = c(1.5, 9))
dev.off()


### Figure S9
### Heatmap_marker
markers <- read.table("DEGs_celltype_heart.txt", sep = '\t', row.names = 1, header = T)
marker_top30 <- data.frame(markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))

pdf("heart_celltype_heatmap.pdf", width = 5, height = 9)
AverageHeatmap(object = heart.harmony,
               markerGene = marker_top30$gene,
               annoCol = TRUE, showRowNames = F,
               row_title = "Top30 marker genes",
               myanCol = c('#cb4936', '#ff0000', '#984EA3', '#1616c8', '#0D63A5',
                           '#A12568', '#F7C394', '#F499C1', '#20b2aa', '#ade87c',
                           '#EE4590', '#0d6c0d', '#6A5ACD'))
dev.off()


### 3.Kidney
### load data
load("kidney.harmony_singlet_round2.RData")

### Figure S1C
### Vlnplot of qc
Idents(kidney.harmony) <- "orig.ident"
color_sample <- c('#4b6aa8', '#3ca0cf', '#c376a7', '#ad98c3', '#cea5c7',
                  '#53738c', '#a5a9b0', '#a78982', '#696a6c', '#92699e',
                  '#d69971', '#df5734', '#6c408e', '#ac6894', '#d4c2db',
                  '#537eb7', '#83ab8e', '#ece399', '#405993', '#cc7f73',
                  '#b95055', '#d5bb72', '#bc9a7f', '#e6b884')

pdf("kidney_QC_nFeature.pdf", width = 8)
VlnPlot(kidney.harmony, features = c("nFeature_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Genes") + xlab("") + labs("nFeature_RNA")
dev.off()

pdf("kidney_QC_nCount.pdf", width = 8)
VlnPlot(kidney.harmony, features = c("nCount_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Counts") + xlab("") + labs("nCount_RNA")
dev.off()


pdf("kidney_QC_mt.pdf", width = 8)
VlnPlot(kidney.harmony, features = c("percent.mt"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("MT-genes(%)") + xlab("") + labs("percent.mt")
dev.off()

### Figure S2C
### DotPlot
gene <- c("Plpp1", "Egfl7", "Flt1", "Abca8a", "Arhgap24",
          "Pdgfra", "Myh11", "Kcnab1", "Nphs1", "Nphs2",
          "Slc5a2", "Slc5a12", "Gatm", "Cyp2e1", "Cyp4b1",
          "Slc13a3", "Slc7a13", "Cyp7b1", "Atp11a", "Slc12a1",
          "Umod", "Egf", "Slc12a3", "Slc8a1", "S100g",
          "Klk1", "Atp6v1g3", "Atp6v0d2", "Aqp2", "Fxyd4",
          "Hsd11b2", "Insrr", "Cd74", "H2-Ab1", "H2-Eb1")
marker_df <- data.frame(gene = gene, cluster =  c(rep("EC", 3), rep("Fib", 3), rep("SMC", 2), rep("Podo", 2), rep("PT-S1", 3), rep("PT-S2", 3), rep("PT-S3", 3), rep("LOH", 3), rep("DCT1", 1), rep("DCT2", 1), rep("CNT", 2),rep("CD-IC", 2), rep("CD-PC", 3), rep("CD-trans", 1), rep("Mac1", 3)))

p5 <- DotPlot(kidney.harmony, features = marker_df$gene) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  RotatedAxis() +
  theme(
       panel.border = element_rect(color = "black"),
       panel.spacing = unit(1, "mm"),
       axis.title = element_blank(),
       axis.text.y = element_blank(),
       axis.text.x = element_text(face = "italic"),
       strip.text.x = element_text(angle = 90)
   ) +
  scale_size(range=c(2.1, 10))


df <- data.frame(x = 0, y =rev(levels(kidney.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = rev(df$y))

p6 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = color_kidney) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
        plot.margin = margin(r = 0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

pdf("kidney_marker_dotplot.pdf", width = 15.8, height = 6.3)
plot_grid(p6, p5, align ="h", axis="bt", rel_widths = c(1.5, 9))
dev.off()


### Figure S10
### Heatmap_marker
markers <- read.table("DEGs_celltype_kidney.txt", sep = '\t', row.names = 1, header = T)
marker_top30 <- data.frame(markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))

pdf("kidney_celltype_heatmap.pdf", width = 5, height = 9)
AverageHeatmap(object = kidney.harmony,
               markerGene = marker_top30$gene,
               annoCol = TRUE, showRowNames = F, row_title = "Top30 marker genes",
               myanCol = c('#CBDF9AFF', '#cb4936', '#1616c8', '#007ABA', '#8B58A4',
                           '#DF75AE', '#00B7CA', '#A8C7E9', '#E91E25', '#F37121',
                           '#FBB36E', '#F58D93', '#00A163', '#8CCA7C', '#EE4590'))
dev.off()



### 4.Liver
### load data
load("liver.harmony_singlet_round2.RData")

### Figure S1D
### Vlnplot of qc
Idents(liver.harmony) <- "orig.ident"
color_sample <- c('#da6f6d', '#ebb1a4', '#a44e89', '#a9c2cb', '#b85292',
                  '#6d6fa0', '#8d689d', '#c8c7e1', '#d25774', '#c49abc',
                  '#927c9a', '#3674a2', '#9f8d89', '#72567a', '#63a3b8',
                  '#c4daec', '#61bada', '#b7deea', '#e29eaf', '#4490c4',
                  '#e6e2a3', '#de8b36', '#c4612f', '#9a70a8')

pdf("liver_QC_nFeature.pdf", width = 8)
VlnPlot(liver.harmony, features = c("nFeature_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Genes") + xlab("") + labs("nFeature_RNA")
dev.off()

pdf("liver_QC_nCount.pdf", width = 8)
VlnPlot(liver.harmony, features = c("nCount_RNA"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("Counts") + xlab("") + labs("nCount_RNA")
dev.off()

pdf("liver_QC_mt.pdf", width = 8)
VlnPlot(liver.harmony, features = c("percent.mt"), pt.size = 0, cols = color_sample) +
  theme(legend.text = element_text(colour = "black", size = 11, face = "bold"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) +
  ylab("MT-genes(%)") + xlab("") + labs("percent.mt")
dev.off()


### Figure S2D
### DotPlot
gene <- c("Ank3", "Nrxn1", "Reln", "Ptprb", "Egfl7",
          "Stab2", "Thsd4", "Pkhd1", "Cd5l", "Fyb",
          "Adgb", "Apob", "Mug1", "Cyp2e1", "Ces3a",
          "Hamp", "Cyp2f2", "Hal", "Gls2")
marker_df <- data.frame(gene = gene, cluster =  c(rep("HSC", 3), rep("EC", 3), rep("Cho", 2), rep("Kup", 3), rep("PC-Hep", 4), rep("MZ-Hep", 1), rep("PP-Hep", 3)))

p5 <- DotPlot(liver.harmony, features = marker_df$gene) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  RotatedAxis() +
  theme(
       panel.border = element_rect(color = "black"),
       panel.spacing = unit(1, "mm"),
       axis.title = element_blank(),
       axis.text.y = element_blank(),
       axis.text.x = element_text(face = "italic"),
       strip.text.x = element_text(angle = 90)
   ) +
  scale_size(range=c(2.4, 10))


df <- data.frame(x = 0, y = rev(levels(liver.harmony)), stringsAsFactors = F)
df$y <- factor(df$y, levels = rev(df$y))

p6 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = color_liver) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
        plot.margin = margin(r = 0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

pdf("liver_marker_dotplot.pdf", width = 15.8, height = 6.3)
plot_grid(p6, p5, align = "h", axis = "bt", rel_widths = c(1.5, 9))
dev.off()


### Figure S12
### Heatmap_marker
markers <- read.table("DEGs_celltype_liver.txt", sep = '\t', row.names = 1, header = T)
marker_top30 <- data.frame(markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))

pdf("liver_celltype_heatmap.pdf", width = 5, height = 9)
AverageHeatmap(object = liver.harmony,
               markerGene = marker_top30$gene,
               annoCol = TRUE, showRowNames = F,row_title = "Top30 marker genes",
               myanCol =  c("#3CB371", "#CBDF9AFF", "#4682B4", "#483D8B", "#008B8B",
                            "#FF8C00", "#FD6467"))
dev.off()