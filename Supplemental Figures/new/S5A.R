### Figure
### S5A

library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)

### Brain
### load data
load("brain.harmony_singlet_round2.RData")

### get genes from GOBP_RIBOSOME_ASSEMBLY
gmtfile <- 'reference/msigdb.v2023.1.Mm.symbols.gmt'
geneset <- read.gmt(gmtfile)
gene <- geneset[geneset$term == "GOBP_RIBOSOME_ASSEMBLY",]$gene
brain.harmony_OY <- AddModuleScore(brain.harmony_OY, features = list(gene), name = 'RLR')

pdf("ribosome.pdf", width = 10)
VlnPlot(brain.harmony_OY,
        features= "RLR1",
        pt.size = 0,
        cols = rep(c("#0073C2FF", "#EEC000FF", "#CD534CFF", "#00CC33FF", "#984EA3",
                     "#8F7700FF", "#FCCDE5", "#C46DA0FF", "#F0027F", "#374E55FF",
                     "#FD6467", "#0B775E", "#FF59FF", "#FDC086", "#5B1A18",
                     "#69C8ECFF", "#3da8e0", "#CBDF9AFF", "#802268FF"), each = 2)) +
  theme_classic() +
  theme(text = element_text(size = 13, colour = "black"),
        panel.background = element_rect(size = 1, color = "black")) +
  RotatedAxis() +
  labs(title = "GOBP_RIBOSOME_ASSEMBLY") +
  theme(axis.line = element_blank(), legend.position = "", axis.title.x = element_blank()) +
  stat_summary(fun.data = "median_hilow", fun.args = list(conf.int = 0.5), geom = "pointrange", width = 0.2, size = 0.3)
dev.off()
