### Figure 4
### D

library(Seurat)
library(ggplot2)
library(dplyr)

### load data
load("brain.harmony_singlet_round2.RData")

### only retain Young, Old, MET, NR, D+Q and SPD groups
brain.harmony_dan <- subset(brain.harmony, subset = origin %in% c("Young", "Old", "MET", "NR", "D+Q", "SPD"))

### set origin color
color_nmn <- c("#d8c91c", "#47a3bc", "#a3c8dc", "#a91f24", "#f09594",
               "#8e4b98", "#F2B71F", "#0160A6")
names(color_nmn) <- c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")

### vlnplot
Idents(brain.harmony_dan) <- "origin"
VlnPlot(brain.harmony_dan, features = "Mef2c", pt.size = 0, cols = color_nmn) +
  theme_classic() +
  theme(text = element_text(size = 13, colour = "black"),
        panel.background = element_rect(size = 1, color = "black")) +
  labs(title = "Mef2c expression") +
  theme(axis.line = element_blank(),
        legend.position = "",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_summary(fun = mean,
               geom = "pointrange",
               size = 0.3)
ggsave("Mef2c_expression.pdf")