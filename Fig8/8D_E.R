### Figure 8
### DE

### Brain
### load data
load("brain.harmony_singlet_round2.RData")

### set origin color
color_nmn <- c("#47a3bc", "#d8c91c", "#a3c8dc", "#a91f24", "#f09594",
               "#8e4b98", "#F2B71F", "#0160A6")
names(color_nmn) <- c("Old", "Young", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")

### Transcriptom noise
educlean_cal <- function(a, b){
  sqrt(sum((a - b)^2))
}

### for each sample
all_tn <- c()
for (i in levels(brain.harmony@meta.data$orig.ident)){
  brain.harmony_i <- subset(brain.harmony, subset = orig.ident == i)
  Idents(brain.harmony_i) <- "celltype"
  aver_expr <- AverageExpression(brain.harmony_i)$RNA
  expr_i <- as.matrix(brain.harmony_i@assays$RNA@data)
  all_tn_i <- c()
  for(j in 1:(dim(expr_i)[2])){
    all_tn_i <- c(all_tn_i, educlean_cal(as.numeric(expr_i[, j]), as.numeric(aver_expr[, brain.harmony_i@meta.data[j, "celltype"]])))
  }
  all_tn_i <- data.frame(sample = i,
                         origin = unique(brain.harmony_i@meta.data$origin),
                         tn = all_tn_i
                         )
  all_tn <- rbind(all_tn, all_tn_i)
}


### plot for MET+NR
ggplot(data = all_tn %>% filter(origin %in% c("Young", "Old", "MET", "NR", "MET+NR")),
       aes(x = origin, y = tn, color = origin)) +
  geom_boxplot() +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "none"
  ) +
  labs(x = "", y = "Transcriptom Noise")
ggsave("TN_MN.pdf", width = 6, height = 6)

### plot for MET+SPD
ggplot(data = all_tn %>% filter(origin %in% c("Young", "Old", "MET", "SPD", "MET+SPD")),
       aes(x = origin, y = tn, color = origin)) +
  geom_boxplot() +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "none"
  ) +
  labs(x = "", y = "Transcriptom Noise")
ggsave("TN_MS.pdf", width = 6, height = 6)