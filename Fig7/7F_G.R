### Figure 7
### FG

library(dplyr)
library(pheatmap)
library(mixOmics)

### load metabolomics data
expression_metabolism <- read.table("metabolites_full_table.txt",
                                    sep = "\t", quote = "", comment.char = "", header = T, row.names = 1)

### Figure 7G
### PLS-DA
expr <- as.data.frame(t(expression_metabolism[, 2:37]))
pos_expr <- as.data.frame(t(expression_metabolism[2080:3795, 2:37]))
neg_expr <- as.data.frame(t(expression_metabolism[1:2079, 2:37]))
class <- rep(c("MET",  "NR", "Old", "D+Q", "SPD", "Young"), each = 6)

### build PLS-DA model
plsda.srbct <- plsda(expr, class, ncomp = 10)

set.seed(123)
perf.plsda.srbct <- mixOmics::perf(plsda.srbct, validation = "Mfold", folds = 3, progressBar = F, nrepeat = 10)
plot(perf.plsda.srbct, sd = T, legend.posotion = "horizontal")

### ncomp = 3
final.plsda.srbct <- plsda(expr, class, ncomp = 3)

### set origin color
color_nmn <- c("#47a3bc", "#d8c91c", "#a3c8dc", "#a91f24", "#f09594",
               "#8e4b98", "#F2B71F", "#0160A6")
names(color_nmn) <- c("Old", "Young", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")

### plot
plotIndiv(final.plsda.srbct, ind.names = FALSE, legend = TRUE,
          comp = c(1, 2), ellipse = TRUE,
          title = 'PLS-DA on SRBCT comp 1-2',
          col.per.group = color_nmn[c("D+Q", "MET", "NR", "Old", "SPD", "Young")]
)
ggsave("PLS-DA.pdf", width = 8, height = 6)



### Figure 7E
metabolite_par_anti <- data.frame(group = c("MET", "NR", "D+Q", "SPD", "MET", "NR", "D+Q", "SPD"),
                                  count = c(-191, -298, -475, -495, 263 + 184, 194 + 200, 190 + 123, 179 + 140))
metabolite_par_anti$group <- factor(metabolite_par_anti$group, levels = c("SPD", "D+Q", "NR", "MET"))
ggplot(data = metabolite_par_anti, aes(x = count, y = group, color = group)) +
  geom_segment(aes(x = 0, xend = count, y = group, yend = group),
               linetype = "solid",
               linewidth = 2) +
  geom_point(aes(size = abs(count))) +
  scale_size_continuous(range = c(5, 10)) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13)
  ) +
  labs(x = "Count", y = "")
ggsave("metabolite_size_bangbangtang.pdf", width = 6, height = 3)