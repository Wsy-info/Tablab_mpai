### Figure 10

library(Seurat)
library(dplyr)
library(ggplot2)
library(gggibbous)
library(reshape2)

pathway_select <- read.csv("heatmap.csv")

color_origin <- c("#d8c91c", "#47a3bc", "#a3c8dc", "#a91f24", "#f09594",
                  "#8e4b98", "#fdbd10", "#0066b2")
origin_order <- c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
col_origin <- setNames(color_origin, origin_order)


### load heatmap data
### 1. Brain
pathway_select_brain <- pathway_select %>% filter(origin == "Brain")
rownames(pathway_select_brain) <- pathway_select_brain$Pathway_name

load("brain_heatmap.RData")
df_brain <- melt(df_final)
colnames(df_brain) <- c("drug_celltype", "value")
df_brain$Pathway_name <- rownames(df_final)
df_brain$pathway_type <- pathway_select_brain[df_brain$Pathway_name, "pathway_type"]
df_brain$direction <- pathway_select_brain[df_brain$Pathway_name, "direction"]
df_brain$origin <- apply(df_brain, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
df_brain$celltype <- apply(df_brain, 1, function(x){
  gsub("(Old_)|(MET_)|(NR_)|(D\\+Q_)|(SPD_)|(MET\\+NR_)|(MET\\+SPD_)", "", x[1])
})

moon_df <- matrix(0, 6, length(unique(pathway_select_brain$pathway_type)))
rownames(moon_df) <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
colnames(moon_df) <- unique(pathway_select_brain$pathway_type)
for(i in unique(pathway_select_brain$pathway_type)){
  df_brain_i <- df_brain %>% filter(pathway_type == i)
  df_brain_i_old <- df_brain_i %>% filter(origin == "Old")
  df_brain_i_met <- df_brain_i %>% filter(origin == "MET")
  df_brain_i_nr <- df_brain_i %>% filter(origin == "NR")
  df_brain_i_dq <- df_brain_i %>% filter(origin == "D+Q")
  df_brain_i_spd <- df_brain_i %>% filter(origin == "SPD")
  df_brain_i_mn <- df_brain_i %>% filter(origin == "MET+NR")
  df_brain_i_ms <- df_brain_i %>% filter(origin == "MET+SPD")

  moon_df["MET", i] <- sum((df_brain_i_met$value * df_brain_i_old$value) < 0) / length(df_brain_i_old$value)
  moon_df["NR", i] <- sum((df_brain_i_nr$value * df_brain_i_old$value) < 0) / length(df_brain_i_old$value)
  moon_df["D+Q", i] <- sum((df_brain_i_dq$value * df_brain_i_old$value) < 0) / length(df_brain_i_old$value)
  moon_df["SPD", i] <- sum((df_brain_i_spd$value * df_brain_i_old$value) < 0) / length(df_brain_i_old$value)
  moon_df["MET+NR", i] <- sum((df_brain_i_mn$value * df_brain_i_old$value) < 0) / length(df_brain_i_old$value)
  moon_df["MET+SPD", i] <- sum((df_brain_i_ms$value * df_brain_i_old$value) < 0) / length(df_brain_i_old$value)
}

rowMeans(moon_df)

moon_df_plot <- melt(moon_df)
colnames(moon_df_plot) <- c("origin", "pathway_type", "ratio")
moon_df_plot$pathway_type <- factor(moon_df_plot$pathway_type, levels = unique(pathway_select_brain$pathway_type))
moon_df_plot$origin <- factor(moon_df_plot$origin, levels = rev(c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")))

ggplot(moon_df_plot) +
  geom_point(aes(x = pathway_type, y = origin), size = 10, color = "grey", alpha = 0.6) +
  geom_moon(aes(x = pathway_type, y = origin, ratio = ratio), color = NA, right = F, fill = "#2386E8") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "None"
  ) +
  labs(x = "", y = "")
ggsave("brain_ratio.pdf", width = 8, height = 8)


### 2. Heart
pathway_select_heart <- pathway_select %>% filter(origin == "Heart")
rownames(pathway_select_heart) <- pathway_select_heart$Pathway_name

load("heart_heatmap.RData")
df_heart <- melt(df_final)
colnames(df_heart) <- c("drug_celltype", "value")
df_heart$Pathway_name <- rownames(df_final)
df_heart$pathway_type <- pathway_select_heart[df_heart$Pathway_name, "pathway_type"]
df_heart$direction <- pathway_select_heart[df_heart$Pathway_name, "direction"]
df_heart$origin <- apply(df_heart, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
df_heart$celltype <- apply(df_heart, 1, function(x){
  gsub("(Old_)|(MET_)|(NR_)|(D\\+Q_)|(SPD_)", "", x[1])
})

moon_df <- matrix(0, 4, length(unique(pathway_select_heart$pathway_type)))
rownames(moon_df) <- c("MET", "NR", "D+Q", "SPD")
colnames(moon_df) <- unique(pathway_select_heart$pathway_type)
for(i in unique(pathway_select_heart$pathway_type)){
  df_heart_i <- df_heart %>% filter(pathway_type == i)
  df_heart_i_old <- df_heart_i %>% filter(origin == "Old")
  df_heart_i_met <- df_heart_i %>% filter(origin == "MET")
  df_heart_i_nr <- df_heart_i %>% filter(origin == "NR")
  df_heart_i_dq <- df_heart_i %>% filter(origin == "D+Q")
  df_heart_i_spd <- df_heart_i %>% filter(origin == "SPD")

  moon_df["MET", i] <- sum((df_heart_i_met$value * df_heart_i_old$value) < 0) / length(df_heart_i_old$value)
  moon_df["NR", i] <- sum((df_heart_i_nr$value * df_heart_i_old$value) < 0) / length(df_heart_i_old$value)
  moon_df["D+Q", i] <- sum((df_heart_i_dq$value * df_heart_i_old$value) < 0) / length(df_heart_i_old$value)
  moon_df["SPD", i] <- sum((df_heart_i_spd$value * df_heart_i_old$value) < 0) / length(df_heart_i_old$value)
}

rowMeans(moon_df)

moon_df_plot <- melt(moon_df)
colnames(moon_df_plot) <- c("origin", "pathway_type", "ratio")
moon_df_plot$pathway_type <- factor(moon_df_plot$pathway_type, levels = unique(pathway_select_heart$pathway_type))
moon_df_plot$origin <- factor(moon_df_plot$origin, levels = rev(c("MET", "NR", "D+Q", "SPD")))

ggplot(moon_df_plot) +
  geom_point(aes(x = pathway_type, y = origin), size = 10, color = "grey", alpha = 0.6) +
  geom_moon(aes(x = pathway_type, y = origin, ratio = ratio), color = NA, right = F, fill = "#9456DB") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "None"
  ) +
  labs(x = "", y = "")
ggsave("heart_ratio.pdf", width = 8, height = 6)


### 3. Kidney
pathway_select_kidney <- pathway_select %>% filter(origin == "Kidney")
rownames(pathway_select_kidney) <- pathway_select_kidney$Pathway_name

load("kidney_heatmap.RData")
df_kidney <- melt(df_final)
colnames(df_kidney) <- c("drug_celltype", "value")
df_kidney$Pathway_name <- rownames(df_final)
df_kidney$pathway_type <- pathway_select_kidney[df_kidney$Pathway_name, "pathway_type"]
df_kidney$direction <- pathway_select_kidney[df_kidney$Pathway_name, "direction"]
df_kidney$origin <- apply(df_kidney, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
df_kidney$celltype <- apply(df_kidney, 1, function(x){
  gsub("(Old_)|(MET_)|(NR_)|(D\\+Q_)|(SPD_)", "", x[1])
})

moon_df <- matrix(0, 4, length(unique(pathway_select_kidney$pathway_type)))
rownames(moon_df) <- c("MET", "NR", "D+Q", "SPD")
colnames(moon_df) <- unique(pathway_select_kidney$pathway_type)[c(2, 1, 4, 5, 6, 3, 7)]
for(i in unique(pathway_select_kidney$pathway_type)){
  df_kidney_i <- df_kidney %>% filter(pathway_type == i)
  df_kidney_i_old <- df_kidney_i %>% filter(origin == "Old")
  df_kidney_i_met <- df_kidney_i %>% filter(origin == "MET")
  df_kidney_i_nr <- df_kidney_i %>% filter(origin == "NR")
  df_kidney_i_dq <- df_kidney_i %>% filter(origin == "D+Q")
  df_kidney_i_spd <- df_kidney_i %>% filter(origin == "SPD")

  moon_df["MET", i] <- sum((df_kidney_i_met$value * df_kidney_i_old$value) < 0) / length(df_kidney_i_old$value)
  moon_df["NR", i] <- sum((df_kidney_i_nr$value * df_kidney_i_old$value) < 0) / length(df_kidney_i_old$value)
  moon_df["D+Q", i] <- sum((df_kidney_i_dq$value * df_kidney_i_old$value) < 0) / length(df_kidney_i_old$value)
  moon_df["SPD", i] <- sum((df_kidney_i_spd$value * df_kidney_i_old$value) < 0) / length(df_kidney_i_old$value)
}

rowMeans(moon_df)

moon_df_plot <- melt(moon_df)
colnames(moon_df_plot) <- c("origin", "pathway_type", "ratio")
moon_df_plot$pathway_type <- factor(moon_df_plot$pathway_type, levels = unique(pathway_select_kidney$pathway_type)[c(2, 1, 4, 5, 6, 3, 7)])
moon_df_plot$origin <- factor(moon_df_plot$origin, levels = rev(c("MET", "NR", "D+Q", "SPD")))

ggplot(moon_df_plot) +
  geom_point(aes(x = pathway_type, y = origin), size = 10, color = "grey", alpha = 0.6) +
  geom_moon(aes(x = pathway_type, y = origin, ratio = ratio), color = NA, right = F, fill = "#ED4ABB") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "None"
  ) +
  labs(x = "", y = "")
ggsave("kidney_ratio.pdf", width = 8, height = 6)


### 4. Liver
pathway_select_liver <- pathway_select %>% filter(origin == "Liver")
rownames(pathway_select_liver) <- pathway_select_liver$Pathway_name

load("liver_heatmap.RData")
df_liver <- melt(df_final)
colnames(df_liver) <- c("drug_celltype", "value")
df_liver$Pathway_name <- rownames(df_final)
df_liver$pathway_type <- pathway_select_liver[df_liver$Pathway_name, "pathway_type"]
df_liver$direction <- pathway_select_liver[df_liver$Pathway_name, "direction"]
df_liver$origin <- apply(df_liver, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
df_liver$celltype <- apply(df_liver, 1, function(x){
  gsub("(Old_)|(MET_)|(NR_)|(D\\+Q_)|(SPD_)", "", x[1])
})

moon_df <- matrix(0, 4, length(unique(pathway_select_liver$pathway_type)))
rownames(moon_df) <- c("MET", "NR", "D+Q", "SPD")
colnames(moon_df) <- unique(pathway_select_liver$pathway_type)
for(i in unique(pathway_select_liver$pathway_type)){
  df_liver_i <- df_liver %>% filter(pathway_type == i)
  df_liver_i_old <- df_liver_i %>% filter(origin == "Old")
  df_liver_i_met <- df_liver_i %>% filter(origin == "MET")
  df_liver_i_nr <- df_liver_i %>% filter(origin == "NR")
  df_liver_i_dq <- df_liver_i %>% filter(origin == "D+Q")
  df_liver_i_spd <- df_liver_i %>% filter(origin == "SPD")

  moon_df["MET", i] <- sum((df_liver_i_met$value * df_liver_i_old$value) < 0) / length(df_liver_i_old$value)
  moon_df["NR", i] <- sum((df_liver_i_nr$value * df_liver_i_old$value) < 0) / length(df_liver_i_old$value)
  moon_df["D+Q", i] <- sum((df_liver_i_dq$value * df_liver_i_old$value) < 0) / length(df_liver_i_old$value)
  moon_df["SPD", i] <- sum((df_liver_i_spd$value * df_liver_i_old$value) < 0) / length(df_liver_i_old$value)
}

rowMeans(moon_df)

moon_df_plot <- melt(moon_df)
colnames(moon_df_plot) <- c("origin", "pathway_type", "ratio")
moon_df_plot$pathway_type <- factor(moon_df_plot$pathway_type, levels = unique(pathway_select_liver$pathway_type))
moon_df_plot$origin <- factor(moon_df_plot$origin, levels = rev(c("MET", "NR", "D+Q", "SPD")))

ggplot(moon_df_plot) +
  geom_point(aes(x = pathway_type, y = origin), size = 10, color = "grey", alpha = 0.6) +
  geom_moon(aes(x = pathway_type, y = origin, ratio = ratio), color = NA, right = F, fill = "#06948E") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "None"
  ) +
  labs(x = "", y = "")
ggsave("liver_ratio.pdf", width = 8, height = 6)



