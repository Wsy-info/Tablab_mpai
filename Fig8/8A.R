### Figure 8
### A

library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)
library(reshape2)

load("brain.harmony_singlet_round2.RData")
load("heart.harmony_singlet_round2.RData")
load("kidney.harmony_singlet_round2.RData")
load("liver.harmony_singlet_round2.RData")

### get some pathways not in heatmap for some organs
### load heatmap selected pathways
pathway_heatmap <- read.csv("heatmap.csv")
pathway_order <- read.csv("heatmap2.csv", header = F)
colnames(pathway_order) <- c("pathway_type", "Pathway_name", "direction")
pathway_order <- unique(pathway_order)

brain_old_young <- read.csv("../matrix/brain_old_young.csv")
heart_old_young <- read.csv("../matrix/heart_old_young.csv")
kidney_old_young <- read.csv("../matrix/kidney_old_young.csv")
liver_old_young <- read.csv("../matrix/liver_old_young.csv")

rownames(pathway_order) <- pathway_order$Pathway_name
pathway_order$organ <- apply(pathway_order, 1, function(x){
  organ_label <- ""
  if(x[3] == "up"){
    if(x[2] %in% (brain_old_young %>% filter(UP_freq > DOWN_freq))$X){
      organ_label <- "Brain"
    }
    if(x[2] %in% (heart_old_young %>% filter(UP_freq > DOWN_freq))$X){
      organ_label <- ifelse(organ_label == "", "Heart", paste0(organ_label, ";", "Heart"))
    }
    if(x[2] %in% (kidney_old_young %>% filter(UP_freq > DOWN_freq))$X){
      organ_label <- ifelse(organ_label == "", "Kideny", paste0(organ_label, ";", "Kidney"))
    }
    if(x[2] %in% (liver_old_young %>% filter(UP_freq > DOWN_freq))$X){
      organ_label <- ifelse(organ_label == "", "Liver", paste0(organ_label, ";", "Liver"))
    }
  } else {
    if(x[2] %in% (brain_old_young %>% filter(UP_freq < DOWN_freq))$X){
      organ_label <- "Brain"
    }
    if(x[2] %in% (heart_old_young %>% filter(UP_freq < DOWN_freq))$X){
      organ_label <- ifelse(organ_label == "", "Heart", paste0(organ_label, ";", "Heart"))
    }
    if(x[2] %in% (kidney_old_young %>% filter(UP_freq < DOWN_freq))$X){
      organ_label <- ifelse(organ_label == "", "Kideny", paste0(organ_label, ";", "Kidney"))
    }
    if(x[2] %in% (liver_old_young %>% filter(UP_freq < DOWN_freq))$X){
      organ_label <- ifelse(organ_label == "", "Liver", paste0(organ_label, ";", "Liver"))
    }
  }
  return(organ_label)
})

pathway_order_df <- c()
for(i in 1:nrow(pathway_order)){
  organ_i <- strsplit(pathway_order$organ[i], split = ";")[[1]]
  pathway_order_df_i <- data.frame(pathway_type = pathway_order$pathway_type[i],
                                   Pathway_name = pathway_order$Pathway_name[i],
                                   organ = organ_i,
                                   direction = pathway_order$direction[i]
  )
  pathway_order_df <- rbind(pathway_order_df, pathway_order_df_i)
}
colnames(pathway_heatmap)[3] <- "organ"
pathway_heatmap <- pathway_heatmap %>% filter(Pathway_name %in% pathway_order_df$Pathway_name)

pathway_order_df <- rbind(pathway_order_df, pathway_heatmap)
pathway_order_df <- unique(pathway_order_df[, -1])
pathway_order_df$pathway_type <- pathway_order[pathway_order_df$Pathway_name, "pathway_type"]


pathway_select_up <- pathway_order_df %>% filter(direction == "up")
pathway_select_down <- pathway_order_df %>% filter(direction == "down")

### Brain
load("brain_heatmap.RData")
### up
pathway_select_brain <- pathway_select_up %>% filter(organ == "Brain")

### calculate mean_value and ratio
df_final_brain <- df_final[unique(pathway_select_brain$Pathway_name), ]
rownames(pathway_select_brain) <- unique(pathway_select_brain$Pathway_name)

df_brain <- melt(df_final_brain)
df_brain$Pathway_name <- rownames(df_final_brain)
df_brain$pathway_type <- pathway_select_brain[df_brain$Pathway_name, "pathway_type"]
df_brain$direction <- pathway_select_brain[df_brain$Pathway_name, "direction"]
df_brain$drug <- apply(df_brain, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_brain)[1:2] <- c("drug_celltype", "value")

df_brain_ratio_up <- df_brain %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value < 0), mean_score = mean(value)) %>%
  mutate(organ = "Brain", direction = "up")


### down
pathway_select_brain <- pathway_select_down %>% filter(organ == "Brain")

### calculate mean_value and ratio
df_final_brain <- df_final[pathway_select_brain$Pathway_name, ]
rownames(pathway_select_brain) <- pathway_select_brain$Pathway_name

df_brain <- melt(df_final_brain)
df_brain$Pathway_name <- rownames(df_final_brain)
df_brain$pathway_type <- pathway_select_brain[df_brain$Pathway_name, "pathway_type"]
df_brain$direction <- pathway_select_brain[df_brain$Pathway_name, "direction"]
df_brain$drug <- apply(df_brain, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_brain)[1:2] <- c("drug_celltype", "value")

df_brain_ratio_down <- df_brain %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value > 0), mean_score = mean(value)) %>%
  mutate(organ = "Brain", direction = "down")

df_brain_ratio <- rbind(df_brain_ratio_up, df_brain_ratio_down)
save(df_brain_ratio, file = "df_brain_ratio.RData")





### Heart
load("heart_heatmap.RData")
### up
pathway_select_heart <- pathway_select_up %>% filter(organ == "Heart")

### calculate mean_value and ratio
df_final_heart <- df_final[pathway_select_heart$Pathway_name, ]
rownames(pathway_select_heart) <- pathway_select_heart$Pathway_name

df_heart <- melt(df_final_heart)
df_heart$Pathway_name <- rownames(df_final_heart)
df_heart$pathway_type <- pathway_select_heart[df_heart$Pathway_name, "pathway_type"]
df_heart$direction <- pathway_select_heart[df_heart$Pathway_name, "direction"]
df_heart$drug <- apply(df_heart, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_heart)[1:2] <- c("drug_celltype", "value")

df_heart_ratio_up <- df_heart %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value < 0), mean_score = mean(value)) %>%
  mutate(organ = "Heart", direction = "up")


### down
pathway_select_heart <- pathway_select_down %>% filter(organ == "Heart")

### calculate mean_value and ratio
df_final_heart <- df_final[pathway_select_heart$Pathway_name, ]
rownames(pathway_select_heart) <- pathway_select_heart$Pathway_name

df_heart <- melt(df_final_heart)
df_heart$Pathway_name <- rownames(df_final_heart)
df_heart$pathway_type <- pathway_select_heart[df_heart$Pathway_name, "pathway_type"]
df_heart$direction <- pathway_select_heart[df_heart$Pathway_name, "direction"]
df_heart$drug <- apply(df_heart, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_heart)[1:2] <- c("drug_celltype", "value")

df_heart_ratio_down <- df_heart %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value > 0), mean_score = mean(value)) %>%
  mutate(organ = "Heart", direction = "down")

df_heart_ratio <- rbind(df_heart_ratio_up, df_heart_ratio_down)
save(df_heart_ratio, file = "df_heart_ratio.RData")



### Kidney
load("kidney_heatmap.RData")
### up
pathway_select_kidney <- pathway_select_up %>% filter(organ == "Kidney")

### calculate mean_value and ratio
df_final_kidney <- df_final[pathway_select_kidney$Pathway_name, ]
rownames(pathway_select_kidney) <- pathway_select_kidney$Pathway_name

df_kidney <- melt(df_final_kidney)
df_kidney$Pathway_name <- rownames(df_final_kidney)
df_kidney$pathway_type <- pathway_select_kidney[df_kidney$Pathway_name, "pathway_type"]
df_kidney$direction <- pathway_select_kidney[df_kidney$Pathway_name, "direction"]
df_kidney$drug <- apply(df_kidney, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_kidney)[1:2] <- c("drug_celltype", "value")

df_kidney_ratio_up <- df_kidney %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value < 0), mean_score = mean(value)) %>%
  mutate(organ = "Kidney", direction = "up")


### down
pathway_select_kidney <- pathway_select_down %>% filter(organ == "Kidney")

### calculate mean_value and ratio
df_final_kidney <- df_final[pathway_select_kidney$Pathway_name, ]
rownames(pathway_select_kidney) <- pathway_select_kidney$Pathway_name

df_kidney <- melt(df_final_kidney)
df_kidney$Pathway_name <- rownames(df_final_kidney)
df_kidney$pathway_type <- pathway_select_kidney[df_kidney$Pathway_name, "pathway_type"]
df_kidney$direction <- pathway_select_kidney[df_kidney$Pathway_name, "direction"]
df_kidney$drug <- apply(df_kidney, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_kidney)[1:2] <- c("drug_celltype", "value")

df_kidney_ratio_down <- df_kidney %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value > 0), mean_score = mean(value)) %>%
  mutate(organ = "Kidney", direction = "down")

df_kidney_ratio <- rbind(df_kidney_ratio_up, df_kidney_ratio_down)
save(df_kidney_ratio, file = "df_kidney_ratio.RData")


### Liver
load("liver_heatmap.RData")
### up
pathway_select_liver <- pathway_select_up %>% filter(organ == "Liver")

### calculate mean_value and ratio
df_final_liver <- df_final[pathway_select_liver$Pathway_name, ]
rownames(pathway_select_liver) <- pathway_select_liver$Pathway_name

df_liver <- melt(df_final_liver)
df_liver$Pathway_name <- rownames(df_final_liver)
df_liver$pathway_type <- pathway_select_liver[df_liver$Pathway_name, "pathway_type"]
df_liver$direction <- pathway_select_liver[df_liver$Pathway_name, "direction"]
df_liver$drug <- apply(df_liver, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_liver)[1:2] <- c("drug_celltype", "value")

df_liver_ratio_up <- df_liver %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value < 0), mean_score = mean(value)) %>%
  mutate(organ = "Liver", direction = "up")


### down
pathway_select_liver <- pathway_select_down %>% filter(organ == "Liver")

### calculate mean_value and ratio
df_final_liver <- df_final[pathway_select_liver$Pathway_name, ]
rownames(pathway_select_liver) <- pathway_select_liver$Pathway_name

df_liver <- melt(df_final_liver)
df_liver$Pathway_name <- rownames(df_final_liver)
df_liver$pathway_type <- pathway_select_liver[df_liver$Pathway_name, "pathway_type"]
df_liver$direction <- pathway_select_liver[df_liver$Pathway_name, "direction"]
df_liver$drug <- apply(df_liver, 1, function(x){
  strsplit(x[1], split = "_")[[1]][1]
})
colnames(df_liver)[1:2] <- c("drug_celltype", "value")

df_liver_ratio_down <- df_liver %>% group_by(drug, Pathway_name, pathway_type) %>% summarise(ratio = mean(value > 0), mean_score = mean(value)) %>%
  mutate(organ = "Liver", direction = "down")

df_liver_ratio <- rbind(df_liver_ratio_up, df_liver_ratio_down)
save(df_liver_ratio, file = "df_liver_ratio.RData")


load("df_brain_ratio.RData")
load("df_heart_ratio.RData")
load("df_kidney_ratio.RData")
load("df_liver_ratio.RData")
### combine
df_all <- rbind(df_brain_ratio, df_heart_ratio, df_kidney_ratio, df_liver_ratio)
df_all_use <- df_all %>% filter(drug %in% c("MET", "NR", "D+Q", "SPD"))
df_all_use <- df_all_use %>% filter((direction == "up" & mean_score < 0) | (direction == "down" & mean_score > 0))

### get pathway order
all_sort <- c()
for(i in unique(pathway_order$pathway_type)[c(2, 5, 1, 3, 4, 8:10)]){
  all_sort_i <- df_all_use %>% filter(pathway_type == i)
  all_sort <- rbind(all_sort, all_sort_i %>% group_by(Pathway_name) %>% mutate(drug_time = n_distinct(drug), organ_time = n_distinct(organ)) %>% arrange(desc(drug_time), desc(organ_time)) %>% as.data.frame())
}

df_all_use$Pathway_name <- factor(df_all_use$Pathway_name, levels = rev(unique(all_sort$Pathway_name)))


df_all_use$drug_organ <- paste0(tolower(df_all_use$drug), "|", df_all_use$organ)
df_all_use$drug_organ <- factor(df_all_use$drug_organ, levels = c("met|Brain", "met|Heart", "met|Kidney", "met|Liver",
                                                          "nr|Brain", "nr|Heart", "nr|Kidney", "nr|Liver",
                                                          "d+q|Brain", "d+q|Heart", "d+q|Kidney", "d+q|Liver",
                                                          "spd|Brain", "spd|Heart", "spd|Kidney", "spd|Liver"
))


### plot
p1 <- ggplot(data = df_all_use %>% filter(direction == "up" & mean_score < 0), aes(x = drug_organ, y = Pathway_name, size = ratio)) +
  geom_point(data = df_all_use %>% filter(direction == "up" & mean_score < -0.05), color = "#0051b0") +
  geom_point(data = df_all_use %>% filter(direction == "up" & mean_score >= -0.05 & mean_score < 0), aes(color = mean_score)) +
  scale_color_gradient2(high = "#ff6a6f", low = "#0051b0") +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), color = "grey", alpha = 0.75, linetype = 2) +
  geom_hline(yintercept = c(6.5, 12.5, 19.5, 34.5), color = "grey", alpha = 0.75, linetype = 2) +
  scale_size_continuous(range = c(2, 8)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 12, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")
        ) +
  scale_x_discrete(labels = c("met|Brain" = "Brain", "met|Heart" = "Heart", "met|Kidney" = "Kidney", "met|Liver" = "Liver",
                              "nr|Brain" = "Brain", "nr|Heart" = "Heart", "nr|Kidney" = "Kidney", "nr|Liver" = "Liver",
                              "d+q|Brain" = "Brain", "d+q|Heart" = "Heart", "d+q|Kidney" = "Kidney", "d+q|Liver" = "Liver",
                              "spd|Brain" = "Brain", "spd|Heart" = "Heart", "spd|Kidney" = "Kidney", "spd|Liver" = "Liver")) +
  labs(x = "", y = "", title = "")


p2 <- ggplot(df_all_use %>% filter(direction == "down" & mean_score > 0), aes(x = drug_organ, y = Pathway_name, color = mean_score, size = ratio)) +
  geom_point(data = df_all_use %>% filter(direction == "down" & mean_score > 0.05), color = "#ff6a6f") +
  geom_point(data = df_all_use %>% filter(direction == "down" & mean_score <= 0.05  & mean_score > 0), aes(color = mean_score)) +
  scale_color_gradient2(high = "#ff6a6f", low = "#0051b0") +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), color = "grey", alpha = 0.75, linetype = 2) +
  geom_hline(yintercept = c(6.5, 14.5), color = "grey", alpha = 0.75, linetype = 2) +
  scale_size_continuous(range = c(2, 8)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 12, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")
        ) +
  scale_x_discrete(labels = c("met|Brain" = "Brain", "met|Heart" = "Heart", "met|Kidney" = "Kidney", "met|Liver" = "Liver",
                              "nr|Brain" = "Brain", "nr|Heart" = "Heart", "nr|Kidney" = "Kidney", "nr|Liver" = "Liver",
                              "d+q|Brain" = "Brain", "d+q|Heart" = "Heart", "d+q|Kidney" = "Kidney", "d+q|Liver" = "Liver",
                              "spd|Brain" = "Brain", "spd|Heart" = "Heart", "spd|Kidney" = "Kidney", "spd|Liver" = "Liver")) +
  labs(x = "", y = "", title = "")
