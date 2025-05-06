### Figure 8
### B

### load rev- and pro-aging DEGs data
load("brain_rev_pro.RData")
load("heart_rev_pro_EC.RData")
load("kidney_rev_pro.RData")
load("liver_rev_pro.RData")

### prepare data
### EC, Mac1, PT, SMC, Fib
data_brain$origin <- "Brain"
data_heart$origin <- "Heart"
data_kidney$origin <- "Kidney"
data_liver$origin <- "Liver"

data_all_organ <- rbind(data_brain, data_heart, data_kidney, data_liver)
data_all_organ$ratio <- data_all_organ$value / data_all_organ$sum
data_all_organ$sum[which(is.na(data_all_organ$sum))] <- 0
data_all_organ$ratio[which(is.na(data_all_organ$ratio))] <- 0
colnames(data_all_organ) <- c("celltype", "group", "count", "label", "sum", "origin", "ratio")

### get ratio
get_matrix_ratio <- function(origin_label){
  data_all_organ_i <- data_all_organ %>% filter(group == origin_label)
  i_matrix <- matrix(0, 8, 48)
  rownames(i_matrix) <- c("Brain-pro", "Heart-pro", "Kidney-pro", "Liver-pro", "Brain-rev", "Heart-rev", "Kidney-rev", "Liver-rev")
  colnames(i_matrix) <- c('EC', 'Mac1', 'PC', 'SMC', 'Fib',
                          levels(brain.harmony@meta.data$celltype)[c(1:4, 6:17, 19)],
                          levels(heart.harmony@meta.data$celltype)[c(2, 5:10, 12:13)],
                          levels(kidney.harmony@meta.data$celltype)[c(4:14)],
                          levels(liver.harmony@meta.data$celltype)[c(1, 3:7)]
  )

  for(i in 1:nrow(data_all_organ_i)){
    if(data_all_organ_i$label[i] == "Pro-aging"){
      i_matrix[paste0(data_all_organ_i$origin[i], "-pro"), as.character(data_all_organ_i$celltype[i])] <- data_all_organ_i$ratio[i]
    } else {
      i_matrix[paste0(data_all_organ_i$origin[i], "-rev"), as.character(data_all_organ_i$celltype[i])] <- data_all_organ_i$ratio[i]
    }
  }

  i_matrix <- i_matrix[, names(sort(colSums(i_matrix), decreasing = T))]
  i_df <- melt(i_matrix)
  i_df$group <- origin_label
  return(i_df)
}

met_df <- get_matrix_ratio("MET")
nr_df <- get_matrix_ratio("NR")
dq_df <- get_matrix_ratio("D+Q")
spd_df <- get_matrix_ratio("SPD")


bar_df <- rbind(met_df, nr_df, dq_df, spd_df)
colnames(bar_df) <- c("label", "celltype", "ratio", "group")
bar_df$origin <- apply(bar_df, 1, function(x){
  strsplit(as.character(x[1]), split = "-")[[1]][1]
})

bar_df$type <- apply(bar_df, 1, function(x){
  strsplit(as.character(x[1]), split = "-")[[1]][2]
})


### for commom cell
color_origin <- c("#d8c91c", "#47a3bc", "#a3c8dc", "#a91f24", "#f09594",
                  "#8e4b98", "#fdbd10", "#0066b2")
origin_order <- c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
col_origin <- setNames(color_origin, origin_order)
### for barplot
plots <- list()
color_rev_pro <- data.frame(color = c("#326582", "#5F1115", "#B81452", "#64346A", as.character(col_origin[3:6])),
                            label = c("MET", "NR", "D+Q", "SPD", "MET", "NR", "D+Q", "SPD"))

bar_df_common <- bar_df
bar_df_common_i <- bar_df_common %>% filter(celltype == "EC")

bar_df_common_i$origin <- factor(bar_df_common_i$origin, levels = c((bar_df_common_i %>% group_by(origin) %>% reframe(sum = sum(ratio)) %>% arrange(-sum))$origin))
bar_df_common_i$group <- factor(bar_df_common_i$group, levels = c("MET", "NR", "D+Q", "SPD"))
for(j in c("MET", "NR", "D+Q", "SPD")){
    color_rev_pro_j <- color_rev_pro %>% filter(label == j)
    color_rev_pro_j <- as.character(color_rev_pro_j$color)
    names(color_rev_pro_j) <- c("pro", "rev")

    bar_df_common_i_j <- bar_df_common_i %>% filter(group == j)
    plots[[paste0(i, "_", j)]] <- ggplot() +
      geom_col(data = bar_df_common_i_j, aes(x = origin, y = 1), width = 0.75, color = NA, fill = "grey", alpha = 0.5) +
      geom_col(data = bar_df_common_i_j, aes(x = origin, y = ratio, fill = type), width = 0.75, color = NA) +
      theme_classic() +
      scale_fill_manual(values = color_rev_pro_j) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title = element_text(size = 12)) +
      labs(x = "", y = "DEG ratio")
  }


for(i in c("Fib", "Mac1", "SMC", "PC")){
  bar_df_common_i <- bar_df_common %>% filter(celltype == i)
  bar_df_common_i$origin <- factor(bar_df_common_i$origin, levels = c((bar_df_common_i %>% group_by(origin) %>% reframe(sum = sum(ratio)) %>% arrange(-sum))$origin))
  bar_df_common_i$group <- factor(bar_df_common_i$group, levels = c("MET", "NR", "D+Q", "SPD"))

  for(j in c("MET", "NR", "D+Q", "SPD")){
    color_rev_pro_j <- color_rev_pro %>% filter(label == j)
    color_rev_pro_j <- as.character(color_rev_pro_j$color)
    names(color_rev_pro_j) <- c("pro", "rev")

    bar_df_common_i_j <- bar_df_common_i %>% filter(group == j)
    plots[[paste0(i, "_", j)]] <- ggplot() +
      geom_col(data = bar_df_common_i_j, aes(x = origin, y = 1), width = 0.75, color = NA, fill = "grey", alpha = 0.5) +
      geom_col(data = bar_df_common_i_j, aes(x = origin, y = ratio, fill = type), width = 0.75, color = NA) +
      theme_classic() +
      scale_fill_manual(values = color_rev_pro_j) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
            axis.title = element_text(size = 12),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank()
      ) +
      labs(x = "", y = "")
  }
}


plots[[1]] + plots[[5]] + plots[[9]] + plots[[13]] + plots[[17]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "MET_rev_pro_barplot_ratio.pdf", width = 8, height = 2)
plots[[2]] + plots[[6]] + plots[[10]] + plots[[14]] + plots[[18]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "NR_rev_pro_barplot_ratio.pdf", width = 8, height = 2)
plots[[3]] + plots[[7]] + plots[[11]] + plots[[15]] + plots[[19]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "DQ_rev_pro_barplot_ratio.pdf", width = 8, height = 2)
plots[[4]] + plots[[8]] + plots[[12]] + plots[[16]] + plots[[20]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "SPD_rev_pro_barplot_ratio.pdf", width = 8, height = 2)








### get number
get_matrix_number <- function(origin_label){
  data_all_organ_i <- data_all_organ %>% filter(group == origin_label)
  i_matrix <- matrix(0, 12, 48)
  rownames(i_matrix) <- c("Brain-pro", "Heart-pro", "Kidney-pro", "Liver-pro", "Brain-rev", "Heart-rev", "Kidney-rev", "Liver-rev", "Brain-aging", "Heart-aging", "Kidney-aging", "Liver-aging")
  colnames(i_matrix) <- c('EC', 'Mac1', 'PC', 'SMC', 'Fib',
                          levels(brain.harmony@meta.data$celltype)[c(1:4, 6:17, 19)],
                          levels(heart.harmony@meta.data$celltype)[c(2, 5:10, 12:13)],
                          levels(kidney.harmony@meta.data$celltype)[c(4:14)],
                          levels(liver.harmony@meta.data$celltype)[c(1, 3:7)]
  )

  for(i in 1:nrow(data_all_organ_i)){
    i_matrix[paste0(data_all_organ_i$origin[i], "-aging"), as.character(data_all_organ_i$celltype[i])] <- data_all_organ_i$sum[i]
    if(data_all_organ_i$label[i] == "Pro-aging"){
      i_matrix[paste0(data_all_organ_i$origin[i], "-pro"), as.character(data_all_organ_i$celltype[i])] <- data_all_organ_i$count[i]
    } else {
      i_matrix[paste0(data_all_organ_i$origin[i], "-rev"), as.character(data_all_organ_i$celltype[i])] <- data_all_organ_i$count[i]
    }
  }

  i_matrix <- i_matrix[, names(sort(colSums(i_matrix), decreasing = T))]
  i_df <- melt(i_matrix)
  i_df$group <- origin_label
  return(i_df)
}

met_df <- get_matrix_number("MET")
nr_df <- get_matrix_number("NR")
dq_df <- get_matrix_number("D+Q")
spd_df <- get_matrix_number("SPD")

bar_df <- rbind(met_df, nr_df, dq_df, spd_df)
colnames(bar_df) <- c("label", "celltype", "count", "group")
bar_df$origin <- apply(bar_df, 1, function(x){
  strsplit(as.character(x[1]), split = "-")[[1]][1]
})

bar_df$type <- apply(bar_df, 1, function(x){
  strsplit(as.character(x[1]), split = "-")[[1]][2]
})


### for commom cell
color_origin <- c("#d8c91c", "#47a3bc", "#a3c8dc", "#a91f24", "#f09594",
                  "#8e4b98", "#fdbd10", "#0066b2")
origin_order <- c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
col_origin <- setNames(color_origin, origin_order)

### for barplot
plots <- list()
color_rev_pro <- data.frame(color = c("#326582", "#5F1115", "#B81452", "#64346A", as.character(col_origin[3:6])),
                            label = c("MET", "NR", "D+Q", "SPD", "MET", "NR", "D+Q", "SPD"))

bar_df_common <- bar_df %>% filter(celltype %in% c("EC", "Fib", "Mac1", "SMC", "PC"))
max_y <- 170
bar_df_common_i <- bar_df_common %>% filter(celltype == "EC")

bar_df_common_i$origin <- factor(bar_df_common_i$origin, levels = c((bar_df_common_i %>% group_by(origin) %>% reframe(sum = sum(count)) %>% arrange(-sum))$origin))
bar_df_common_i$group <- factor(bar_df_common_i$group, levels = c("MET", "NR", "D+Q", "SPD"))
for(j in c("MET", "NR", "D+Q", "SPD")){
    color_rev_pro_j <- color_rev_pro %>% filter(label == j)
    color_rev_pro_j <- c(as.character(color_rev_pro_j$color))
    names(color_rev_pro_j) <- c("pro", "rev")

    bar_df_common_i_j <- bar_df_common_i %>% filter(group == j)
    plots[[paste0("EC_", j)]] <- ggplot() +
      geom_col(data = bar_df_common_i_j %>% filter(!type == "aging"), aes(x = origin, y = count, fill = type), width = 0.75, color = NA) +
      theme_classic() +
      scale_fill_manual(values = color_rev_pro_j) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max_y)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title = element_text(size = 12)) +
      labs(x = "", y = "DEG number")
  }


for(i in c("Fib", "PC", "Mac1", "SMC")){
  bar_df_common_i <- bar_df_common %>% filter(celltype == i)
  bar_df_common_i$origin <- factor(bar_df_common_i$origin, levels = c((bar_df_common_i %>% group_by(origin) %>% reframe(sum = sum(count)) %>% arrange(-sum))$origin))
  bar_df_common_i$group <- factor(bar_df_common_i$group, levels = c("MET", "NR", "D+Q", "SPD"))

  for(j in c("MET", "NR", "D+Q", "SPD")){
    color_rev_pro_j <- color_rev_pro %>% filter(label == j)
    color_rev_pro_j <- c(as.character(color_rev_pro_j$color))
    names(color_rev_pro_j) <- c("pro", "rev")

    bar_df_common_i_j <- bar_df_common_i %>% filter(group == j)
    plots[[paste0(i, "_", j)]] <- ggplot() +
      geom_col(data = bar_df_common_i_j %>% filter(!type == "aging"), aes(x = origin, y = count, fill = type), width = 0.75, color = NA) +
      theme_classic() +
      scale_fill_manual(values = color_rev_pro_j) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max_y + 2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
            axis.title = element_text(size = 12),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank()
      ) +
      labs(x = "", y = "")
  }
}


plots[[1]] + plots[[5]] + plots[[9]] + plots[[13]] + plots[[17]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "MET_rev_pro_barplot_number.pdf", width = 8, height = 2)
plots[[2]] + plots[[6]] + plots[[10]] + plots[[14]] + plots[[18]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "NR_rev_pro_barplot_number.pdf", width = 8, height = 2)
plots[[3]] + plots[[7]] + plots[[11]] + plots[[15]] + plots[[19]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "DQ_rev_pro_barplot_number.pdf", width = 8, height = 2)
plots[[4]] + plots[[8]] + plots[[12]] + plots[[16]] + plots[[20]] + plot_layout(nrow = 1, guides = "collect", width = c(5, 4, 4, 4, 4))
ggsave(filename = "SPD_rev_pro_barplot_number.pdf", width = 8, height = 2)
