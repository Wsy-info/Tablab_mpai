### Figure 7
### H

library(fgsea)
library(dplyr)
library(ggplot2)

### load KEGG information
kegg_all_pathway <- read.csv("keggAllCompoundReshapedData2-2024-11-02.csv")
kegg_all_pathway$pathwayName[which(kegg_all_pathway$pathwayName == "Glycolysis / Gluconeogenesis")] <- "Glycolysis or Gluconeogenesis"

### prepare data
kegg_all_pathway_new <- c()
for(i in 1:12164){
  compound <- strsplit(kegg_all_pathway[i, 2], split = "\\; ")[[1]]
  df <- data.frame(CID = kegg_all_pathway[i, 1],
                   CompoundName = compound,
                   pathwayId = kegg_all_pathway[i, 3],
                   pathwayName = kegg_all_pathway[i, 4]
  )
  kegg_all_pathway_new <- rbind(kegg_all_pathway_new, df)
}

### load metabolomics data
expression_metabolism <- read.table("metabolites_full_table.txt",
                                    sep = "\t", quote = "", comment.char = "", header = T, row.names = 1)

### get all pathways
met_up <- read.csv("met_up.csv") %>% filter(Hits > 1 & Impact > 0)
nr_up <- read.csv("nr_up.csv") %>% filter(Hits > 1 & Impact > 0)
dq_up <- read.csv("dq_up.csv") %>% filter(Hits > 1 & Impact > 0)
spd_up <- read.csv("spd_up.csv") %>% filter(Hits > 1 & Impact > 0)

pathway_rev_up <- rbind(data.frame(met_up, group = "MET"),
                        data.frame(nr_up, group = "NR"),
                        data.frame(dq_up, group = "D+Q"),
                        data.frame(spd_up, group = "SPD")
)

### down
met_down <- read.csv("met_down.csv") %>% filter(Hits > 1 & Impact > 0)
nr_down <- read.csv("nr_down.csv") %>% filter(Hits > 1 & Impact > 0)
dq_down <- read.csv("dq_down.csv") %>% filter(Hits > 1 & Impact > 0)
spd_down <- read.csv("spd_down.csv") %>% filter(Hits > 1 & Impact > 0)

pathway_rev_down <- rbind(data.frame(met_down, group = "MET"),
                        data.frame(nr_down, group = "NR"),
                        data.frame(dq_down, group = "D+Q"),
                        data.frame(spd_down, group = "SPD")
)

all_pathway <- union(as.character(pathway_rev_up$X), as.character(pathway_rev_down$X))


### MSEA
all_fgsea_result_group <- c()

### Old vs. Young
FCgenelist <- expression_metabolism$C_vs_F_log2FC
FCgenelist <- - FCgenelist * -1
names(FCgenelist) <- expression_metabolism$name
FCgenelist <- sort(FCgenelist, decreasing = T)
all_fgsea_result <- c()
for(i in 1:(length(all_pathway))){
  pathway_i <- all_pathway[i]
  pathway_list_i <- list(i = kegg_all_pathway_new[which(kegg_all_pathway_new$pathwayName == pathway_i), 2])
  names(pathway_list_i) <- pathway_i
  fgsea_result <- fgsea(pathways = pathway_list_i,
                        stats = FCgenelist)
  all_fgsea_result <- rbind(all_fgsea_result, fgsea_result)
}
all_fgsea_result$group <- "Young"
all_fgsea_result_group <- rbind(all_fgsea_result_group, all_fgsea_result)

### MET vs. Old
FCgenelist <- expression_metabolism$A_vs_C_log2FC
names(FCgenelist) <- expression_metabolism$name
FCgenelist <- sort(FCgenelist, decreasing = T)
all_fgsea_result <- c()
for(i in 1:(length(all_pathway))){
  pathway_i <- all_pathway[i]
  pathway_list_i <- list(i = kegg_all_pathway_new[which(kegg_all_pathway_new$pathwayName == pathway_i), 2])
  names(pathway_list_i) <- pathway_i
  fgsea_result <- fgsea(pathways = pathway_list_i,
                        stats = FCgenelist)
  all_fgsea_result <- rbind(all_fgsea_result, fgsea_result)
}
all_fgsea_result$group <- "MET"
all_fgsea_result_group <- rbind(all_fgsea_result_group, all_fgsea_result)

### NR vs. Old
FCgenelist <- expression_metabolism$B_vs_C_log2FC
names(FCgenelist) <- expression_metabolism$name
FCgenelist <- sort(FCgenelist, decreasing = T)
all_fgsea_result <- c()
for(i in 1:(length(all_pathway))){
  pathway_i <- all_pathway[i]
  pathway_list_i <- list(i = kegg_all_pathway_new[which(kegg_all_pathway_new$pathwayName == pathway_i), 2])
  names(pathway_list_i) <- pathway_i
  fgsea_result <- fgsea(pathways = pathway_list_i,
                        stats = FCgenelist)
  all_fgsea_result <- rbind(all_fgsea_result, fgsea_result)
}
all_fgsea_result$group <- "NR"
all_fgsea_result_group <- rbind(all_fgsea_result_group, all_fgsea_result)

### D+Q vs. Old
FCgenelist <- expression_metabolism$C_vs_D_log2FC
names(FCgenelist) <- expression_metabolism$name
FCgenelist <- sort(FCgenelist, decreasing = T)
all_fgsea_result <- c()
for(i in 1:(length(all_pathway))){
  pathway_i <- all_pathway[i]
  pathway_list_i <- list(i = kegg_all_pathway_new[which(kegg_all_pathway_new$pathwayName == pathway_i), 2])
  names(pathway_list_i) <- pathway_i
  fgsea_result <- fgsea(pathways = pathway_list_i,
                        stats = FCgenelist)
  all_fgsea_result <- rbind(all_fgsea_result, fgsea_result)
}
all_fgsea_result$group <- "D+Q"
all_fgsea_result_group <- rbind(all_fgsea_result_group, all_fgsea_result)

### SPD vs. Old
FCgenelist <- expression_metabolism$C_vs_E_log2FC
names(FCgenelist) <- expression_metabolism$name
FCgenelist <- sort(FCgenelist, decreasing = T)
all_fgsea_result <- c()
for(i in 1:(length(all_pathway))){
  pathway_i <- all_pathway[i]
  pathway_list_i <- list(i = kegg_all_pathway_new[which(kegg_all_pathway_new$pathwayName == pathway_i), 2])
  names(pathway_list_i) <- pathway_i
  fgsea_result <- fgsea(pathways = pathway_list_i,
                        stats = FCgenelist)
  all_fgsea_result <- rbind(all_fgsea_result, fgsea_result)
}
all_fgsea_result$group <- "SPD"
all_fgsea_result_group <- rbind(all_fgsea_result_group, all_fgsea_result)

all_fgsea_result_group$group <- as.character(all_fgsea_result_group$group)
gsea_matrix <- reshape2::dcast(all_fgsea_result_group, pathway ~ group, value.var = "NES")
rownames(gsea_matrix) <- gsea_matrix$pathway

gsea_matrix <- gsea_matrix[, -1]


### filter pathway
GSEA_pathway_index <- c("Purine metabolism", "Propanoate metabolism", "Arginine and proline metabolism", "Amino sugar and nucleotide sugar metabolism", ### up
                        "Metabolism of xenobiotics by cytochrome P450", "Glycine, serine and threonine metabolism", "Retinol metabolism" ### down
)
all_fgsea_result_group_draw <- all_fgsea_result_group %>% filter(pathway %in% GSEA_pathway_index)
all_fgsea_result_group_draw$pathway <- factor(all_fgsea_result_group_draw$pathway, levels = rev(GSEA_pathway_index))
all_fgsea_result_group_draw$group <- factor(all_fgsea_result_group_draw$group, levels = c("Young", "MET", "NR", "D+Q", "SPD"), labels = c("Old/Young", "MET/Old", "NR/Old", "D+Q/Old", "SPD/Old"))

### dotplot
ggplot(all_fgsea_result_group_draw, aes(x = group, y = pathway)) +
  geom_point(aes(color = NES, size = abs(NES))) +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_gradient2(low = "#0051b0", high = "#ff6a6f", mid = "white") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 11, color = "black")) +
  labs(x = "", y = "")
ggsave("GSEA_result_dotplot_new.pdf", width = 7, height = 3)