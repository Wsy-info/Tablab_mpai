### Figure 6
### G

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

### use data from Fig6/6E_F.R
### find downregulated genes
met_down <- c()
nr_down <- c()
dq_down <- c()
for(i in 1:length(list1)){
  met_down <- c(met_down, list1[[i]]$`Young_Old(MET)`)
  nr_down <- c(nr_down, list1[[i]]$`Young_Old(NR)`)
  dq_down <- c(dq_down, list1[[i]]$`Young_Old(D+Q)`)
}

met_down <- unique(met_down)
nr_down <- unique(nr_down)
dq_down <- unique(dq_down)

### find upregulated genes
met_up <- c()
nr_up <- c()
dq_up <- c()
for(i in 1:length(list1)){
  met_up <- c(met_up, list2[[i]]$Old_MET)
  nr_up <- c(nr_up, list2[[i]]$Old_NR)
  dq_up <- c(dq_up, list2[[i]]$`Old_D+Q`)
}

met_up <- unique(met_up)
nr_up <- unique(nr_up)
dq_up <- unique(dq_up)

get_pathway <- function(genes, filename){
  entrez_ids <- bitr(genes, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
  go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
  met_result <- simplify(go_results)
  met_result <- met_result@result %>% filter(pvalue < 0.05 & Count > 2) %>% arrange(pvalue)

  ### prepare data
  met_result$geneID <- gsub("\\/", ",", met_result$geneID)
  met_result <- data.frame(met_result[, c(1, 2, 5, 6)],
                           Phenotype = 1,
                           met_result[, 8])
  write.table(met_result, file = filename, sep = "\t", quote = F, row.names = F)
  return(met_result)
}


### save results for network(Cytoscape_v3.10.2)
met_up_resluts <- get_pathway(met_up, "met_up.txt")
nr_up_resluts <- get_pathway(nr_up, "nr_up.txt")
dq_up_resluts <- get_pathway(dq_up, "dq_up.txt")

met_down_resluts <- get_pathway(met_down, "met_down.txt")
nr_down_resluts <- get_pathway(nr_down, "nr_down.txt")
dq_down_resluts <- get_pathway(dq_down, "dq_down.txt")