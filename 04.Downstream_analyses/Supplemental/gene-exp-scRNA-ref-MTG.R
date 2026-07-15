#Analysis of neuronal gene-expression profiles from the MTG of 166,868 nuclei from 5 donors from SEA-AD study (2022)
# https://brain-map.org/consortia/sea-ad/human-mtg-10x-sea-ad

library(Seurat)
library(dplyr)
seu <- readRDS("/home/AD/vbidhan/study232-missionbio_TDP-C/Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds")
seu_neurons <- subset(seu, subset = class_label %in% c("Neuronal: GABAergic", "Neuronal: Glutamatergic"))

genes <- c("TARDBP", "UNC13A", "TET2", "OPTN", "TBK1", "TMEM106B", "GRN") 
genes <- genes[genes %in% rownames(seu_neurons)]
avg_exp <- rowMeans(GetAssayData(seu_neurons, assay = "RNA", layer = "data")[genes, ])

# Total number of genes
print(nrow(seu_neurons))
print(ncol(seu_neurons))
#Calculate deciles for all genes
avg_all <- rowMeans(GetAssayData(seu_neurons, assay = "RNA", layer = "data"))
avg_df <- data.frame(gene = names(avg_all),avg_exp = avg_all) %>% mutate(decile = ntile(avg_exp, 10))
avg_df[genes, ]