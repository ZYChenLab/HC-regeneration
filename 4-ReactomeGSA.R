setwd("~/Desktop/WTvsDox-CellRanger-Outs")
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(MAST)
library(ReactomeGSA)
library(limma)
library(edgeR)

sample.integrated <- readRDS("~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Seurat_Combined.rds")

new.cluster.ids <- c("IdC", "CCOS","RMC","DC","CCOS","KO","CCOS","HeC","SC","RMC","TBRM","TBRM","DC","HC")
names(new.cluster.ids) <- levels(sample.integrated)
sample.integrated <- RenameIdents(sample.integrated, `0` = "IdC", `1` = "CCOS", `2` = "RMC", 
                                  `3` = "DC", `4` = "CCOS", `5` = "KO", `6` = "CCOS", `7` = "HeC", 
                                  `8` = "SC", `9` = "RMC", `10` = "TBRM", `11` = "TBRM", 
                                  `12` = "DC", `13` = "HC")

DefaultAssay(sample.integrated) <- "RNA"
sample.integrated <- NormalizeData(sample.integrated, verbose = FALSE)

sample.integrated$celltype.sample <- paste(Idents(sample.integrated), sample.integrated$sample, sep = "_")
sample.integrated$celltype <- Idents(sample.integrated)
Idents(sample.integrated) <- "celltype.sample"
levels(sample.integrated)

gsva_result <- analyse_sc_clusters(sample.integrated, verbose = TRUE)

pathway_expression <- pathways(gsva_result)

colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min

max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

Notch_path <- c("R-HSA-2122948", "R-HSA-2122947", "R-HSA-2979096", "R-HSA-2197563", "R-HSA-9013507", "R-HSA-9013508",
                "R-HSA-9013700", "R-HSA-9013695", "R-HSA-9017802", "R-HSA-350054", "R-HSA-157118", "R-HSA-1980143", 
                "R-HSA-2691230", "R-HSA-2894858", "R-HSA-2644602", "R-HSA-2644603", "R-HSA-2660825", "R-HSA-1980145", 
                "R-HSA-9012852", "R-HSA-9013694", "R-HSA-9017802", "R-HSA-8941856")     

PKA_path <- c("R-HSA-163615", "R-HSA-164378", "R-HSA-111931", "R-HSA-163358")    

mTOR_path <- c("R-HSA-380972", "R-HSA-9639288", "R-HSA-165159", "R-HSA-166208", "R-HSA-9639288")

WNT_path <- c("R-HSA-195721", "R-HSA-3238698", "R-HSA-5340573", "R-HSA-201688", "R-HSA-9673324",
              "R-HSA-5140745", "R-HSA-5099900", "R-HSA-3858494", "R-HSA-196299", "R-HSA-4641265",
              "R-HSA-8951430") 

#cAMP, Ca, CaM, Reelin, Connexin, ERK/AMP, Nodal, Hedge, Hippo, 
Signaling_path2 <- c("R-HSA-170660", "R-HSA-170670", "R-HSA-4086398", "R-HSA-111997", "R-HSA-8866376",
                     "R-HSA-190827", "R-HSA-198753", "R-HSA-1433617", "R-HSA-1181150", "R-HSA-5358351", "R-HSA-2028269")

Signaling_pathways <- c("R-HSA-1502540", "R-HSA-201451", "R-HSA-6802952", "R-HSA-177929", "R-HSA-1227986",
                        "R-HSA-1236394", "R-HSA-9006335","R-HSA-190236","R-HSA-5654736","R-HSA-5654738",
                        "R-HSA-5654741","R-HSA-5654743","R-HSA-372790","R-HSA-5358351","R-HSA-2028269",
                        "R-HSA-74752","R-HSA-449147","R-HSA-2586552","R-HSA-6806834","R-HSA-8852405",
                        "R-HSA-1181150","R-HSA-157118","R-HSA-1980143","R-HSA-1980145","R-HSA-9012852",
                        "R-HSA-9013694","R-HSA-187037","R-HSA-9006115","R-HSA-9034015","R-HSA-166520",
                        "R-HSA-9006927","R-HSA-9006931","R-HSA-186797","R-HSA-8848021","R-HSA-376176",
                        "R-HSA-9006934","R-HSA-5362517","R-HSA-194315","R-HSA-1433557","R-HSA-170834",
                        "R-HSA-9006936","R-HSA-2404192","R-HSA-194138","R-HSA-195721","R-HSA-983705",
                        "R-HSA-198765","R-HSA-187687","R-HSA-167044","R-HSA-198745","R-HSA-187706")

plot_gsva_heatmap(gsva_result, 
                  pathway_ids = PKA_path, #(or Notch, mTOR, Signaling_paths etc)
                  dendrogram = "col",
                  scale = "row",
                  key = F,
                  lwid=c(0.1,4),
                  truncate_names = F)