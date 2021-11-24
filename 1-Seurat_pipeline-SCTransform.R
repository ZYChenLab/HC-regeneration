setwd("~/Desktop/WTvsDox-CellRanger-Outs")
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(MAST)

###
Dox1.data <- Read10X(data.dir = "~/Desktop/WTvsDox-CellRanger-Outs/A1-filtered_feature_bc_matrix")
Control.data <- Read10X(data.dir = "~/Desktop/WTvsDox-CellRanger-Outs/A2-filtered_feature_bc_matrix")
Dox2.data <- Read10X(data.dir = "~/Desktop/WTvsDox-CellRanger-Outs/A3-filtered_feature_bc_matrix")

Dox1 <- CreateSeuratObject(counts = Dox1.data, project = "Dox1", min.cells = 3, min.features = 200)
Dox1[["percent.mt"]] <- PercentageFeatureSet(Dox1, pattern = "^mt-")
VlnPlot(Dox1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Dox1 <- subset(Dox1, subset = nFeature_RNA > 200 & percent.mt < 25)

Control <- CreateSeuratObject(counts = Control.data, project = "Control", min.cells = 3, min.features = 200)
Control[["percent.mt"]] <- PercentageFeatureSet(Control, pattern = "^mt-")
VlnPlot(Control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Control <- subset(Control, subset = nFeature_RNA > 200 & percent.mt < 25)

Dox2 <- CreateSeuratObject(counts = Dox2.data, project = "Dox2", min.cells = 3, min.features = 200)
Dox2[["percent.mt"]] <- PercentageFeatureSet(Dox2, pattern = "^mt-")
VlnPlot(Dox2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Dox2 <- subset(Dox2, subset = nFeature_RNA > 200 & percent.mt < 25)

Dox1$stim <- "Dox"
Control$stim <- "Control"
Dox2$stim <- "Dox"
Dox1$sample <- "Dox1"
Control$sample <- "Control"
Dox2$sample <- "Dox2"

sample.integrated <- merge(x = Control, y = list(Dox1, Dox2), add.cell.ids = c("Control", "Dox1", "Dox2"), project = "MyProject")

all.list <- SplitObject(sample.integrated, split.by = "sample")

for (i in names(all.list)) {
  all.list[[i]] <- SCTransform(all.list[[i]], verbose = FALSE)
}

sample.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 4000)
sample.list <- PrepSCTIntegration(object.list = all.list, anchor.features = sample.features)
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT", 
                                         anchor.features = sample.features)
sample.integrated <- IntegrateData(anchorset = sample.anchors, normalization.method = "SCT")
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE)
sample.integrated <- RunUMAP(object = sample.integrated, dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:20)
sample.integrated <- FindClusters(sample.integrated, resolution = 0.5)

saveRDS(sample.integrated, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Seurat_Combined.rds")

p1 <- DimPlot(sample.integrated, reduction = "umap", group.by = "orig.ident", pt.size = 0.2)
p2 <- DimPlot(sample.integrated, reduction = "umap", label = T,   label.size = 6, repel = T, pt.size = 0.2)
plot_grid(p1, p2)
DimPlot(sample.integrated, reduction = "umap", split.by = "stim")

all.markers <- FindAllMarkers(sample.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(sample.integrated, features = top10$gene) + NoLegend()

write.csv(all.markers, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/FindAllMarkers.csv")

VlnPlot(sample.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "stim")
VlnPlot(sample.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")

DefaultAssay(sample.integrated) <- "RNA"
cons0 <- FindConservedMarkers(sample.integrated, ident.1 = 0, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons1 <- FindConservedMarkers(sample.integrated, ident.1 = 1, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons2 <- FindConservedMarkers(sample.integrated, ident.1 = 2, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons3 <- FindConservedMarkers(sample.integrated, ident.1 = 3, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons4 <- FindConservedMarkers(sample.integrated, ident.1 = 4, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons5 <- FindConservedMarkers(sample.integrated, ident.1 = 5, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons6 <- FindConservedMarkers(sample.integrated, ident.1 = 6, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons7 <- FindConservedMarkers(sample.integrated, ident.1 = 7, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons8 <- FindConservedMarkers(sample.integrated, ident.1 = 8, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons9 <- FindConservedMarkers(sample.integrated, ident.1 = 9, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons10 <- FindConservedMarkers(sample.integrated, ident.1 = 10, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons11 <- FindConservedMarkers(sample.integrated, ident.1 = 11, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")
cons12 <- FindConservedMarkers(sample.integrated, ident.1 = 12, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "stim")

write.csv(cons0, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons0.csv")
write.csv(cons1, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons1.csv")
write.csv(cons2, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons2.csv")
write.csv(cons3, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons3.csv")
write.csv(cons4, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons4.csv")
write.csv(cons5, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons5.csv")
write.csv(cons6, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons6.csv")
write.csv(cons7, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons7.csv")
write.csv(cons8, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons8.csv")
write.csv(cons9, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons9.csv")
write.csv(cons10, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons10.csv")
write.csv(cons11, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons11.csv")
write.csv(cons12, file = "~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Cons12.csv")

new.cluster.ids <- c("IdC", "CCOS","RMC","DC","CCOS","KO","CCOS","HeC","SC","RMC","TBRM","TBRM","DC","HC")
names(new.cluster.ids) <- levels(sample.integrated)
sample.integrated <- RenameIdents(sample.integrated, `0` = "IdC", `1` = "CCOS", `2` = "RMC", 
                                  `3` = "DC", `4` = "CCOS", `5` = "KO", `6` = "CCOS", `7` = "HeC", 
                                  `8` = "SC", `9` = "RMC", `10` = "TBRM", `11` = "TBRM", 
                                  `12` = "DC", `13` = "HC")
