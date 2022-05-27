sample.integrated <- readRDS("~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Seurat_Combined.rds")

new.cluster.ids <- c("IdC", "CCOS","RMC","DC","CCOS","KO","CCOS","HeC","uSC","RMC","MLC","MLC","DC","HC")
names(new.cluster.ids) <- levels(sample.integrated)
sample.integrated <- RenameIdents(sample.integrated, `0` = "IdC", `1` = "CCOS", `2` = "RMC", 
                                  `3` = "DC", `4` = "CCOS", `5` = "KO", `6` = "CCOS", `7` = "HeC", 
                                  `8` = "uSC", `9` = "RMC", `10` = "MLC", `11` = "MLC", 
                                  `12` = "DC", `13` = "HC")

#DefaultAssay(sample.integrated) <- "RNA"
#sample.integrated <- NormalizeData(sample.integrated, verbose = FALSE)

# Notch1 and Myc
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Notch1", "Myc"),
            label = F,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = TRUE,
            pt.size = 0.2,
            split.by = "stim",
)

VlnPlot(sample.integrated, features = c("Notch1", "Myc"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, split.by = "stim")


#IdC Interdental cells
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Galm", "Nfix", "Epcam", "Nfib"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Galm", "Nfix", "Epcam", "Nfib"), 
        pt.size = 0, assay = "RNA", adjust = 1, stack = T, flip = T, fill.by = "ident")

#CCOS Claudius cells/outer sulcus cells
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Cp", "Col3a1", "Fstl1", "Lhfp"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Cp", "Col3a1", "Fstl1", "Lhfp"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, fill.by = "ident")



#RMC Reisner's membrane cells
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Vmo1", "Stim2", "Meis2", "Fut9"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Vmo1", "Stim2", "Meis2", "Fut9"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, fill.by = "ident")

#DC Deiter's cells
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Otog", "Slitrk6", "Emid1", "Nupr1l"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Otog", "Slitrk6", "Emid1", "Nupr1l"), 
         pt.size = 0, assay = "integrated", adjust = 1, log = F, stack = T, flip = T, y.max = 2, fill.by = "ident")


#KO Kolliker's organ
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Epyc", "Clu", "Slc39a8", "Lum"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Epyc", "Clu", "Slc39a8", "Lum"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, fill.by = "ident")

#HeC Hensen's cells
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Plp1", "Pmp22", "Sostdc1", "Art3"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Plp1", "Pmp22", "Sostdc1", "Art3"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, fill.by = "ident")

#Unclassified supporting cells
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Ciart", "Hspa1a", "Hspa1b", "Trp53rkb"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Ciart", "Hspa1a", "Hspa1b", "Trp53rkb"), 
        pt.size = 0, assay = "integrated", adjust = 2, stack = T, flip = T, fill.by = "ident")

#TBRM Tympanic border resident macrophages
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Cx3cr1", "Emilin2", "Lyz2", "Fcgr1"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Cx3cr1", "Emilin2", "Lyz2", "Fcgr1"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, fill.by = "ident")


#HC Hair cells
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Pvalb", "Atp8a2", "Faim2", "Dnm3"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6
)


VlnPlot(sample.integrated, features = c("Pvalb", "Atp8a2", "Faim2", "Dnm3", "Nefl", "Ablim2",
                                        "Dner", "Calm2", "Calb2", "Ormdl3", "Apba1", "Espn",
                                        "Slc26a5", "Myo7a"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, fill.by = "ident")


#Mitotic markers
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
cochlea <- CellCycleScoring(sample.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(cochlea[[]])
p3 <- DimPlot(cochlea)
p4 <- DimPlot(cochlea, split.by = "sample")
plot_grid(p3, p4)

#G2M genes
VlnPlot(cochlea, features = c("Hmgb2", "Cdk1", "Ube2c", "Birc5", "Tpx2", "Top2a",
                                        "Mki67", "Cenpf", "Aurkb", "Cenpa"), 
        pt.size = 0, assay = "RNA", adjust = 1, stack = T, flip = T, fill.by = "ident")

VlnPlot(sample.integrated, features = c("Hmgb2", "Cdk1", "Ube2c", "Birc5", "Tpx2", "Top2a",
                              "Mki67", "Cenpf", "Aurkb", "Cenpa"), 
        pt.size = 0, assay = "RNA", adjust = 1, stack = T, flip = T, fill.by = "ident")

#S genes
VlnPlot(cochlea, features = c("Pcna", "Mcm4", "Rrm1", "Ung", "Dtl", "Hells",
                                        "Ubr7", "Rad51", "Usp1", "Pola1", "Dscc1", "E2f8"), 
        pt.size = 0, assay = "RNA", adjust = 1, stack = T, flip = T, fill.by = "ident")

VlnPlot(sample.integrated, features = c("Pcna", "Mcm4", "Rrm1", "Ung", "Dtl", "Hells",
                                        "Ubr7", "Rad51", "Usp1", "Pola1", "Dscc1", "E2f8"), 
        pt.size = 0, assay = "RNA", adjust = 1, stack = T, flip = T, fill.by = "ident")

#Mitotic G2M + M genes in clusters
VlnPlot(sample.integrated, features = c( "Ube2c", "Birc5", "Tpx2", "Top2a", 
                                        "Rad51", "Pcna", "Dscc1", "Mcm4"), 
        pt.size = 0, assay = "RNA", adjust = 1, stack = T, flip = T, fill.by = "ident",
        idents = c("IdC", "RMC"))

#Heatmap from cluster averages
DoHeatmap(cluster.averages, features = c("Galm", "Nfix", "Epcam", "Tspan8",
                                         "Cp", "Col3a1", "Fstl1", "Lhfp",
                                         "Vmo1", "Stim2", "Meis2", "Fut9",
                                         "S100a1", "Cep41", "Rbp7", "Ceacam16",
                                         "Epyc", "Clu", "Slc39a8", "Lum",
                                         "Plp1", "Pmp22", "Sostdc1", "Art3",
                                         "Nefm", "Stmn3", "Hspa1b", "Kcnn2",
                                         "Cx3cr1", "Emilin2", "Lyz2", "Fcgr1",
                                         "Pvalb", "Atp8a2", "Faim2", "Dnm3"), size = 4, 
          disp.max = 6, draw.lines = FALSE)

#Heatmap from sample integrated with integrated assay
DoHeatmap(sample.integrated, features = c("Galm", "Nfix", "Epcam", "Tspan8",
                                          "Cp", "Col3a1", "Fstl1", "Lhfp",
                                          "Vmo1", "Stim2", "Meis2", "Fut9",
                                          "S100a1", "Cep41", "Rbp7", "Ceacam16",
                                          "Epyc", "Clu", "Slc39a8", "Lum",
                                          "Plp1", "Pmp22", "Sostdc1", "Art3",
                                          "Nefm", "Stmn3", "Hspa1b", "Kcnn2",
                                          "Cx3cr1", "Emilin2", "Lyz2", "Fcgr1",
                                          "Pvalb", "Atp8a2", "Faim2", "Dnm3"), 
          assay = "integrated", size = 4, disp.max = 6, draw.lines = FALSE)

#Expression analysis of RT-qPCR targets
VlnPlot(sample.integrated, features = c("Notch1", "Six1", "Gata3", "Eya1", "Hes5", "Hes1", "Cdkn1b", "Myc", "Foxg1",
                                        "Prox1", "Dlx5", "Sox2", "Alpl", "Fut4"), 
        pt.size = 0, assay = "RNA", stack = T, log = T, flip = F, fill.by = "ident")


VlnPlot(sample.integrated, features = c("Tmprss3", "Pvalb", "Pxdc1", "Veph1",
                                        "Ormdl3", "Pla2g7", "Fbxo2", "Pcdh15"), 
        pt.size = 0, assay = "RNA", adjust = 1, stack = T, flip = T, fill.by = "ident", split.by = "sample")


VlnPlot(cluster.averages, features = c("Lamp1", "Cst3", "Mical1", "Fabp3"), 
        pt.size = 0, assay = "RNA")

#Paper figures
FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Galm", "Col3a1", "Vmo1"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6,
            ncol = 3
)

FeaturePlot(sample.integrated, 
                  reduction = "umap", 
                  features = c("Ceacam16", "Epyc", "Plp1"),
                  label = T,
                  cols = c("lightgray", "red"),
                  coord.fixed = TRUE,
                  min.cutoff = 'q10',
                  max.cutoff = 'q90',
                  repel = T,
                  pt.size = 0.6,
                  ncol = 3
)

FeaturePlot(sample.integrated, 
                  reduction = "umap", 
                  features = c("Gdf15", "Cx3cr1", "Pvalb"),
                  label = T,
                  cols = c("lightgray", "red"),
                  coord.fixed = TRUE,
                  min.cutoff = 'q10',
                  max.cutoff = 'q90',
                  repel = T,
                  pt.size = 0.6,
                  ncol = 3
)

VlnPlot(sample.integrated, features = c("Galm", "Col3a1", "Vmo1", "Car14", "Epyc", "Plp1", "Gdf15", "Cx3cr1", "Pvalb"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, flip = T, fill.by = "ident")

VlnPlot(sample.integrated, features = c("Galm", "Col3a1", "Vmo1", "Ceacam16", "Epyc", "Plp1", "Gdf15", "Cx3cr1", "Pvalb"), 
        pt.size = 0, assay = "RNA", adjust = 2, stack = T, log = F, flip = F, fill.by = "ident", y.max = 100)

VlnPlot(sample.integrated, features = c("Galm", "Col3a1", "Vmo1", "Ceacam16", "Epyc", "Plp1", "Gdf15", "Cx3cr1", "Pvalb"), 
        pt.size = 0, assay = "SCT", adjust = 2, stack = T, log = F, flip = T, fill.by = "ident", y.max = 100)

FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = "Myo7a",
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.8,
            split.by = "stim"
)

FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = "Espn",
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            pt.size = 0.6,
            split.by = "stim"
)

VlnPlot(sample.integrated, features = "Myo7a", pt.size = 0, 
        assay = "RNA", adjust = 1, split.by = "sample", log = F, flip = F, fill.by = "ident", y.max = 10)

VlnPlot(sample.integrated, features = "Espn", pt.size = 0, 
        assay = "RNA", adjust = 1, split.by = "sample", log = F, flip = F, fill.by = "ident", y.max = 20)

FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Galm", "Nfix", "Epcam", "Cp", "Col3a1", "Fstl1", "Vmo1", "Meis2", "Fut9"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            ncol = 3,
            pt.size = 0.4
)


FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Bace2", "Ceacam16", "Car14", "Epyc", "Slc39a8", "Lum", "Plp1", "Pmp22", "Sostdc1"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            ncol = 3,
            pt.size = 0.4
)


FeaturePlot(sample.integrated, 
            reduction = "umap", 
            features = c("Gdf15", "Ciart", "Hspa1a", "Cx3cr1", "Emilin2", "Fcgr1", "Pvalb", "Atp8a2", "Nefl"),
            label = T,
            cols = c("lightgray", "red"),
            coord.fixed = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            repel = T,
            ncol = 3,
            pt.size = 0.4
)

VlnPlot(sample.integrated, features = c("Nefl", "Atp8a2", "Pvalb",
                                        "Fcgr1", "Emilin2","Cx3cr1", 
                                        "Hspa1a", "Ciart", "Gdf15",
                                        "Sostdc1", "Pmp22" , "Plp1",
                                        "Lum", "Slc39a8", "Epyc",
                                        "Car14", "Ceacam16", "Bace2",
                                        "Fut9", "Meis2", "Vmo1",
                                        "Fstl1", "Col3a1" , "Cp",
                                        "Epcam", "Nfix","Galm"), 
        pt.size = 0, assay = "RNA", log = T, adjust = 1, stack = T, flip = F, fill.by = "ident")
