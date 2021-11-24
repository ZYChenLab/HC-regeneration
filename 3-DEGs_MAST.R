# DEGs between samples
sample.integrated <- readRDS("~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Seurat_Combined.rds")

new.cluster.ids <- c("IdC", "CCOS","RMC","DC","CCOS","KO","CCOS","HeC","SC","RMC","TBRM","TBRM","DC","HC")
names(new.cluster.ids) <- levels(sample.integrated)
sample.integrated <- RenameIdents(sample.integrated, `0` = "IdC", `1` = "CCOS", `2` = "RMC", 
                                  `3` = "DC", `4` = "CCOS", `5` = "KO", `6` = "CCOS", `7` = "HeC", 
                                  `8` = "SC", `9` = "RMC", `10` = "TBRM", `11` = "TBRM", 
                                  `12` = "DC", `13` = "HC")

DefaultAssay(sample.integrated) <- "RNA"
sample.integrated <- NormalizeData(sample.integrated, verbose = FALSE)
Idents(sample.integrated) <- "sample"
sample.integrated@active.assay <- "RNA"
DEG_D1C_MAST <- FindMarkers(sample.integrated, ident.1 = "Dox1", ident.2 = "Control", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_D1C_MAST <- DEG_D1C_MAST[which(DEG_D1C_MAST$p_val_adj < 0.05),]
DEG_D1C_MAST <- DEG_D1C_MAST[is.finite(rowSums(DEG_D1C_MAST)),]

write.csv(DEG_D1C_MAST, file = "~/Desktop/DOX1up_MAST.csv")

DEG_CD1_MAST <- FindMarkers(sample.integrated, ident.1 = "Control", ident.2 = "Dox1", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CD1_MAST <- DEG_CD1_MAST[which(DEG_CD1_MAST$p_val_adj < 0.05),]
DEG_CD1_MAST <- DEG_CD1_MAST[is.finite(rowSums(DEG_CD1_MAST)),]

write.csv(DEG_CD1_MAST, file = "~/Desktop/MAST/CTRL1up_MAST.csv")

DEG_D2C_MAST <- FindMarkers(sample.integrated, ident.1 = "Dox2", ident.2 = "Control", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_D2C_MAST <- DEG_D2C_MAST[which(DEG_D2C_MAST$p_val_adj < 0.05),]
DEG_D2C_MAST <- DEG_D2C_MAST[is.finite(rowSums(DEG_D2C_MAST)),]

write.csv(DEG_D2C_MAST, file = "~/Desktop/MAST/DOX2up_MAST.csv")

DEG_CD2_MAST <- FindMarkers(sample.integrated, ident.1 = "Control", ident.2 = "Dox2", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CD2_MAST <- DEG_CD2_MAST[which(DEG_CD2_MAST$p_val_adj < 0.05),]
DEG_CD2_MAST <- DEG_CD2_MAST[is.finite(rowSums(DEG_CD2_MAST)),]

write.csv(DEG_CD2_MAST, file = "~/Desktop/MAST/CTRL2up_MAST.csv")

dev.off()

#DEGs between individual clusters 
sample.integrated <- readRDS("~/Desktop/WTvsDox-CellRanger-Outs/Data_combined_all/SCT/SCT-III/Seurat_Combined.rds")

new.cluster.ids <- c("IdC", "CCOS","RMC","DC","CCOS","KO","CCOS","HeC","SC","RMC","TBRM","TBRM","DC","HC")
names(new.cluster.ids) <- levels(sample.integrated)
sample.integrated <- RenameIdents(sample.integrated, `0` = "IdC", `1` = "CCOS", `2` = "RMC", 
                                  `3` = "DC", `4` = "CCOS", `5` = "KO", `6` = "CCOS", `7` = "HeC", 
                                  `8` = "SC", `9` = "RMC", `10` = "TBRM", `11` = "TBRM", 
                                  `12` = "DC", `13` = "HC")

DefaultAssay(sample.integrated) <- "RNA"
sample.integrated <- NormalizeData(sample.integrated, verbose = FALSE)

sample.integrated$celltype.stim <- paste(Idents(sample.integrated), sample.integrated$stim, sep = "_")
sample.integrated$celltype <- Idents(sample.integrated)
Idents(sample.integrated) <- "celltype.stim"
levels(sample.integrated)

DEG_DoxIdC <- FindMarkers(sample.integrated, ident.1 = "IdC_Dox", ident.2 = "IdC_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxIdC <- DEG_DoxIdC[which(DEG_DoxIdC$p_val_adj < 0.05),]
DEG_DoxIdC <- DEG_DoxIdC[is.finite(rowSums(DEG_DoxIdC)),]
write.csv(DEG_DoxIdC, file = "~/Desktop/MAST/DEG_DoxIdC.csv")

DEG_CtrlIdC <- FindMarkers(sample.integrated, ident.1 = "IdC_Control", ident.2 = "IdC_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlIdC <- DEG_CtrlIdC[which(DEG_CtrlIdC$p_val_adj < 0.05),]
DEG_CtrlIdC <- DEG_CtrlIdC[is.finite(rowSums(DEG_CtrlIdC)),]
write.csv(DEG_CtrlIdC, file = "~/Desktop/MAST/DEG_CtrlIdC.csv")

DEG_DoxCCOS <- FindMarkers(sample.integrated, ident.1 = "CCOS_Dox", ident.2 = "CCOS_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxCCOS <- DEG_DoxCCOS[which(DEG_DoxCCOS$p_val_adj < 0.05),]
DEG_DoxCCOS <- DEG_DoxCCOS[is.finite(rowSums(DEG_DoxCCOS)),]
write.csv(DEG_DoxCCOS, file = "~/Desktop/MAST/DEG_DoxCCOS.csv")

DEG_CtrlCCOS <- FindMarkers(sample.integrated, ident.1 = "CCOS_Control", ident.2 = "CCOS_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlCCOS <- DEG_CtrlCCOS[which(DEG_CtrlCCOS$p_val_adj < 0.05),]
DEG_CtrlCCOS <- DEG_CtrlCCOS[is.finite(rowSums(DEG_CtrlCCOS)),]
write.csv(DEG_CtrlCCOS, file = "~/Desktop/MAST/DEG_CtrlCCOS.csv")

DEG_DoxHeC <- FindMarkers(sample.integrated, ident.1 = "HeC_Dox", ident.2 = "HeC_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxHeC <- DEG_DoxHeC[which(DEG_DoxHeC$p_val_adj < 0.05),]
DEG_DoxHeC <- DEG_DoxHeC[is.finite(rowSums(DEG_DoxHeC)),]
write.csv(DEG_DoxHeC, file = "~/Desktop/MAST/DEG_DoxHeC.csv")

DEG_CtrlHeC <- FindMarkers(sample.integrated, ident.1 = "HeC_Control", ident.2 = "HeC_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlHeC <- DEG_CtrlHeC[which(DEG_CtrlHeC$p_val_adj < 0.05),]
DEG_CtrlHeC <- DEG_CtrlHeC[is.finite(rowSums(DEG_CtrlHeC)),]
write.csv(DEG_CtrlHeC, file = "~/Desktop/MAST/DEG_CtrlHeC.csv")

DEG_DoxDC <- FindMarkers(sample.integrated, ident.1 = "DC_Dox", ident.2 = "DC_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxDC <- DEG_DoxDC[which(DEG_DoxDC$p_val_adj < 0.05),]
DEG_DoxDC <- DEG_DoxDC[is.finite(rowSums(DEG_DoxDC)),]
write.csv(DEG_DoxDC, file = "~/Desktop/MAST/DEG_DoxDC.csv")

DEG_CtrlDC <- FindMarkers(sample.integrated, ident.1 = "DC_Control", ident.2 = "DC_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlDC <- DEG_CtrlDC[which(DEG_CtrlDC$p_val_adj < 0.05),]
DEG_CtrlDC <- DEG_CtrlDC[is.finite(rowSums(DEG_CtrlDC)),]
write.csv(DEG_CtrlDC, file = "~/Desktop/MAST/DEG_CtrlDC.csv")

DEG_DoxSC <- FindMarkers(sample.integrated, ident.1 = "SC_Dox", ident.2 = "SC_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxSC <- DEG_DoxSC[which(DEG_DoxSC$p_val_adj < 0.05),]
DEG_DoxSC <- DEG_DoxSC[is.finite(rowSums(DEG_DoxSC)),]
write.csv(DEG_DoxSC, file = "~/Desktop/MAST/DEG_DoxSC.csv")

DEG_CtrlSC <- FindMarkers(sample.integrated, ident.1 = "SC_Control", ident.2 = "SC_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlSC <- DEG_CtrlSC[which(DEG_CtrlSC$p_val_adj < 0.05),]
DEG_CtrlSC <- DEG_CtrlSC[is.finite(rowSums(DEG_CtrlSC)),]
write.csv(DEG_CtrlSC, file = "~/Desktop/MAST/DEG_CtrlSC.csv")

DEG_DoxTBRM <- FindMarkers(sample.integrated, ident.1 = "TBRM_Dox", ident.2 = "TBRM_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxTBRM <- DEG_DoxTBRM[which(DEG_DoxTBRM$p_val_adj < 0.05),]
DEG_DoxTBRM <- DEG_DoxTBRM[is.finite(rowSums(DEG_DoxTBRM)),]
write.csv(DEG_DoxTBRM, file = "~/Desktop/MAST/DEG_DoxTBRM.csv")

DEG_CtrlTBRM <- FindMarkers(sample.integrated, ident.1 = "TBRM_Control", ident.2 = "TBRM_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlTBRM <- DEG_CtrlTBRM[which(DEG_CtrlTBRM$p_val_adj < 0.05),]
DEG_CtrlTBRM <- DEG_CtrlTBRM[is.finite(rowSums(DEG_CtrlTBRM)),]
write.csv(DEG_CtrlTBRM, file = "~/Desktop/MAST/DEG_CtrlTBRM.csv")

DEG_DoxRMC <- FindMarkers(sample.integrated, ident.1 = "RMC_Dox", ident.2 = "RMC_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxRMC <- DEG_DoxRMC[which(DEG_DoxRMC$p_val_adj < 0.05),]
DEG_DoxRMC <- DEG_DoxRMC[is.finite(rowSums(DEG_DoxRMC)),]
write.csv(DEG_DoxRMC, file = "~/Desktop/MAST/DEG_DoxRMC.csv")

DEG_CtrlRMC <- FindMarkers(sample.integrated, ident.1 = "RMC_Control", ident.2 = "RMC_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlRMC <- DEG_CtrlRMC[which(DEG_CtrlRMC$p_val_adj < 0.05),]
DEG_CtrlRMC <- DEG_CtrlRMC[is.finite(rowSums(DEG_CtrlRMC)),]
write.csv(DEG_CtrlRMC, file = "~/Desktop/MAST/DEG_CtrlRMC.csv")

DEG_DoxKO <- FindMarkers(sample.integrated, ident.1 = "KO_Dox", ident.2 = "KO_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxKO <- DEG_DoxKO[which(DEG_DoxKO$p_val_adj < 0.05),]
DEG_DoxKO <- DEG_DoxKO[is.finite(rowSums(DEG_DoxKO)),]
write.csv(DEG_DoxKO, file = "~/Desktop/MAST/DEG_DoxKO.csv")

DEG_CtrlKO <- FindMarkers(sample.integrated, ident.1 = "KO_Control", ident.2 = "KO_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlKO <- DEG_CtrlKO[which(DEG_CtrlKO$p_val_adj < 0.05),]
DEG_CtrlKO <- DEG_CtrlKO[is.finite(rowSums(DEG_CtrlKO)),]
write.csv(DEG_CtrlKO, file = "~/Desktop/MAST/DEG_CtrlKO.csv")

DEG_DoxHC <- FindMarkers(sample.integrated, ident.1 = "HC_Dox", ident.2 = "HC_Control", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_DoxHC <- DEG_DoxHC[which(DEG_DoxHC$p_val_adj < 0.05),]
DEG_DoxHC <- DEG_DoxHC[is.finite(rowSums(DEG_DoxHC)),]
write.csv(DEG_DoxHC, file = "~/Desktop/MAST/DEG_DoxHC.csv")

DEG_CtrlHC <- FindMarkers(sample.integrated, ident.1 = "HC_Control", ident.2 = "HC_Dox", assay = "RNA", logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_CtrlHC <- DEG_CtrlHC[which(DEG_CtrlHC$p_val_adj < 0.05),]
DEG_CtrlHC <- DEG_CtrlHC[is.finite(rowSums(DEG_CtrlHC)),]
write.csv(DEG_CtrlHC, file = "~/Desktop/MAST/DEG_CtrlHC.csv")