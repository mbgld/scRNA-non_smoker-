# set up
library(Seurat); library(dplyr)
YS_ITG_MK.outdir = "I:/Merged/2_Assign/Results/"

# Get previous dataset
# load("/media/yschang/W/Merged/1_Combi_Integ/Except_KBY_mito20/Results/1_QC_Combi_ITG_Except_KBY_mito20.RData") # leave YS_Itg.suj only.
load(paste0(YS_ITG_MK.outdir, "YS_ITG_MK.RData"))

YS_Itg.list <- SplitObject(YS_Itg.suj, split.by = "orig.ident")

# Define useful fxs
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_logFC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_logFC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}

# Hiracheal clustering and ClusterTree
DimPlot(YS_Itg.suj, reduction = "umap", label = T)
PlotClusterTree(BuildClusterTree(YS_Itg.suj))

# Keep previous labels
# cluster.id
YS_Itg.suj$subcluster.id <- Idents(YS_Itg.suj)

# Tissue.id
## Tu cell barcodes
YS_Tu.bc <- c(colnames(YS_Itg.list$KJA_Tu),colnames(YS_Itg.list$LYJ_Tu), colnames(YS_Itg.list$SSS_Tu), colnames(YS_Itg.list$NKH_Tu), colnames(YS_Itg.list$YJO_Tu), colnames(YS_Itg.list$LDS_Tu))
## NL cell barcodes
YS_NL.bc <- c(colnames(YS_Itg.list$KJA_NL),colnames(YS_Itg.list$LYJ_NL), colnames(YS_Itg.list$NKH_NL), colnames(YS_Itg.list$YJO_NL), colnames(YS_Itg.list$LDS_NL))
## Assign
Idents(YS_Itg.suj, cells = YS_NL.bc) <- "NL"
Idents(YS_Itg.suj, cells = YS_Tu.bc) <- "Tu"
YS_Itg.suj$tissue.id <- Idents(YS_Itg.suj)

# Reset subcluster.id
Idents(YS_Itg.suj) <- YS_Itg.suj$subcluster.id

#
DefaultAssay(YS_Itg.suj) <- 'RNA'

# Get Conserved Markers 
for (i in levels(YS_Itg.suj)) {
  YS_FCM_i <- FindConservedMarkers(YS_Itg.suj, ident.1 = i, grouping.var = 'orig.ident', assay = 'RNA')
  write.table(x = YS_FCM_i, file = paste0(YS_ITG_MK.outdir, "YS_ITG_", i, "_cnv.mk"), quote = F, row.names = TRUE, sep = "\t")
}

# 5.2 FindAllmarkers
YS_ITG_All.mk <- FindAllMarkers(YS_Itg.suj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
write.table(x = YS_ITG_All.mk, file = paste0(YS_ITG_MK.outdir, "/Integ_YS_ITG_All.mk"), quote = F, row.names = FALSE, sep = "\t")

YS_ITG_All.pw <- FindAllMarkers(YS_Itg.suj, min.pct = 0.25, only.pos = TRUE , logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
YS_ITG_All <- merge(YS_ITG_All.mk, YS_ITG_All.pw, by="row.names", all.x = F, all.y = F)
YS_ITG_All <- subset(YS_ITG_All, select = c('Row.names', 'cluster.x', 'avg_logFC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
YS_ITG_All <- YS_ITG_All[rev(order(YS_ITG_All$cluster.x, YS_ITG_All$avg_logFC.x, YS_ITG_All$power)), ]
write.table(x = YS_ITG_All, file = paste0(YS_ITG_MK.outdir, "/Integ_YS_ITG_All.mp"), quote = F, sep = "\t")

# 5.3 Find each cluster markers. 
for (idx in levels(YS_Itg.suj)) {
  a <- FindMarkers(YS_Itg.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE)
  b <- FindMarkers(YS_Itg.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE, test.use = "roc")
  c <- MK_merge(a, b)
  assign(paste0("YS_ITG_",idx,".mp"), c)
  write.table(x = c, file = paste0(YS_ITG_MK.outdir, "YS_ITG_", idx, ".mp"), quote = F, row.names = FALSE, sep = "\t")
  rm(a, b, c)
}

################################################################################
################################################################################
################################################################################
# assign celltype id
YS_Itg.sfj <- YS_Itg.suj

# Cancer (CA) and Premalignant (PM) barcode
# KJA
KJA_CA.bc <- readLines('/media/yschang/T/KJA/Tumor_Results/2_InferCNV/subcluster/KJA_CA.bc')
KJA_PM.bc <- readLines('/media/yschang/T/KJA/Tumor_Results/2_InferCNV/subcluster/KJA_PM.bc')

# LYJ
LYJ_CA.bc <- readLines('/media/yschang/T/LYJ/Combi_Results/2_InferCNV/subcluster/LYJ_CA.bc')
LYJ_PM.bc <- readLines('/media/yschang/T/LYJ/Combi_Results/2_InferCNV/subcluster/LYJ_PM.bc')

# YJO
YJO_CA.bc <- readLines('/media/yschang/S/YJO/Combi_Results/2_InferCNV/Results/subcluster/YJO_CA.bc')
YJO_PM.bc <- readLines('/media/yschang/S/YJO/Combi_Results/2_InferCNV/Results/subcluster/YJO_PM.bc')

# SSS
SSS_CA1.bc <- readLines('/media/yschang/S/SSS/Tumor_Results/2_InferCNV/subcluster/SSS_CA1.bc')
SSS_CA2.bc <- readLines('/media/yschang/S/SSS/Tumor_Results/2_InferCNV/subcluster/SSS_CA2.bc')
SSS_CA.bc <- c(SSS_CA1.bc, SSS_CA2.bc)

# NKH
NKH_CA.bc <- readLines('/media/yschang/S/NKH/Combi_Results/2_InferCNV/subcluster/NKH_CA.bc')

# LDS
LDS_CA.bc <- readLines('/media/yschang/S/LDS/Combi_Results/2_InferCNV/Results/subcluster/LDS_CA.bc')
LDS_PM.bc <- readLines('/media/yschang/S/LDS/Combi_Results/2_InferCNV/Results/subcluster/LDS_PM.bc')

# Cancer cell barcodes
YS_CA.bc<- c(KJA_CA.bc, LYJ_CA.bc, YJO_CA.bc, SSS_CA.bc, NKH_CA.bc, LDS_CA.bc)
YS_PM.bc<- c(KJA_PM.bc, LYJ_PM.bc, YJO_PM.bc, LDS_PM.bc)
YS_NEO.bc <- c(YS_CA.bc, YS_PM.bc)

# CA(cancer) and PM(premalignant)
YS_Itg.sfj <- SetIdent(object = YS_Itg.sfj, cells = YS_CA.bc, value = "CA")
YS_Itg.sfj <- SetIdent(object = YS_Itg.sfj, cells = YS_PM.bc, value = "PM")

# NE(Neoplastic) = CA + PM
YS_Itg.sfj <- SetIdent(object = YS_Itg.sfj, cells = c(YS_CA.bc, YS_PM.bc), value = "NE")

# EP
YS_EP_temp.bc <- WhichCells(subset(x = YS_Itg.sfj, idents = c('2', '3', '6', '10', '12', '13', '17')))
YS_EP.bc <- setdiff(EP_temp.bc, YS_NEO.bc)
YS_Itg.sfj <- SetIdent(object = YS_Itg.sfj, cells = YS_EP.bc, value = "EP")

# TC
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, '0' = 'TC', '1' = 'TC')

DimPlot(object = YS_Itg.sfj, reduction = "umap", label = T, label.size = 4)

# BC
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, '14' = 'BC')

# MY
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, '4' = 'MY', '5' = 'MY', '7' = 'MY', '15' = 'MY', '16' ='MY')

# MA
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, '8' = 'MA')

# EC
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, '9' = 'EC')

# FB
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, '11' = 'FB')

###########
# Visualize
###########
DimPlot(object = YS_Itg.sfj, cells.highlight= YS_CA.bc, cols.highlight = "#FF0000", reduction = "umap", label = T, label.size = 4)
DimPlot(YS_Itg.sfj, reduction = "umap", label = T, pt.size = 1.0)

YS_Itg.sfj$celltype.id <- Idents(YS_Itg.sfj)

## Assign cell name
# NE (CA + PM)
jpeg(filename = paste0(YS_ITG_MK.outdir, "NE.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c('NE')), cols.highlight = "#FF66CC", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
Plt <- DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c('NE')), cols.highlight = "#FF66CC", pt.size = 0.1, label = F)
postscript(paste0(YS_ITG_MK.outdir, "NE.eps"))
Plt
dev.off()

# EP
jpeg(filename = paste0(YS_ITG_MK.outdir, "EP.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EP")), cols.highlight = "#6699FF",  pt.size = 0.1, label = F)
dev.off()
# Save as EPS
Plt <- DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EP")), cols.highlight = "#6699FF",  pt.size = 0.1, label = F)
postscript(paste0(YS_ITG_MK.outdir, "EP.eps"))
Plt
dev.off()

# FB
jpeg(filename = paste0(YS_ITG_MK.outdir, "FB.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("FB")), cols.highlight = "#FF6600", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
Plt <- DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("FB")), cols.highlight = "#FF6600", pt.size = 0.1, label = F)
postscript(paste0(YS_ITG_MK.outdir, "FB.eps"))
Plt
dev.off()

# EC
jpeg(filename = paste0(YS_ITG_MK.outdir, "EC.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EC")), cols.highlight = "#CC9933", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
Plt <-DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EC")), cols.highlight = "#CC9933", pt.size = 0.1, label = F)
postscript(paste0(YS_ITG_MK.outdir, "EC.eps"))
Plt
dev.off()

# TC
jpeg(filename = paste0(YS_ITG_MK.outdir, "TC.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("TC")), cols.highlight = "#00CCFF", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
Plt <- DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("TC")), cols.highlight = "#00CCFF", pt.size = 0.1, label = F)
postscript(paste0(YS_ITG_MK.outdir, "TC.eps"))
Plt
dev.off()

# BC
jpeg(filename = paste0(YS_ITG_MK.outdir, "BC.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("BC")), cols.highlight = "#66CC99", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
Plt <- DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("BC")), cols.highlight = "#66CC99", pt.size = 0.1, label = F)
postscript(paste0(YS_ITG_MK.outdir, "BC.eps"))
Plt
dev.off()

# MY
jpeg(filename = paste0(YS_ITG_MK.outdir, "MY.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MY")), cols.highlight = "#009900", pt.size = 0.1,label = F)
dev.off()
# Save as EPS
Plt <- DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MY")), cols.highlight = "#009900", pt.size = 0.1,label = F)
postscript(paste0(YS_ITG_MK.outdir, "MY.eps"))
Plt
dev.off()

# MA
jpeg(filename = paste0(YS_ITG_MK.outdir, "MA.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MA")), cols.highlight = "#999900", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
Plt<-DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MA")), cols.highlight = "#999900", pt.size = 0.1, label = F)
postscript(paste0(YS_ITG_MK.outdir, "MA.eps"))
Plt
dev.off()

# visualize
levels(YS_Itg.sfj) <- c("NE", "FB", "EC", "MA", "MY", "BC", "TC", "EP")
ColorPalate = c('#FF66CC', '#FF6600', '#CC9933', '#999900', '#009900', '#66CC99', '#00CCFF', '#6699FF')

# Major DimPlot
jpeg(filename = paste0(YS_ITG_MK.outdir, "Major_DimPlot.jpeg"), width=5000, height=3500, res = 1200)
DimPlot(YS_Itg.sfj, reduction = "umap", cols = ColorPalate, label = T, label.size = 0.5)
dev.off()
# Save as EPS
Plt<-DimPlot(YS_Itg.sfj, reduction = "umap", cols = ColorPalate, label = T, label.size = 0.5)
postscript(paste0(YS_ITG_MK.outdir, "Major_DimPlot.eps"))
Plt
dev.off()

# Split by tissue.id
jpeg(filename = paste0(YS_ITG_MK.outdir, "DimPlot_by_tissue.jpeg"), width=10000, height=6000, res = 1200)
DimPlot(YS_Itg.sfj, reduction = "umap", cols = ColorPalate, label = F, label.size = 0.5, split.by = 'tissue.id', ncol = 1)
dev.off()
# Save as EPS
Plt<-DimPlot(YS_Itg.sfj, reduction = "umap", cols = ColorPalate, label = F, label.size = 0.5, split.by = 'tissue.id', ncol = 2)
postscript(paste0(YS_ITG_MK.outdir, "DimPlot_by_tissue.eps"))
Plt
dev.off()

# DimPlot of each tissue
unique(Idents(YS_Itg.sfj))
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value='tissue.id')


Tu_Itg.sfj <- subset(YS_Itg.sfj, idents = 'Tu')

# NL dimplot
NL_Itg.sfj <- subset(YS_Itg.sfj, idents = 'NL')
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value='tissue.id')
jpeg(filename = paste0(YS_ITG_MK.outdir, "NL_DimPlot.jpeg"), width=10000, height=6000, res = 1200)
DimPlot(NL_Itg.sfj, reduction = "umap", label = T, label.size = 0.5, ncol = 1)
dev.off()


# dotplot
markers.to.plot <- YS_ITG_Top3.mk$gene
DotPlot(YS_Itg.sfj, features = markers.to.plot, cols = c('RdYlBu'), dot.scale = 8) + RotatedAxis()

# Save
save.image(paste0(YS_ITG_MK.outdir, "YS_ITG_MK.RData"))
saveRDS(YS_Itg.averages, paste0(YS_ITG_MK.outdir,"YS_Itg.averages.rds"))

rm(list=ls());gc()
q()
