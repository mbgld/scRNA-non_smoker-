# set-up
library(pheatmap); library(Seurat); library(dplyr)
YS_ITG_MK.outdir = "I:/Merged/2_Assign/Results/"
load(paste0(YS_ITG_MK.outdir, "ITG_MK_heatmap.RData"))

# Take out marker genes
YS_ITG_All.mk <- FindAllMarkers(YS_Itg.sfj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
YS_ITG_Top3.mk <- YS_ITG_All.mk %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

# get average matrix
YS_Itg1.avg <- AverageExpression(YS_Itg.sfj)

# generate scaled RNA expression matrix and get primitive graph 
YS_Itg2.avg <- as.data.frame(YS_Itg1.avg[["RNA"]][YS_ITG_Top3.mk$gene, ]); mode(YS_Itg2.avg)

# To exaggerate 
YS_Itg3.avg  <- t(scale(t(YS_Itg2.avg)))

########################################################
# End of number matrix
########################################################

# 1. Legends for cluster
column.lgd <- as.data.frame(colnames(YS_Itg3.avg), row.names = colnames(YS_Itg3.avg))
colnames(column.lgd) <- 'Cell type'

# 2. Row legends
df_row.lgd <- as.data.frame(row.names(YS_Itg3.avg), row.names = row.names(YS_Itg3.avg))
colnames(df_row.lgd) <- 'Marker gene'

# 2.1 Marker genes
FB.gen <- c('DCN','MGP','LUM')
EC.gen <- c('SPARCL1','VWF','EPAS1')
MA.gen <- c('TPSB2','TPSAB1','CPA3')
MY.gen <- c('LYZ','APOC1','C1QA')
BC.gen <- c('IGKC','IGHA1','IGLC2')
TC.gen <- c('CCL5','NKG7','GNLY')
EP.gen <- c('SCGB1A1','CAPS','SLPI')
NE.gen <- c('SFTPB','NAPSA','SPINK1')

# 2.2 MARKER GENE GROUPING FOR ANNOTATION
df_row.lgd$'Marker genes' <- 
  ifelse(df_row.lgd$`Marker gene` %in% FB.gen,"FB",
    ifelse(df_row.lgd$`Marker gene` %in% EC.gen,"EC",
      ifelse(df_row.lgd$`Marker gene` %in% MA.gen,"MA",
        ifelse(df_row.lgd$`Marker gene` %in% MY.gen,"MY",
          ifelse(df_row.lgd$`Marker gene` %in% BC.gen,"BC",
            ifelse(df_row.lgd$`Marker gene` %in% TC.gen,"TC",
              ifelse(df_row.lgd$`Marker gene` %in% EP.gen,"EP",
                "NE")))))))

# Extraction for annotation because options require data.frame format! Stupid Work!
row.lgd <- df_row.lgd %>% dplyr::select('Marker genes')

# 3. Set color
my_color <- list(
  'Marker genes' = c('FB'='#F87C74', 'EC'='#D79D19', 'MA'='#92A900',
                     'MY'='#0EBE43', 'BC'='#00C19F', 'TC'='#06BBE3',
                     'EP'='#90A7FF', 'NE'='#FF7ACC'),
  'Cell type' =    c('FB'='#F87C74', 'EC'='#D79D19', 'MA'='#92A900',
                     'MY'='#0EBE43', 'BC'='#00C19F', 'TC'='#06BBE3',
                     'EP'='#90A7FF', 'NE'='#FF7ACC'))

# End products
# Figure 1

jpeg(filename = paste0("C:/Users/User/Desktop/heatmap_Tcell_1.jpeg"), width = 2000, height =4000, res = 600)

pheatmap(Figure1.averages, 
               annotation_col = my_col,
               # annotation_row = my_row,
               annotation_colors = my_color,
               cutree_cols = F, cluster_rows = F, cluster_cols = F,
               # gaps_col = c(1,2,3,4,5,6,7),
               # gaps_row = c(3,6,9,12,15,18,21),
               labels_row = as.expression(newnames),
               annotation_legend = F
               )

dev.off()


# Figure 1-A
jpeg(filename = paste0(YS_ITG_MK.outdir, "heatmap_Major_cell_Marker.jpeg"), width = 3500, height =4000, res = 600)
pheatmap(YS_Itg3.avg, 
         annotation_col = column.lgd,
         annotation_row = row.lgd,
         annotation_colors = my_color,
         cutree_cols = F, cluster_rows = F, cluster_cols = F,
         # gaps_col = c(1,2,3,4,5,6,7),
         gaps_row = c(3,6,9,12,15,18,21),
         # labels_row = as.expression(row.lgd),
         annotation_legend = T
)
dev.off()

save.image(paste0(YS_ITG_MK.outdir, "ITG_MK_heatmap.RData"))
rm(list=ls()); gc(); q()