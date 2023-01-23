options(stringsAsFactors = FALSE)
options(scipen = 100)
.libPaths("~/anaconda3/envs/R4.1/lib/R/library")

library(Matrix,lib.loc = "~/anaconda3/envs/R4.1/lib/R/library")
library(BuenColors)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(VennDiagram)
library(Rmisc)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
# library(chromVAR)
# library(terra, lib.loc = "/usr/local/lib/R/site-library")
# library(cicero)
# library(ArchR)
# library(SnapATAC)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm10)


DMmouse_CRE <- readRDS("~/zjw/20220307DMmouse/Seurat_merge_CRE/results/combined_cca_40.RDS")
genes_remove <- rownames(DMmouse_CRE@assays$RNA@counts)[rowSums(DMmouse_CRE@assays$RNA@counts)==0]
DMmouse_CRE <- subset(DMmouse_CRE, features = rownames(DMmouse_CRE@assays$RNA@counts)[!rownames(DMmouse_CRE@assays$RNA@counts) %in% genes_remove])

CRE <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/cm.aggregate/annotate/DM/bed/DM.CRE.annot.bed.gz",sep="\t")
CRE <- CRE[!CRE$V4 %in% genes_remove,]



DMmouse_CRE.snap <- readRDS("~/zjw/20220307DMmouse/DMmouse_CRE.snap.RDS")




mmat <- data.frame(DMmouse_CRE.snap@mmat,row.names = DMmouse_CRE.snap@barcode)
mmat <- mmat[rownames(DMmouse_CRE@meta.data),]

DMmouse_CRE@meta.data <- cbind(DMmouse_CRE@meta.data,mmat)

temp <- CreateAssayObject(counts = t(mmat))

DMmouse_CRE[["mmat"]] <- temp

DMmouse_CRE.CTRL <- subset(DMmouse_CRE, subset = var == "CTRL")


a <- DMmouse_CRE.CTRL@meta.data

t.test_motif_CTRL <- data.frame()
for (i in 10:395) {
  for (celltype in unique(a$celltype)) {
    
    t.test_motif_CTRL <- rbind(t.test_motif_CTRL,
                               data.frame(p.val=t.test(a[a$celltype==celltype,i],a[!a$celltype==celltype,i])$p.value,
                                          ave_delta=mean(a[a$celltype==celltype,i])-mean(a[!a$celltype==celltype,i]),
                                          mean_1=mean(a[a$celltype==celltype,i]),
                                          mean_2=mean(a[!a$celltype==celltype,i]),
                                          cluster=celltype,
                                          gene=colnames(a)[i]
                               )
    )
  }
}

t.test_motif_CTRL$adj.P.Val <- p.adjust(t.test_motif_CTRL$p.val,method = "BH") 

t.test_motif_CTRL$cluster <- factor(t.test_motif_CTRL$cluster,levels = c("Rod","Cone","Rod BC","Cone BC","MG","AC","HC","RGC","Micro","VEC","Pericyte"))

t.test_motif_CTRL_DE <- t.test_motif_CTRL[t.test_motif_CTRL$p.val<0.05 & t.test_motif_CTRL$ave_delta>0.02,]

# t.test_motif_CTRL_DE_sub <- t.test_motif_CTRL_DE[t.test_motif_CTRL_DE$gene %in% names(table(t.test_motif_CTRL_DE$gene)[table(t.test_motif_CTRL_DE$gene)==1]),]

motif_celltype_df <- data.frame(t.test_motif_CTRL$mean_1[t.test_motif_CTRL$cluster==levels(t.test_motif_CTRL$cluster)[1]])
colnames(motif_celltype_df) <- levels(t.test_motif_CTRL$cluster)[1]

for (i in levels(t.test_motif_CTRL$cluster)[-1]) {
  temp <- data.frame(t.test_motif_CTRL$mean_1[t.test_motif_CTRL$cluster==i])
  motif_celltype_df <- cbind(motif_celltype_df,temp)
  
}
colnames(motif_celltype_df)[-1] <- levels(t.test_motif_CTRL$cluster)[-1]

rownames(motif_celltype_df) <- unique(t.test_motif_CTRL$gene)

motif_celltype_df <- motif_celltype_df[,levels(Idents(DMmouse_CRE))]

gene_list <- c()
for (i in colnames(motif_celltype_df)) {
  temp <- t.test_motif_CTRL[t.test_motif_CTRL$cluster==i,]
  temp <- temp[order(temp$ave_delta,decreasing = T),]
  gene_list <- c(gene_list,temp$gene)
}


df <- data.frame()
for (i in gene_list) {
  df <- rbind(df,
              motif_celltype_df[rownames(motif_celltype_df)==i,])
}

colnames(df) <- names(motif_celltype_df[1,])
# rownames(df) <- gene_list


# df <- ave_exp[rownames(ave_exp) %in% c("KRT14","KRT15","KRT3","KRT12","ZNF281"),colnames(ave_exp) %in% c("LSC","LEC","ETAC","CTAC","CEC")]



library(matrixStats)
df <- (df-rowMeans(df))/(rowSds(as.matrix(df)))[row(df)]



gap <- table(t.test_motif_CTRL_DE$cluster)
# gap <- gap[ord]
for (i in 2:length(gap)) {
  gap[i] <- gap[i]+gap[i-1]
}

df_new <- data.frame()

for (i in 1:length(gap)) {
  if (i==1) {
    temp <- df[1:gap[i],]
  } else {
    temp <- df[(gap[i-1]+1):gap[i],]
  }
  
  temp <- temp[order(temp[,i],decreasing = T),]
  
  df_new <- rbind(df_new,temp)
}


pheatmap(df_new,
         scale = "none",
         cluster_cols = F,
         cluster_rows = F,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         show_rownames = F,
         show_colnames = F,
         treeheight_row = F,
         treeheight_col = F,
         legend = F,
         border = F,
         gaps_row = gap)


# p <- DoHeatmap(DMmouse_CRE.CTRL,
#                features = rownames(DMmouse_CRE.CTRL@assays$mmat@counts),
#                cells = NULL,
#                # group.by = "ident",
#                group.bar = TRUE,
#                group.colors = NULL,
#                disp.min = -2.5,
#                disp.max = NULL,
#                slot = "data",
#                assay = "mmat",
#                label = TRUE,
#                size = 5.5,
#                hjust = 0,
#                angle = 45,
#                raster = TRUE,
#                draw.lines = TRUE,
#                lines.width = NULL,
#                group.bar.height = 0.02,
#                combine = TRUE
# )
# 
# 
# ggplot2::ggsave(p, filename = "temp.png", height = 30, width = 20, limitsize = FALSE)


temp <- t.test_motif_CTRL_DE[t.test_motif_CTRL_DE$cluster=="Rod",]

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

DoHeatmap()


FeaturePlot(DMmouse_CRE,features = c("MA0090.2_TEAD1"), min.cutoff = "q1", max.cutoff = "q99",order = F)+
  scale_color_gradientn(colours = SpatialColors(n = 100)[51:100])
  

runChromVAR