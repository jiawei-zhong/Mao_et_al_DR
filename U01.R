.libPaths("/disk2/user/jizhon/anaconda3/envs/R4.1/lib/R/library")


library(Matrix)
library(BuenColors)
library(Seurat)
library(SeuratData)
library(ggplot2)heat
library(patchwork)
library(dplyr)
library(stringr)
library(VennDiagram)
library(Rmisc)
library(ggpubr)
library(pheatmap)



# DMmouse <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge/results/combined_cca_40.RDS")
# DMmouse_old <- readRDS("~/zjw/20220307_DMmouse/combined_cca_40_old_normal.CDS")
# # DimPlot(DMmouse_old,label = T)
# 
# # DMmouse_meta <- DMmouse@meta.data
# 
# # 
# # temp <- Read10X("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/count/SRR12443300/matrix")
# # dim(temp)
# # colnames(temp) <- gsub(pattern = "-",replacement = ".",x = colnames(temp))
# # temp <- temp[,colnames(temp) %in% (rownames(DMmouse_meta)[DMmouse_meta$multi=="DM1"] %>% gsub(pattern = "DM1_",replacement = ""))]
# # dim(temp)
# # write.table(temp,"~/zjw/20220307_DMmouse/Seurat_merge_CRE/DM1.tsv",quote = F,row.names = T,col.names = T,sep = "\t")
# # 
# # 
# # temp <- Read10X("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/count/SRR12443301/matrix")
# # dim(temp)
# # colnames(temp) <- gsub(pattern = "-",replacement = ".",x = colnames(temp))
# # temp <- temp[,colnames(temp) %in% (rownames(DMmouse_meta)[DMmouse_meta$multi=="DM2"] %>% gsub(pattern = "DM2_",replacement = ""))]
# # dim(temp)
# # write.table(temp,"~/zjw/20220307_DMmouse/Seurat_merge_CRE/DM2.tsv",quote = F,row.names = T,col.names = T,sep = "\t")
# # 
# # 
# # temp <- Read10X("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/count/SRR12443302/matrix")
# # dim(temp)
# # colnames(temp) <- gsub(pattern = "-",replacement = ".",x = colnames(temp))
# # temp <- temp[,colnames(temp) %in% (rownames(DMmouse_meta)[DMmouse_meta$multi=="DM3"] %>% gsub(pattern = "DM3_",replacement = ""))]
# # dim(temp)
# # write.table(temp,"~/zjw/20220307_DMmouse/Seurat_merge_CRE/DM3.tsv",quote = F,row.names = T,col.names = T,sep = "\t")
# # 
# # 
# # temp <- Read10X("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/count/SRR12443303/matrix")
# # dim(temp)
# # colnames(temp) <- gsub(pattern = "-",replacement = ".",x = colnames(temp))
# # temp <- temp[,colnames(temp) %in% (rownames(DMmouse_meta)[DMmouse_meta$multi=="CTRL1"] %>% gsub(pattern = "CTRL1_",replacement = ""))]
# # dim(temp)
# # write.table(temp,"~/zjw/20220307_DMmouse/Seurat_merge_CRE/CTRL1.tsv",quote = F,row.names = T,col.names = T,sep = "\t")
# # 
# # 
# # temp <- Read10X("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/count/SRR12443304/matrix")
# # dim(temp)
# # colnames(temp) <- gsub(pattern = "-",replacement = ".",x = colnames(temp))
# # temp <- temp[,colnames(temp) %in% (rownames(DMmouse_meta)[DMmouse_meta$multi=="CTRL2"] %>% gsub(pattern = "CTRL2_",replacement = ""))]
# # dim(temp)
# # write.table(temp,"~/zjw/20220307_DMmouse/Seurat_merge_CRE/CTRL2.tsv",quote = F,row.names = T,col.names = T,sep = "\t")
# # 
# # 
# # temp <- Read10X("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/count/SRR12443305/matrix")
# # dim(temp)
# # colnames(temp) <- gsub(pattern = "-",replacement = ".",x = colnames(temp))
# # temp <- temp[,colnames(temp) %in% (rownames(DMmouse_meta)[DMmouse_meta$multi=="CTRL3"] %>% gsub(pattern = "CTRL3_",replacement = ""))]
# # dim(temp)
# # write.table(temp,"~/zjw/20220307_DMmouse/Seurat_merge_CRE/CTRL3.tsv",quote = F,row.names = T,col.names = T,sep = "\t")
# 
# 
# 
# # DMmouse <-  RunUMAP(object = DMmouse, reduction = "pca", dims = 1:30)
# # DMmouse <-  FindNeighbors(object = DMmouse, reduction = "pca", dims = 1:30)
# # DMmouse <-  FindClusters(DMmouse, resolution = 0.4)
# # DimPlot(DMmouse,label = T)
# 
# # DefaultAssay(DMmouse) <- "RNA"
# 
# # DotPlot(object = DMmouse,features = c("Rho","Gnat1","Scgn","Sebox","Glul","Rlbp1","Slc32a1","Gad1","Pou4f1","Sncg","Opn1mw","Opn1sw","Aif1","Cx3cr1","Cdh5","Pecam1","Kcnj8","Pdgfrb","Lhx1"))
# 
# DMmouse_old <- RenameIdents(object = DMmouse_old,
#                             `0` = "Rod",              
#                             `1` = "Rod",                
#                             `2` = "Rod",               
#                             `3` = "MG",       
#                             `4` = "Rod BC",                  
#                             `5` = "Cone",                   
#                             `6` = "AC",              
#                             `7` = "MG",                       
#                             `8` = "Cone BC",                   
#                             `9` = "Cone BC",                   
#                             `10` = "Cone BC",                   
#                             `11` = "Cone BC",                   
#                             `12` = "AC",                   
#                             `13` = "Cone BC",                   
#                             `14` = "Cone BC",                   
#                             `15` = "Cone BC",                   
#                             `16` = "Micro",                   
#                             `17` = "AC",                   
#                             `18` = "VEC",                   
#                             `19` = "HC",                   
#                             `20` = "Rod BC",                   
#                             `21` = "MG",                   
#                             `22` = "Pericyte",                   
#                             `23` = "MG",                   
#                             `24` = "RGC",                   
#                             `25` = "Rod")
# 
# 
# DMmouse <- RenameIdents(object = DMmouse,
#                             `0` = "Rod",              
#                             `1` = "Rod",                
#                             `2` = "MG",               
#                             `3` = "Rod BC",       
#                             `4` = "Cone",                  
#                             `5` = "AC",                   
#                             `6` = "Cone BC",              
#                             `7` = "Cone BC",                       
#                             `8` = "Cone BC",                   
#                             `9` = "Cone BC",                   
#                             `10` = "AC",                   
#                             `11` = "MG",                   
#                             `12` = "Cone BC",                   
#                             `13` = "Cone BC",                   
#                             `14` = "Cone BC",                   
#                             `15` = "Micro",                   
#                             `16` = "AC",                   
#                             `17` = "Pericyte",                   
#                             `18` = "HC",                   
#                             `19` = "AC",                   
#                             `20` = "MG",                   
#                             `21` = "MG")  
# 
# 
# DimPlot(object = DMmouse_old,label = T)+
#   scale_color_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B","CLP","MEP","HSC")))
# 
# 
# 
# 
# c_temp <- Idents(DMmouse) %>% as.character()
# names(c_temp) <- names(Idents(DMmouse))
# c_temp[names(c_temp) %in% names(Idents(DMmouse_old))[Idents(DMmouse_old)=="RGC"]] <- "RGC"
# c_temp[names(c_temp) %in% names(Idents(DMmouse_old))[Idents(DMmouse_old)=="VEC"]] <- "VEC"
# # c_temp[names(c_temp) %in% names(Idents(DMmouse_old))[Idents(DMmouse_old)=="Pericyte"]] <- "Pericyte"
# c_temp <- factor(c_temp,levels = c("Rod", "Cone", "Rod BC", "Cone BC", "MG", "AC", "HC", "RGC", "Micro", "VEC", "Pericyte"))
# Idents(DMmouse) <- c_temp
# DMmouse$celltype <- Idents(DMmouse)
# 
# DimPlot(object = DMmouse,label = T)+
#   scale_color_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B","CLP","MEP","HSC")))
# 
# 
# 
# # rods [rhodopsin (Rho) and G protein subunit α transducin 1 (Gnat1)], 
# # cone bipolar cells [secretagogin (Scgn)], 
# # rod bipolar cells [Sebox], 
# # Müller glia [glutamate-ammonia ligase (Glul) and Rlbp1], 
# # amacrine cells [solute carrier family 32 member 1 (Slc32a1) and glutamate decarboxylase 1 (Gad1)], 
# # retinal ganglion cells [RGCs] [POU class 4 homeobox 1 (Pou4f1) and synuclein γ (Sncg)], 
# # cones [opsin 1 (cone pigments) medium-wave-sensitive (color blindness, deutan) (Opn1mw) and opsin 1 (cone pigments) short-wave-sensitive (color blindness, tritan) (Opn1sw)], 
# # microglia [allograft inflammatory factor 1 (Aif1) and C-X3-C motif chemokine receptor 1 (Cx3cr1)], 
# # vascular endothelial cells [cadherin 5 (Cdh5) and platelet and endothelial cell adhesion molecule 1 (Pecam1)], 
# # pericytes [potassium voltage-gated channel subfamily J member 8 (Kcnj8) and platelet-derived growth factor receptor β (Pdgfrb)], and 
# # horizontal cells [LIM homeobox 1 (Lhx1)]
# 
# 
# 
# 
# 
DMmouse_CRE <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/combined_cca_40.RDS")

# Idents(DMmouse_CRE) <- Idents(DMmouse)
DMmouse_CRE$celltype <- Idents(DMmouse_CRE)
DimPlot(object = DMmouse_CRE,label = T)+
  scale_color_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B","CLP","MEP","HSC")))


p <- DimPlot(object = DMmouse_CRE) +
  NoLegend() +
  NoAxes() +
  scale_color_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B","CLP","MEP","HSC")))

ggsave(filename="Dimplot_all.png", 
       plot=p, 
       width = 5, height = 5, 
       units = 'in', 
       dpi = 1000,
       device = "png")

p <- DimPlot(object = DMmouse_CRE,split.by = "var",group.by = "var") +
  NoLegend() +
  NoAxes() + 
  scale_color_manual(values=c("blue","red"))


ggsave(filename="Dimplot_split.png", 
       plot=p, 
       width = 10, height = 5.715, 
       units = 'in', 
       dpi = 1000,
       device = "png")


DMmouse_CRE.CTRL <- subset(DMmouse_CRE, subset = var == "CTRL")

mk.CTRL <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/CTRL.marker.RDS")
mk.CTRL <- mk.CTRL[mk.CTRL$avg_log2FC>0.2 & mk.CTRL$p_val_adj<0.05,]

ave_exp <- AverageExpression(DMmouse_CRE.CTRL,assays = "RNA")[[1]]


gene_list <- c()
for (i in unique(mk.CTRL$cluster)) {
  temp <- mk.CTRL[mk.CTRL$cluster==i,]
  temp <- temp[order(temp$avg_log2FC,decreasing = T),]
  gene_list <- c(gene_list,temp[1:50,]$gene)
}


df <- data.frame()
for (i in gene_list) {
  df <- rbind(df,
              ave_exp[rownames(ave_exp)==i,])
}

colnames(df) <- names(ave_exp[1,])
# rownames(df) <- gene_list


# df <- ave_exp[rownames(ave_exp) %in% c("KRT14","KRT15","KRT3","KRT12","ZNF281"),colnames(ave_exp) %in% c("LSC","LEC","ETAC","CTAC","CEC")]



pheatmap(df,
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         show_rownames = F,
         treeheight_row = F,
         treeheight_col = F,
         border = F)

temp <- Idents(DMmouse_CRE.CTRL) %>% table() %>% as.data.frame()

ggplot(data=temp, mapping=aes(x=.,y=Freq)) +
  geom_point(stat= "identity",aes(size=Freq),show.legend = TRUE) +
  scale_fill_continuous(limits=c(0, 20000), breaks=seq(0,20000,by=4000))

ggdotplot(temp, x = ".", y = "Freq",size="Freq")


DMmouse_CRE <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/combined_cca_40.RDS")

# Idents(DMmouse_CRE) <- Idents(DMmouse)
DMmouse_CRE$celltype <- Idents(DMmouse_CRE)

Idents(DMmouse_CRE) <- "multi"

ave_exp <- AverageExpression(DMmouse_CRE,assays = "RNA")[[1]]

ave_exp <- ave_exp[!rowSums(ave_exp)==0,]

cor_mtx <- cor(ave_exp)

# 
# 
# rm(list = "DMmouse")
# rm(list = "DMmouse_old")


# DMmouse_CRE <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/combined_cca_40.RDS")
# 
# DMmouse_CRE@assays$RNA@scale.data <- matrix(nrow = 0,ncol = 0)
# 
# DMmouse_CRE <- DietSeurat(DMmouse_CRE,
#                           counts = TRUE,
#                           data = TRUE,
#                           scale.data = FALSE,
#                           features = NULL,
#                           assays = NULL,
#                           dimreducs = NULL,
#                           graphs = NULL
# )
# 
# DMmouse_CRE.ud <- DMmouse_CRE



# DMmouse_CRE.ud$var[DMmouse_CRE.ud$var=="CTRL"] <- "AGED"
# DMmouse_CRE.ud <- subset(DMmouse_CRE.ud, idents = c("LSC","CTAC", "CEC", "LEC",  "ETAC", "MC"))


# DMmouse_CRE.ud$celltype.stim <- paste(Idents(DMmouse_CRE.ud), DMmouse_CRE.ud$var, sep = "_")
# DMmouse_CRE.ud$celltype <- Idents(DMmouse_CRE.ud)
# Idents(DMmouse_CRE.ud) <- "celltype.stim"
# DefaultAssay(DMmouse_CRE.ud) <- "RNA"

# DMmouse_CRE.ud <- ScaleData(DMmouse_CRE.ud,assay = "RNA")
# DMmouse_CRE.ud <- NormalizeData(DMmouse_CRE.ud,assay = "RNA")


# Rod.diff <-      FindMarkers(DMmouse_CRE.ud, ident.1 = "Rod_DM", ident.2 = "Rod_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0,assay = "RNA")
# Cone.diff <-     FindMarkers(DMmouse_CRE.ud, ident.1 = "Cone_DM", ident.2 = "Cone_DM", verbose = FALSE, logfc.threshold = 0, min.pct = 0,assay = "RNA")
# Rod_BC.diff <-   FindMarkers(DMmouse_CRE.ud, ident.1 = "Rod BC_DM", ident.2 = "Rod BC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0,assay = "RNA")
# Cone_BC.diff <-  FindMarkers(DMmouse_CRE.ud, ident.1 = "Cone BC_DM", ident.2 = "Cone BC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0,assay = "RNA")
# MG.diff <-       FindMarkers(DMmouse_CRE.ud, ident.1 = "MG_DM", ident.2 = "MG_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0,assay = "RNA")
# AC.diff <-       FindMarkers(DMmouse_CRE.ud, ident.1 = "AC_DM", ident.2 = "AC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
# HC.diff <-       FindMarkers(DMmouse_CRE.ud, ident.1 = "HC_DM", ident.2 = "HC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
# RGC.diff <-      FindMarkers(DMmouse_CRE.ud, ident.1 = "RGC_DM", ident.2 = "RGC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
# Micro.diff <-    FindMarkers(DMmouse_CRE.ud, ident.1 = "Micro_DM", ident.2 = "Micro_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
# VEC.diff <-      FindMarkers(DMmouse_CRE.ud, ident.1 = "VEC_DM", ident.2 = "VEC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
# Pericyte.diff <- FindMarkers(DMmouse_CRE.ud, ident.1 = "Pericyte_DM", ident.2 = "Pericyte_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
# 
# 
# saveRDS(Rod.diff,     "Rod.diff.RDS")
# saveRDS(Cond.diff,    "Cond.diff.RDS")
# saveRDS(Rod_BC.diff,  "Rod_BC.diff.RDS")
# saveRDS(Cone_BC.diff, "Cone_BC.diff.RDS")
# saveRDS(MG.diff,      "MG.diff.RDS")
# saveRDS(AC.diff,      "AC.diff.RDS")
# saveRDS(HC.diff,      "HC.diff.RDS")
# saveRDS(RGC.diff,     "RGC.diff.RDS")
# saveRDS(Micro.diff,   "Micro.diff.RDS")
# saveRDS(VEC.diff,     "VEC.diff.RDS")
# saveRDS(Pericyte.diff,"Pericyte.diff.RDS")

Rod.diff      <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/Rod.diff.RDS")
Cone.diff     <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/Cone.diff.RDS")
Rod_BC.diff   <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/Rod_BC.diff.RDS")
Cone_BC.diff  <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/Cone_BC.diff.RDS")
MG.diff       <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/MG.diff.RDS")
AC.diff       <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/AC.diff.RDS")
HC.diff       <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/HC.diff.RDS")
RGC.diff      <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/RGC.diff.RDS")
Micro.diff    <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/Micro.diff.RDS")
VEC.diff      <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/VEC.diff.RDS")
Pericyte.diff <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/Pericyte.diff.RDS")

a <- FindAllMarkers(DMmouse_CRE_temp,logfc.threshold=0,only.pos=F,min.pct=0)

hm_up <- rbind(data.frame(gene=Rod.diff[Rod.diff$p_val_adj<=0.05 & Rod.diff$avg_log2FC>=0.2,] %>% rownames(),               celltype="Rod"),
               data.frame(gene=Cone.diff[Cone.diff$p_val_adj<=0.05 & Cone.diff$avg_log2FC>=0.2,] %>% rownames(),            celltype="Cone"),
               data.frame(gene=Rod_BC.diff[Rod_BC.diff$p_val_adj<=0.05 & Rod_BC.diff$avg_log2FC>=0.2,] %>% rownames(),      celltype="Rod_BC"),
               data.frame(gene=Cone_BC.diff[Cone_BC.diff$p_val_adj<=0.05 & Cone_BC.diff$avg_log2FC>=0.2,] %>% rownames(),   celltype="Cone_BC"),
               data.frame(gene=MG.diff[MG.diff$p_val_adj<=0.05 & MG.diff$avg_log2FC>=0.2,] %>% rownames(),                  celltype="MG"),
               data.frame(gene=AC.diff[AC.diff$p_val_adj<=0.05 & AC.diff$avg_log2FC>=0.2,] %>% rownames(),                  celltype="AC"),
               # data.frame(gene=HC.diff[HC.diff$p_val_adj<=0.05 & HC.diff$avg_log2FC>=0.2,] %>% rownames(),                  celltype="HC"),
               # data.frame(gene=RGC.diff[RGC.diff$p_val_adj<=0.05 & RGC.diff$avg_log2FC>=0.2,] %>% rownames(),               celltype="RGC"),
               data.frame(gene=Micro.diff[Micro.diff$p_val_adj<=0.05 & Micro.diff$avg_log2FC>=0.2,] %>% rownames(),         celltype="Micro")
               # data.frame(gene=VEC.diff[VEC.diff$p_val_adj<=0.05 & VEC.diff$avg_log2FC>=0.2,] %>% rownames(),               celltype="VEC")
               # data.frame(gene=Pericyte.diff[Pericyte.diff$p_val_adj<=0.05 & Pericyte.diff$avg_log2FC>=0.2,] %>% rownames(),celltype="Pericyte")
)



hm_up_new <- data.frame(gene=unique(hm_up$gene))
hm_up_new[,2:8] <- 0
hm_up_new <- hm_up_new[,-1]
rownames(hm_up_new) <- unique(hm_up$gene)
colnames(hm_up_new) <- c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")


for (i in 1:nrow(hm_up)) {
  hm_up_new[rownames(hm_up_new)==hm_up$gene[i],colnames(hm_up_new)==hm_up$celltype[i]] <- 1
}

hm_up_new <-rbind(hm_up_new[rowSums(hm_up_new)==7,],
                  hm_up_new[rowSums(hm_up_new)==6,],
                  hm_up_new[rowSums(hm_up_new)==5,],
                  hm_up_new[rowSums(hm_up_new)==4,],
                  hm_up_new[rowSums(hm_up_new)==3,],
                  hm_up_new[rowSums(hm_up_new)==2,],
                  hm_up_new[rowSums(hm_up_new)==1,])

pheatmap(hm_up_new,
         color = colorRampPalette(c("grey", "red"))(5),cluster_cols = F,cluster_rows = F,show_rownames = F)







hm_dn <- rbind(data.frame(gene=Rod.diff[Rod.diff$p_val_adj<=0.05 & Rod.diff$avg_log2FC<=(-0.2),] %>% rownames(),               celltype="Rod"),
               data.frame(gene=Cone.diff[Cone.diff$p_val_adj<=0.05 & Cone.diff$avg_log2FC<=(-0.2),] %>% rownames(),            celltype="Cone"),
               data.frame(gene=Rod_BC.diff[Rod_BC.diff$p_val_adj<=0.05 & Rod_BC.diff$avg_log2FC<=(-0.2),] %>% rownames(),      celltype="Rod_BC"),
               data.frame(gene=Cone_BC.diff[Cone_BC.diff$p_val_adj<=0.05 & Cone_BC.diff$avg_log2FC<=(-0.2),] %>% rownames(),   celltype="Cone_BC"),
               data.frame(gene=MG.diff[MG.diff$p_val_adj<=0.05 & MG.diff$avg_log2FC<=(-0.2),] %>% rownames(),                  celltype="MG"),
               data.frame(gene=AC.diff[AC.diff$p_val_adj<=0.05 & AC.diff$avg_log2FC<=(-0.2),] %>% rownames(),                  celltype="AC"),
               # data.frame(gene=HC.diff[HC.diff$p_val_adj<=0.05 & HC.diff$avg_log2FC<=(-0.2),] %>% rownames(),                  celltype="HC"),
               # data.frame(gene=RGC.diff[RGC.diff$p_val_adj<=0.05 & RGC.diff$avg_log2FC<=(-0.2),] %>% rownames(),               celltype="RGC"),
               data.frame(gene=Micro.diff[Micro.diff$p_val_adj<=0.05 & Micro.diff$avg_log2FC<=(-0.2),] %>% rownames(),         celltype="Micro")
               # data.frame(gene=VEC.diff[VEC.diff$p_val_adj<=0.05 & VEC.diff$avg_log2FC<=(-0.2),] %>% rownames(),               celltype="VEC")
               # data.frame(gene=Pericyte.diff[Pericyte.diff$p_val_adj<=0.05 & Pericyte.diff$avg_log2FC<=(-0.2),] %>% rownames(),celltype="Pericyte")
)



hm_dn_new <- data.frame(gene=unique(hm_dn$gene))
hm_dn_new[,2:8] <- 0
hm_dn_new <- hm_dn_new[,-1]
rownames(hm_dn_new) <- unique(hm_dn$gene)
colnames(hm_dn_new) <- c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")


for (i in 1:nrow(hm_dn)) {
  hm_dn_new[rownames(hm_dn_new)==hm_dn$gene[i],colnames(hm_dn_new)==hm_dn$celltype[i]] <- 1
}

hm_dn_new <-rbind(hm_dn_new[rowSums(hm_dn_new)==7,],
                  hm_dn_new[rowSums(hm_dn_new)==6,],
                  hm_dn_new[rowSums(hm_dn_new)==5,],
                  hm_dn_new[rowSums(hm_dn_new)==4,],
                  hm_dn_new[rowSums(hm_dn_new)==3,],
                  hm_dn_new[rowSums(hm_dn_new)==2,],
                  hm_dn_new[rowSums(hm_dn_new)==1,])

pheatmap(hm_dn_new,
         color = colorRampPalette(c("grey", "blue"))(5),cluster_cols = F,cluster_rows = F,show_rownames = F)



DMmouse <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge/results/combined_cca_40.RDS")

agg <- AggregateExpression(DMmouse,assays = "RNA",slot = "counts",group.by = "multi")[[1]]
agg <- agg[!rowSums(agg)==0,]

agg <- agg[grep("mt-",rownames(agg),invert = T),]

# a <- melt(cor(agg,method = "spearman"))
# a <- a[!a$value==1,]
# a <- a[order(a$value),]
# a <- a[1:nrow(a)%%2==1,]
# 
# a$group <- "CTRL vs DM"
# a$group[intersect(grep("CTRL",a$Var1),grep("CTRL",a$Var2))] <- "CTRL vs CTRL"
# a$group[intersect(grep("DM",a$Var1),grep("DM",a$Var2))] <- "DM vs DM"
# 
# my_comparison <- list(c("CTRL vs DM","CTRL vs CTRL"),c("CTRL vs DM","DM vs DM"),c("CTRL vs CTRL","DM vs DM"))
# 
# ggboxplot(a, x="group", y="value", fill = "group", palette = "jco",xlab="",add = "jitter")+
#   theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
#         axis.ticks.x = element_blank(),
#         legend.title = element_blank(),
#         axis.text.y = element_text(face = "plain",size = 11),
#         axis.title.x = element_text(face = "plain",size = 13),
#         axis.title.y = element_text(face = "plain",size = 13),
#         legend.position=c(10,10))+
#   # labs(y="log10 (RPM+1)")+
#   #labs(y="log10 (absolute gene expression +1)")+
#   stat_compare_means(comparisons = my_comparison,aes(label=p.sign..),label = "p-value", method = "t.test")+
#   guides(fill = "none")



for (k in 1:ncol(agg)) {
  agg[,k] <- log10(agg[,k]/sum(agg[,k])*1000000+1)
  
}

pheatmap(cor(agg,method = "spearman"),
         # scale = "row",
         cluster_cols = T,cellwidth = 40,cellheight = 40,
         cluster_rows = T,display_numbers = T,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         # col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         col = SpatialColors(n = 100),
         # show_rownames = F,
         # show_colnames = F,
         # treeheight_row = T,
         # treeheight_col = T,
         border = F,
         legend = T)

png("cor_sca_CRE.png", width = 7, height = 7, units = 'in', res = 1000)
heatpairs(agg[,c(1,2,3,4,5,6)],main=NULL)
dev.off()



