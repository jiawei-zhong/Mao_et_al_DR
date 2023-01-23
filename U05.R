### SNN metacell
options(stringsAsFactors = FALSE)
# options(scipen = 100)



.libPaths()

# library(MASS)
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


library(factoextra)


set.seed(19940605)

DMmouse_CRE <- readRDS("~/zjw/20220307DMmouse/Seurat_merge_CRE/results/combined_cca_40.RDS")
DMmouse_CRE.CTRL <- subset(DMmouse_CRE, subset = var == "CTRL")
DMmouse_CRE.DM <- subset(DMmouse_CRE, subset = var == "DM")







# DMmouse_CRE.CTRL.Rod <- subset(DMmouse_CRE.CTRL, subset = celltype == "Rod")
# DMmouse_CRE.CTRL.Rod <- FindNeighbors(object = DMmouse_CRE.CTRL.Rod, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.CTRL.Rod <- kmeans(x = DMmouse_CRE.CTRL.Rod@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.CTRL.Rod)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.CTRL.Rod)
# 
# DMmouse_CRE.CTRL.Cone <- subset(DMmouse_CRE.CTRL, subset = celltype == "Cone")
# DMmouse_CRE.CTRL.Cone <- FindNeighbors(object = DMmouse_CRE.CTRL.Cone, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.CTRL.Cone <- kmeans(x = DMmouse_CRE.CTRL.Cone@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.CTRL.Cone)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.CTRL.Cone)
# 
# DMmouse_CRE.CTRL.Rod_BC <- subset(DMmouse_CRE.CTRL, subset = celltype == "Rod BC")
# DMmouse_CRE.CTRL.Rod_BC <- FindNeighbors(object = DMmouse_CRE.CTRL.Rod_BC, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.CTRL.Rod_BC <- kmeans(x = DMmouse_CRE.CTRL.Rod_BC@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.CTRL.Rod_BC)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.CTRL.Rod_BC)
# 
# DMmouse_CRE.CTRL.Cone_BC <- subset(DMmouse_CRE.CTRL, subset = celltype == "Cone BC")
# DMmouse_CRE.CTRL.Cone_BC <- FindNeighbors(object = DMmouse_CRE.CTRL.Cone_BC, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.CTRL.Cone_BC <- kmeans(x = DMmouse_CRE.CTRL.Cone_BC@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.CTRL.Cone_BC)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.CTRL.Cone_BC)
# 
# DMmouse_CRE.CTRL.MG <- subset(DMmouse_CRE.CTRL, subset = celltype == "MG")
# DMmouse_CRE.CTRL.MG <- FindNeighbors(object = DMmouse_CRE.CTRL.MG, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.CTRL.MG <- kmeans(x = DMmouse_CRE.CTRL.MG@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.CTRL.MG)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.CTRL.MG)
# 
# DMmouse_CRE.CTRL.AC <- subset(DMmouse_CRE.CTRL, subset = celltype == "AC")
# DMmouse_CRE.CTRL.AC <- FindNeighbors(object = DMmouse_CRE.CTRL.AC, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.CTRL.AC <- kmeans(x = DMmouse_CRE.CTRL.AC@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.CTRL.AC)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.CTRL.AC)
# 
# DMmouse_CRE.CTRL.Micro <- subset(DMmouse_CRE.CTRL, subset = celltype == "Micro")
# DMmouse_CRE.CTRL.Micro <- FindNeighbors(object = DMmouse_CRE.CTRL.Micro, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.CTRL.Micro <- kmeans(x = DMmouse_CRE.CTRL.Micro@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.CTRL.Micro)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.CTRL.Micro)
# 
# 
# 
# 
# 
# 
# 
# DMmouse_CRE.DM.Rod <- subset(DMmouse_CRE.DM, subset = celltype == "Rod")
# DMmouse_CRE.DM.Rod <- FindNeighbors(object = DMmouse_CRE.DM.Rod, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.DM.Rod <- kmeans(x = DMmouse_CRE.DM.Rod@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.DM.Rod)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.DM.Rod)
# 
# DMmouse_CRE.DM.Cone <- subset(DMmouse_CRE.DM, subset = celltype == "Cone")
# DMmouse_CRE.DM.Cone <- FindNeighbors(object = DMmouse_CRE.DM.Cone, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.DM.Cone <- kmeans(x = DMmouse_CRE.DM.Cone@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.DM.Cone)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.DM.Cone)
# 
# DMmouse_CRE.DM.Rod_BC <- subset(DMmouse_CRE.DM, subset = celltype == "Rod BC")
# DMmouse_CRE.DM.Rod_BC <- FindNeighbors(object = DMmouse_CRE.DM.Rod_BC, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.DM.Rod_BC <- kmeans(x = DMmouse_CRE.DM.Rod_BC@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.DM.Rod_BC)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.DM.Rod_BC)
# 
# DMmouse_CRE.DM.Cone_BC <- subset(DMmouse_CRE.DM, subset = celltype == "Cone BC")
# DMmouse_CRE.DM.Cone_BC <- FindNeighbors(object = DMmouse_CRE.DM.Cone_BC, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.DM.Cone_BC <- kmeans(x = DMmouse_CRE.DM.Cone_BC@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.DM.Cone_BC)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.DM.Cone_BC)
# 
# DMmouse_CRE.DM.MG <- subset(DMmouse_CRE.DM, subset = celltype == "MG")
# DMmouse_CRE.DM.MG <- FindNeighbors(object = DMmouse_CRE.DM.MG, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.DM.MG <- kmeans(x = DMmouse_CRE.DM.MG@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.DM.MG)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.DM.MG)
# 
# DMmouse_CRE.DM.AC <- subset(DMmouse_CRE.DM, subset = celltype == "AC")
# DMmouse_CRE.DM.AC <- FindNeighbors(object = DMmouse_CRE.DM.AC, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.DM.AC <- kmeans(x = DMmouse_CRE.DM.AC@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.DM.AC)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.DM.AC)
# 
# DMmouse_CRE.DM.Micro <- subset(DMmouse_CRE.DM, subset = celltype == "Micro")
# DMmouse_CRE.DM.Micro <- FindNeighbors(object = DMmouse_CRE.DM.Micro, reduction = "pca", dims = 1:30,compute.SNN = T)
# km.res.DM.Micro <- kmeans(x = DMmouse_CRE.DM.Micro@graphs$integrated_snn, centers = ceiling(dim(DMmouse_CRE.DM.Micro)[2]/200), nstart = 25)
# print("finish")
# gc()
# rm(DMmouse_CRE.DM.Micro)
# 
# for (i in grep("DMmouse",ls(),value = T)) {
#   eval(parse(text=paste0('rm(',i,')')))
# }
# 
# save.image("cluster_20221005.RData")



km.res.CTRL.Rod     <- readRDS("Seurat_merge_CRE/results/km.res.CTRL.Rod.RDS")
km.res.CTRL.Cone    <- readRDS("Seurat_merge_CRE/results/km.res.CTRL.Cone.RDS")
km.res.CTRL.Rod_BC  <- readRDS("Seurat_merge_CRE/results/km.res.CTRL.Rod_BC.RDS")
km.res.CTRL.Cone_BC <- readRDS("Seurat_merge_CRE/results/km.res.CTRL.Cone_BC.RDS")
km.res.CTRL.MG      <- readRDS("Seurat_merge_CRE/results/km.res.CTRL.MG.RDS")
km.res.CTRL.AC      <- readRDS("Seurat_merge_CRE/results/km.res.CTRL.AC.RDS")
km.res.CTRL.Micro   <- readRDS("Seurat_merge_CRE/results/km.res.CTRL.Micro.RDS")


km.res.DM.Rod     <- readRDS("Seurat_merge_CRE/results/km.res.DM.Rod.RDS")
km.res.DM.Cone    <- readRDS("Seurat_merge_CRE/results/km.res.DM.Cone.RDS")
km.res.DM.Rod_BC  <- readRDS("Seurat_merge_CRE/results/km.res.DM.Rod_BC.RDS")
km.res.DM.Cone_BC <- readRDS("Seurat_merge_CRE/results/km.res.DM.Cone_BC.RDS")
km.res.DM.MG      <- readRDS("Seurat_merge_CRE/results/km.res.DM.MG.RDS")
km.res.DM.AC      <- readRDS("Seurat_merge_CRE/results/km.res.DM.AC.RDS")
km.res.DM.Micro   <- readRDS("Seurat_merge_CRE/results/km.res.DM.Micro.RDS")

temp_HC <- c("HC_1",name)


cluster <- rbind(data.frame(barcode=names(km.res.CTRL.Rod$cluster),    cluster=paste0("CTRL.Rod_",    formatC(km.res.CTRL.Rod$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.CTRL.Cone$cluster),   cluster=paste0("CTRL.Cone_",   formatC(km.res.CTRL.Cone$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.CTRL.Rod_BC$cluster), cluster=paste0("CTRL.Rod_BC_", formatC(km.res.CTRL.Rod_BC$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.CTRL.Cone_BC$cluster),cluster=paste0("CTRL.Cone_BC_",formatC(km.res.CTRL.Cone_BC$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.CTRL.MG$cluster),     cluster=paste0("CTRL.MG_",     formatC(km.res.CTRL.MG$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.CTRL.AC$cluster),     cluster=paste0("CTRL.AC_",     formatC(km.res.CTRL.AC$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.CTRL.Micro$cluster),  cluster=paste0("CTRL.Micro_",  formatC(km.res.CTRL.Micro$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(Idents(DMmouse_CRE.CTRL)[Idents(DMmouse_CRE.CTRL)=="HC"]),       cluster=paste0("CTRL.",Idents(DMmouse_CRE.CTRL)[Idents(DMmouse_CRE.CTRL)=="HC"],"_001")),
                 data.frame(barcode=names(Idents(DMmouse_CRE.CTRL)[Idents(DMmouse_CRE.CTRL)=="VEC"]),      cluster=paste0("CTRL.",Idents(DMmouse_CRE.CTRL)[Idents(DMmouse_CRE.CTRL)=="VEC"],"_001")),
                 data.frame(barcode=names(Idents(DMmouse_CRE.CTRL)[Idents(DMmouse_CRE.CTRL)=="Pericyte"]), cluster=paste0("CTRL.",Idents(DMmouse_CRE.CTRL)[Idents(DMmouse_CRE.CTRL)=="Pericyte"],"_001")),
                 data.frame(barcode=names(km.res.DM.Rod$cluster),    cluster=paste0("DM.Rod_",    formatC(km.res.DM.Rod$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.DM.Cone$cluster),   cluster=paste0("DM.Cone_",   formatC(km.res.DM.Cone$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.DM.Rod_BC$cluster), cluster=paste0("DM.Rod_BC_", formatC(km.res.DM.Rod_BC$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.DM.Cone_BC$cluster),cluster=paste0("DM.Cone_BC_",formatC(km.res.DM.Cone_BC$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.DM.MG$cluster),     cluster=paste0("DM.MG_",     formatC(km.res.DM.MG$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.DM.AC$cluster),     cluster=paste0("DM.AC_",     formatC(km.res.DM.AC$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(km.res.DM.Micro$cluster),  cluster=paste0("DM.Micro_",  formatC(km.res.DM.Micro$cluster, flag = 0, width = 3))),
                 data.frame(barcode=names(Idents(DMmouse_CRE.DM)[Idents(DMmouse_CRE.DM)=="HC"]),       cluster=paste0("DM.",Idents(DMmouse_CRE.DM)[Idents(DMmouse_CRE.DM)=="HC"],"_001")),
                 data.frame(barcode=names(Idents(DMmouse_CRE.DM)[Idents(DMmouse_CRE.DM)=="VEC"]),      cluster=paste0("DM.",Idents(DMmouse_CRE.DM)[Idents(DMmouse_CRE.DM)=="VEC"],"_001")),
                 data.frame(barcode=names(Idents(DMmouse_CRE.DM)[Idents(DMmouse_CRE.DM)=="Pericyte"]), cluster=paste0("DM.",Idents(DMmouse_CRE.DM)[Idents(DMmouse_CRE.DM)=="Pericyte"],"_001"))
)

metacell.df <- data.frame()

for (i in unique(cluster$cluster)) {
  
  if (i=="CTRL.Rod_028") {
    metacell.df <- cbind(rowSums(DMmouse_CRE@assays$RNA@counts[,colnames(DMmouse_CRE@assays$RNA@counts) %in% cluster$barcode[cluster$cluster==i]]))
  } else {
    metacell.df <- cbind(metacell.df,
                         rowSums(DMmouse_CRE@assays$RNA@counts[,colnames(DMmouse_CRE@assays$RNA@counts) %in% cluster$barcode[cluster$cluster==i]]))
  }
  
}

colnames(metacell.df) <- unique(cluster$cluster)


metacell <- CreateSeuratObject(counts = metacell.df)


metacell$orig.ident <- strsplit(colnames(metacell),split = "_0|_1") %>% lapply(function(x) x[1]) %>% unlist()

metacell <- ScaleData(metacell)
metacell <- NormalizeData(metacell)

metacell <- RunPCA(object = metacell, features = VariableFeatures(object = metacell))
metacell <- RunUMAP(object = metacell, reduction = "pca", dims = 1:30)
metacell <- FindNeighbors(object = metacell, reduction = "pca", dims = 1:30)
metacell <- FindClusters(object = metacell, resolution = 0.1)
Idents(metacell) <- "orig.ident"
DimPlot(metacell)+
  scale_color_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B","CLP","MEP","HSC")))


temp <- metacell@reductions$umap@cell.embeddings

for (i in 1:nrow(temp)) {
  barcode <- cluster$barcode[cluster$cluster==rownames(temp)[i]]
  temp[i,1] <- mean(DMmouse_CRE.CTRL@reductions$umap@cell.embeddings[,1][names(DMmouse_CRE.CTRL@reductions$umap@cell.embeddings[,1]) %in% barcode])
  temp[i,2] <- mean(DMmouse_CRE.CTRL@reductions$umap@cell.embeddings[,2][names(DMmouse_CRE.CTRL@reductions$umap@cell.embeddings[,2]) %in% barcode])
}

metacell@reductions$umap@cell.embeddings <- temp

