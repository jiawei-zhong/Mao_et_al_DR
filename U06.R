## generate bw
options(stringsAsFactors = FALSE)
options(scipen = 100)
library(CAGEfightR)

# SRR12443300 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443300/bam_to_ctss/SRR12443300/bed/SRR12443300.UMI_CB.ctss.bed.gz")
# SRR12443300$V5 <- 0
# SRR12443300 <- unique(SRR12443300)
# gc()
# 
# 
# SRR12443301 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443301/bam_to_ctss/SRR12443301/bed/SRR12443301.UMI_CB.ctss.bed.gz")
# SRR12443301$V5 <- 0
# SRR12443301 <- unique(SRR12443301)
# gc()
# 
# 
# SRR12443302 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443302/bam_to_ctss/SRR12443302/bed/SRR12443302.UMI_CB.ctss.bed.gz")
# SRR12443302$V5 <- 0
# SRR12443302 <- unique(SRR12443302)
# gc()
# 
# 
# SRR12443303 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443303/bam_to_ctss/SRR12443303/bed/SRR12443303.UMI_CB.ctss.bed.gz")
# SRR12443303$V5 <- 0
# SRR12443303 <- unique(SRR12443303)
# gc()
# 
# 
# SRR12443304 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443304/bam_to_ctss/SRR12443304/bed/SRR12443304.UMI_CB.ctss.bed.gz")
# SRR12443304$V5 <- 0
# SRR12443304 <- unique(SRR12443304)
# gc()
# 
# 
# SRR12443305 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443305/bam_to_ctss/SRR12443305/bed/SRR12443305.UMI_CB.ctss.bed.gz")
# SRR12443305$V5 <- 0
# SRR12443305 <- unique(SRR12443305)
# gc()
# 
# 
# write.table(SRR12443300, "SRR12443300.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# write.table(SRR12443301, "SRR12443301.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# write.table(SRR12443302, "SRR12443302.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# write.table(SRR12443303, "SRR12443303.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# write.table(SRR12443304, "SRR12443304.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# write.table(SRR12443305, "SRR12443305.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# 
# 
# temp <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443305/bam_to_ctss/SRR12443305/bed/SRR12443305.collapse.ctss.bed.gz")
# sum(temp$V5)

# system("mv SRR12443300.bed ~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443300/bam_to_ctss/SRR12443300/bed/SRR12443300.UMI_CB.ctss.unique.bed")
# system("mv SRR12443301.bed ~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443301/bam_to_ctss/SRR12443301/bed/SRR12443301.UMI_CB.ctss.unique.bed")
# system("mv SRR12443302.bed ~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443302/bam_to_ctss/SRR12443302/bed/SRR12443302.UMI_CB.ctss.unique.bed")
# system("mv SRR12443303.bed ~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443303/bam_to_ctss/SRR12443303/bed/SRR12443303.UMI_CB.ctss.unique.bed")
# system("mv SRR12443304.bed ~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443304/bam_to_ctss/SRR12443304/bed/SRR12443304.UMI_CB.ctss.unique.bed")
# system("mv SRR12443305.bed ~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443305/bam_to_ctss/SRR12443305/bed/SRR12443305.UMI_CB.ctss.unique.bed")
# 
# 
# options(stringsAsFactors = FALSE)
# # options(scipen = 100)
# options(bedtools.path = "~/anaconda3/bin/")
# 
# 
# .libPaths()
# 
# # library(MASS)
# library(BuenColors)
# library(Seurat)
# library(SeuratData)
# library(ggplot2)
# library(patchwork)
# library(dplyr)
# library(stringr)
# library(VennDiagram)
# library(Rmisc)
# library(ggpubr)
# library(bedtoolsr)
# library(systemPipeR)
# 
# 
# SRR12443300 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443300/bam_to_ctss/SRR12443300/bed/SRR12443300.UMI_CB.ctss.unique.bed")
# SRR12443301 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443301/bam_to_ctss/SRR12443301/bed/SRR12443301.UMI_CB.ctss.unique.bed")
# SRR12443302 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443302/bam_to_ctss/SRR12443302/bed/SRR12443302.UMI_CB.ctss.unique.bed")
# SRR12443303 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443303/bam_to_ctss/SRR12443303/bed/SRR12443303.UMI_CB.ctss.unique.bed")
# SRR12443304 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443304/bam_to_ctss/SRR12443304/bed/SRR12443304.UMI_CB.ctss.unique.bed")
# SRR12443305 <- read.table("/GPUFS/gzzoc_yjhu_3/zjw/20220307_DMmouse/SCAFE_outs/SRR12443305/bam_to_ctss/SRR12443305/bed/SRR12443305.UMI_CB.ctss.unique.bed")
# 
# 
# SRR12443300$V4 <- SRR12443300$V4 %>% strsplit(split = "_") %>% lapply(function(x) x[1]) %>% unlist()
# SRR12443301$V4 <- SRR12443301$V4 %>% strsplit(split = "_") %>% lapply(function(x) x[1]) %>% unlist()
# SRR12443302$V4 <- SRR12443302$V4 %>% strsplit(split = "_") %>% lapply(function(x) x[1]) %>% unlist()
# SRR12443303$V4 <- SRR12443303$V4 %>% strsplit(split = "_") %>% lapply(function(x) x[1]) %>% unlist()
# SRR12443304$V4 <- SRR12443304$V4 %>% strsplit(split = "_") %>% lapply(function(x) x[1]) %>% unlist()
# SRR12443305$V4 <- SRR12443305$V4 %>% strsplit(split = "_") %>% lapply(function(x) x[1]) %>% unlist()
# 
# DMmouse_CRE <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/combined_cca_40.RDS")
# a <- DMmouse_CRE@meta.data
# a$barcode <- rownames(a) %>% strsplit(split = "_") %>% lapply(function(x) x[2]) %>% unlist() %>% gsub(pattern = "[.]",replacement = "-")
# 
# 
# 
# 
# 
# SRR12443300_Rod <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="Rod"]]
# SRR12443300_Cone <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="Cone"]]
# SRR12443300_Rod_BC <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="Rod BC"]]
# SRR12443300_Cone_BC <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="Cone BC"]]
# SRR12443300_MG <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="MG"]]
# SRR12443300_AC <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="AC"]]
# SRR12443300_HC <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="HC"]]
# SRR12443300_RGC <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="RGC"]]
# SRR12443300_Micro <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="Micro"]]
# SRR12443300_VEC <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="VEC"]]
# SRR12443300_Pericyte <- SRR12443300[SRR12443300$V4 %in% a$barcode[a$multi=="DM1" & a$celltype=="Pericyte"]]
# 
# 
# 
# SRR12443301_Rod <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="Rod"]]
# SRR12443301_Cone <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="Cone"]]
# SRR12443301_Rod_BC <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="Rod BC"]]
# SRR12443301_Cone_BC <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="Cone BC"]]
# SRR12443301_MG <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="MG"]]
# SRR12443301_AC <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="AC"]]
# SRR12443301_HC <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="HC"]]
# SRR12443301_RGC <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="RGC"]]
# SRR12443301_Micro <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="Micro"]]
# SRR12443301_VEC <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="VEC"]]
# SRR12443301_Pericyte <- SRR12443301[SRR12443301$V4 %in% a$barcode[a$multi=="DM2" & a$celltype=="Pericyte"]]
# 
# 
# 
# SRR12443302_Rod <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="Rod"]]
# SRR12443302_Cone <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="Cone"]]
# SRR12443302_Rod_BC <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="Rod BC"]]
# SRR12443302_Cone_BC <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="Cone BC"]]
# SRR12443302_MG <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="MG"]]
# SRR12443302_AC <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="AC"]]
# SRR12443302_HC <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="HC"]]
# SRR12443302_RGC <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="RGC"]]
# SRR12443302_Micro <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="Micro"]]
# SRR12443302_VEC <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="VEC"]]
# SRR12443302_Pericyte <- SRR12443302[SRR12443302$V4 %in% a$barcode[a$multi=="DM3" & a$celltype=="Pericyte"]]
# 
# 
# 
# SRR12443303_Rod <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="Rod"]]
# SRR12443303_Cone <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="Cone"]]
# SRR12443303_Rod_BC <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="Rod BC"]]
# SRR12443303_Cone_BC <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="Cone BC"]]
# SRR12443303_MG <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="MG"]]
# SRR12443303_AC <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="AC"]]
# SRR12443303_HC <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="HC"]]
# SRR12443303_RGC <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="RGC"]]
# SRR12443303_Micro <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="Micro"]]
# SRR12443303_VEC <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="VEC"]]
# SRR12443303_Pericyte <- SRR12443303[SRR12443303$V4 %in% a$barcode[a$multi=="CTRL1" & a$celltype=="Pericyte"]]
# 
# 
# 
# SRR12443304_Rod <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="Rod"]]
# SRR12443304_Cone <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="Cone"]]
# SRR12443304_Rod_BC <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="Rod BC"]]
# SRR12443304_Cone_BC <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="Cone BC"]]
# SRR12443304_MG <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="MG"]]
# SRR12443304_AC <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="AC"]]
# SRR12443304_HC <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="HC"]]
# SRR12443304_RGC <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="RGC"]]
# SRR12443304_Micro <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="Micro"]]
# SRR12443304_VEC <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="VEC"]]
# SRR12443304_Pericyte <- SRR12443304[SRR12443304$V4 %in% a$barcode[a$multi=="CTRL2" & a$celltype=="Pericyte"]]
# 
# 
# 
# SRR12443305_Rod <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="Rod"]]
# SRR12443305_Cone <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="Cone"]]
# SRR12443305_Rod_BC <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="Rod BC"]]
# SRR12443305_Cone_BC <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="Cone BC"]]
# SRR12443305_MG <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="MG"]]
# SRR12443305_AC <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="AC"]]
# SRR12443305_HC <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="HC"]]
# SRR12443305_RGC <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="RGC"]]
# SRR12443305_Micro <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="Micro"]]
# SRR12443305_VEC <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="VEC"]]
# SRR12443305_Pericyte <- SRR12443305[SRR12443305$V4 %in% a$barcode[a$multi=="CTRL3" & a$celltype=="Pericyte"]]
# 
# 
# DM_Rod <- rbind(SRR12443300_Rod,SRR12443301_Rod,SRR12443302_Rod)
# DM_Cone <- rbind(SRR12443300_Cone,SRR12443301_Cone,SRR12443302_Cone)
# DM_Rod_BC <- rbind(SRR12443300_Rod_BC,SRR12443301_Rod_BC,SRR12443302_Rod_BC)
# DM_Cone_BC <- rbind(SRR12443300_Cone_BC,SRR12443301_Cone_BC,SRR12443302_Cone_BC)
# DM_MG <- rbind(SRR12443300_MG,SRR12443301_MG,SRR12443302_MG)
# DM_AC <- rbind(SRR12443300_AC,SRR12443301_AC,SRR12443302_AC)
# DM_HC <- rbind(SRR12443300_HC,SRR12443301_HC,SRR12443302_HC)
# DM_RGC <- rbind(SRR12443300_RGC,SRR12443301_RGC,SRR12443302_RGC)
# DM_Micro <- rbind(SRR12443300_Micro,SRR12443301_Micro,SRR12443302_Micro)
# DM_VEC <- rbind(SRR12443300_VEC,SRR12443301_VEC,SRR12443302_VEC)
# DM_Pericyte <- rbind(SRR12443300_Pericyte,SRR12443301_Pericyte,SRR12443302_Pericyte)
# 
# 
# 
# 
# CTRL_Rod <- rbind(SRR12443303_Rod,SRR12443304_Rod,SRR12443305_Rod)
# CTRL_Cone <- rbind(SRR12443303_Cone,SRR12443304_Cone,SRR12443305_Cone)
# CTRL_Rod_BC <- rbind(SRR12443303_Rod_BC,SRR12443304_Rod_BC,SRR12443305_Rod_BC)
# CTRL_Cone_BC <- rbind(SRR12443303_Cone_BC,SRR12443304_Cone_BC,SRR12443305_Cone_BC)
# CTRL_MG <- rbind(SRR12443303_MG,SRR12443304_MG,SRR12443305_MG)
# CTRL_AC <- rbind(SRR12443303_AC,SRR12443304_AC,SRR12443305_AC)
# CTRL_HC <- rbind(SRR12443303_HC,SRR12443304_HC,SRR12443305_HC)
# CTRL_RGC <- rbind(SRR12443303_RGC,SRR12443304_RGC,SRR12443305_RGC)
# CTRL_Micro <- rbind(SRR12443303_Micro,SRR12443304_Micro,SRR12443305_Micro)
# CTRL_VEC <- rbind(SRR12443303_VEC,SRR12443304_VEC,SRR12443305_VEC)
# CTRL_Pericyte <- rbind(SRR12443303_Pericyte,SRR12443304_Pericyte,SRR12443305_Pericyte)
# 
# 
# rm(SRR12443300)
# rm(SRR12443301)
# rm(SRR12443302)
# rm(SRR12443303)
# rm(SRR12443304)
# rm(SRR12443305)
# 
# 
# rm(DMmouse_CRE)
# 
# 
# save.image("20221020.RData")




load("~/zjw/20220307_DMmouse/20221020.RData")

write.table(DM_Rod,"~/zjw/20220307_DMmouse/bigwig/DM_Rod.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Cone,"~/zjw/20220307_DMmouse/bigwig/DM_Cone.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Rod_BC,"~/zjw/20220307_DMmouse/bigwig/DM_Rod_BC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Cone_BC,"~/zjw/20220307_DMmouse/bigwig/DM_Cone_BC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_MG,"~/zjw/20220307_DMmouse/bigwig/DM_MG.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_AC,"~/zjw/20220307_DMmouse/bigwig/DM_AC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_HC,"~/zjw/20220307_DMmouse/bigwig/DM_HC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_RGC,"~/zjw/20220307_DMmouse/bigwig/DM_RGC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Micro,"~/zjw/20220307_DMmouse/bigwig/DM_Micro.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_VEC,"~/zjw/20220307_DMmouse/bigwig/DM_VEC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Pericyte,"~/zjw/20220307_DMmouse/bigwig/DM_Pericyte.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM,"~/zjw/20220307_DMmouse/bigwig/DM.bed",quote = F,row.names = F,col.names = F,sep = "\t")


write.table(CTRL_Rod,"~/zjw/20220307_DMmouse/bigwig/CTRL_Rod.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Cone,"~/zjw/20220307_DMmouse/bigwig/CTRL_Cone.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Rod_BC,"~/zjw/20220307_DMmouse/bigwig/CTRL_Rod_BC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Cone_BC,"~/zjw/20220307_DMmouse/bigwig/CTRL_Cone_BC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_MG,"~/zjw/20220307_DMmouse/bigwig/CTRL_MG.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_AC,"~/zjw/20220307_DMmouse/bigwig/CTRL_AC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_HC,"~/zjw/20220307_DMmouse/bigwig/CTRL_HC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_RGC,"~/zjw/20220307_DMmouse/bigwig/CTRL_RGC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Micro,"~/zjw/20220307_DMmouse/bigwig/CTRL_Micro.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_VEC,"~/zjw/20220307_DMmouse/bigwig/CTRL_VEC.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Pericyte,"~/zjw/20220307_DMmouse/bigwig/CTRL_Pericyte.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL,"~/zjw/20220307_DMmouse/bigwig/CTRL.bed",quote = F,row.names = F,col.names = F,sep = "\t")

# 
# mm10 <- SeqinfoForUCSCGenome(genome = "mm10")
# 
# 
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_Rod.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_Rod_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_Rod_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_Cone.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_Cone_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_Cone_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_Rod_BC.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_Rod_BC_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_Rod_BC_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_Cone_BC.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_Cone_BC_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_Cone_BC_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_MG.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_MG_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_MG_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_AC.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_AC_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_AC_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_HC.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_HC_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_HC_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_RGC.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_RGC_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_RGC_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_Micro.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_Micro_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_Micro_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_VEC.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_VEC_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_VEC_minus.bw",genome = mm10)
# convertBED2BigWig(input = "~/zjw/20220307_DMmouse/bigwig/DM_Pericyte.bed",outputPlus = "~/zjw/20220307_DMmouse/bigwig/DM_Pericyte_plus.bw",outputMinus = "~/zjw/20220307_DMmouse/bigwig/DM_Pericyte_minus.bw",genome = mm10)



myCAGEexp_DM_Rod <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_Rod.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_Cone <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_Cone.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_Rod_BC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_Rod_BC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_Cone_BC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_Cone_BC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_MG <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_MG.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_AC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_AC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_HC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_HC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_RGC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_RGC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_Micro <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_Micro.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_VEC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_VEC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM_Pericyte <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM_Pericyte.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_DM <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/DM.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")

myCAGEexp_DM_Rod <- getCTSS(myCAGEexp_DM_Rod,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_Cone <- getCTSS(myCAGEexp_DM_Cone,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_Rod_BC <- getCTSS(myCAGEexp_DM_Rod_BC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_Cone_BC <- getCTSS(myCAGEexp_DM_Cone_BC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_MG <- getCTSS(myCAGEexp_DM_MG,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_AC <- getCTSS(myCAGEexp_DM_AC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_HC <- getCTSS(myCAGEexp_DM_HC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_RGC <- getCTSS(myCAGEexp_DM_RGC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_Micro <- getCTSS(myCAGEexp_DM_Micro,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_VEC <- getCTSS(myCAGEexp_DM_VEC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM_Pericyte <- getCTSS(myCAGEexp_DM_Pericyte,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_DM <- getCTSS(myCAGEexp_DM,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)



DM_Rod_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_Rod),CTSScoordinatesGR(myCAGEexp_DM_Rod))
DM_Rod_bedgraph <- cbind(DM_Rod_bedgraph[,2],DM_Rod_bedgraph[,3]-1,DM_Rod_bedgraph[,3],DM_Rod_bedgraph[,c(1,4)])
DM_Rod_bedgraph[,1] <- as.character(DM_Rod_bedgraph[,1])
DM_Rod_bedgraph <- DM_Rod_bedgraph[order(DM_Rod_bedgraph[,2]),]
DM_Rod_bedgraph <- DM_Rod_bedgraph[order(DM_Rod_bedgraph[,1]),]
DM_Rod_bedgraph_plus <-  DM_Rod_bedgraph[DM_Rod_bedgraph$strand=="+",1:4]
DM_Rod_bedgraph_minus <- DM_Rod_bedgraph[DM_Rod_bedgraph$strand=="-",1:4]
write.table(DM_Rod_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_Rod_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Rod_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_Rod_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_Cone_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_Cone),CTSScoordinatesGR(myCAGEexp_DM_Cone))
DM_Cone_bedgraph <- cbind(DM_Cone_bedgraph[,2],DM_Cone_bedgraph[,3]-1,DM_Cone_bedgraph[,3],DM_Cone_bedgraph[,c(1,4)])
DM_Cone_bedgraph[,1] <- as.character(DM_Cone_bedgraph[,1])
DM_Cone_bedgraph <- DM_Cone_bedgraph[order(DM_Cone_bedgraph[,1]),]
DM_Cone_bedgraph_plus <-  DM_Cone_bedgraph[DM_Cone_bedgraph$strand=="+",1:4]
DM_Cone_bedgraph_minus <- DM_Cone_bedgraph[DM_Cone_bedgraph$strand=="-",1:4]
write.table(DM_Cone_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_Cone_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Cone_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_Cone_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_Rod_BC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_Rod_BC),CTSScoordinatesGR(myCAGEexp_DM_Rod_BC))
DM_Rod_BC_bedgraph <- cbind(DM_Rod_BC_bedgraph[,2],DM_Rod_BC_bedgraph[,3]-1,DM_Rod_BC_bedgraph[,3],DM_Rod_BC_bedgraph[,c(1,4)])
DM_Rod_BC_bedgraph[,1] <- as.character(DM_Rod_BC_bedgraph[,1])
DM_Rod_BC_bedgraph <- DM_Rod_BC_bedgraph[order(DM_Rod_BC_bedgraph[,1]),]
DM_Rod_BC_bedgraph_plus <-  DM_Rod_BC_bedgraph[DM_Rod_BC_bedgraph$strand=="+",1:4]
DM_Rod_BC_bedgraph_minus <- DM_Rod_BC_bedgraph[DM_Rod_BC_bedgraph$strand=="-",1:4]
write.table(DM_Rod_BC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_Rod_BC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Rod_BC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_Rod_BC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_Cone_BC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_Cone_BC),CTSScoordinatesGR(myCAGEexp_DM_Cone_BC))
DM_Cone_BC_bedgraph <- cbind(DM_Cone_BC_bedgraph[,2],DM_Cone_BC_bedgraph[,3]-1,DM_Cone_BC_bedgraph[,3],DM_Cone_BC_bedgraph[,c(1,4)])
DM_Cone_BC_bedgraph[,1] <- as.character(DM_Cone_BC_bedgraph[,1])
DM_Cone_BC_bedgraph <- DM_Cone_BC_bedgraph[order(DM_Cone_BC_bedgraph[,1]),]
DM_Cone_BC_bedgraph_plus <-  DM_Cone_BC_bedgraph[DM_Cone_BC_bedgraph$strand=="+",1:4]
DM_Cone_BC_bedgraph_minus <- DM_Cone_BC_bedgraph[DM_Cone_BC_bedgraph$strand=="-",1:4]
write.table(DM_Cone_BC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_Cone_BC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Cone_BC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_Cone_BC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_MG_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_MG),CTSScoordinatesGR(myCAGEexp_DM_MG))
DM_MG_bedgraph <- cbind(DM_MG_bedgraph[,2],DM_MG_bedgraph[,3]-1,DM_MG_bedgraph[,3],DM_MG_bedgraph[,c(1,4)])
DM_MG_bedgraph[,1] <- as.character(DM_MG_bedgraph[,1])
DM_MG_bedgraph <- DM_MG_bedgraph[order(DM_MG_bedgraph[,1]),]
DM_MG_bedgraph_plus <-  DM_MG_bedgraph[DM_MG_bedgraph$strand=="+",1:4]
DM_MG_bedgraph_minus <- DM_MG_bedgraph[DM_MG_bedgraph$strand=="-",1:4]
write.table(DM_MG_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_MG_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_MG_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_MG_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_AC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_AC),CTSScoordinatesGR(myCAGEexp_DM_AC))
DM_AC_bedgraph <- cbind(DM_AC_bedgraph[,2],DM_AC_bedgraph[,3]-1,DM_AC_bedgraph[,3],DM_AC_bedgraph[,c(1,4)])
DM_AC_bedgraph[,1] <- as.character(DM_AC_bedgraph[,1])
DM_AC_bedgraph <- DM_AC_bedgraph[order(DM_AC_bedgraph[,1]),]
DM_AC_bedgraph_plus <-  DM_AC_bedgraph[DM_AC_bedgraph$strand=="+",1:4]
DM_AC_bedgraph_minus <- DM_AC_bedgraph[DM_AC_bedgraph$strand=="-",1:4]
write.table(DM_AC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_AC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_AC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_AC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_HC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_HC),CTSScoordinatesGR(myCAGEexp_DM_HC))
DM_HC_bedgraph <- cbind(DM_HC_bedgraph[,2],DM_HC_bedgraph[,3]-1,DM_HC_bedgraph[,3],DM_HC_bedgraph[,c(1,4)])
DM_HC_bedgraph[,1] <- as.character(DM_HC_bedgraph[,1])
DM_HC_bedgraph <- DM_HC_bedgraph[order(DM_HC_bedgraph[,1]),]
DM_HC_bedgraph_plus <-  DM_HC_bedgraph[DM_HC_bedgraph$strand=="+",1:4]
DM_HC_bedgraph_minus <- DM_HC_bedgraph[DM_HC_bedgraph$strand=="-",1:4]
write.table(DM_HC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_HC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_HC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_HC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_RGC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_RGC),CTSScoordinatesGR(myCAGEexp_DM_RGC))
DM_RGC_bedgraph <- cbind(DM_RGC_bedgraph[,2],DM_RGC_bedgraph[,3]-1,DM_RGC_bedgraph[,3],DM_RGC_bedgraph[,c(1,4)])
DM_RGC_bedgraph[,1] <- as.character(DM_RGC_bedgraph[,1])
DM_RGC_bedgraph <- DM_RGC_bedgraph[order(DM_RGC_bedgraph[,1]),]
DM_RGC_bedgraph_plus <-  DM_RGC_bedgraph[DM_RGC_bedgraph$strand=="+",1:4]
DM_RGC_bedgraph_minus <- DM_RGC_bedgraph[DM_RGC_bedgraph$strand=="-",1:4]
write.table(DM_RGC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_RGC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_RGC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_RGC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_Micro_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_Micro),CTSScoordinatesGR(myCAGEexp_DM_Micro))
DM_Micro_bedgraph <- cbind(DM_Micro_bedgraph[,2],DM_Micro_bedgraph[,3]-1,DM_Micro_bedgraph[,3],DM_Micro_bedgraph[,c(1,4)])
DM_Micro_bedgraph[,1] <- as.character(DM_Micro_bedgraph[,1])
DM_Micro_bedgraph <- DM_Micro_bedgraph[order(DM_Micro_bedgraph[,1]),]
DM_Micro_bedgraph_plus <-  DM_Micro_bedgraph[DM_Micro_bedgraph$strand=="+",1:4]
DM_Micro_bedgraph_minus <- DM_Micro_bedgraph[DM_Micro_bedgraph$strand=="-",1:4]
write.table(DM_Micro_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_Micro_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Micro_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_Micro_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_VEC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_VEC),CTSScoordinatesGR(myCAGEexp_DM_VEC))
DM_VEC_bedgraph <- cbind(DM_VEC_bedgraph[,2],DM_VEC_bedgraph[,3]-1,DM_VEC_bedgraph[,3],DM_VEC_bedgraph[,c(1,4)])
DM_VEC_bedgraph[,1] <- as.character(DM_VEC_bedgraph[,1])
DM_VEC_bedgraph <- DM_VEC_bedgraph[order(DM_VEC_bedgraph[,1]),]
DM_VEC_bedgraph_plus <-  DM_VEC_bedgraph[DM_VEC_bedgraph$strand=="+",1:4]
DM_VEC_bedgraph_minus <- DM_VEC_bedgraph[DM_VEC_bedgraph$strand=="-",1:4]
write.table(DM_VEC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_VEC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_VEC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_VEC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_Pericyte_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM_Pericyte),CTSScoordinatesGR(myCAGEexp_DM_Pericyte))
DM_Pericyte_bedgraph <- cbind(DM_Pericyte_bedgraph[,2],DM_Pericyte_bedgraph[,3]-1,DM_Pericyte_bedgraph[,3],DM_Pericyte_bedgraph[,c(1,4)])
DM_Pericyte_bedgraph[,1] <- as.character(DM_Pericyte_bedgraph[,1])
DM_Pericyte_bedgraph <- DM_Pericyte_bedgraph[order(DM_Pericyte_bedgraph[,1]),]
DM_Pericyte_bedgraph_plus <-  DM_Pericyte_bedgraph[DM_Pericyte_bedgraph$strand=="+",1:4]
DM_Pericyte_bedgraph_minus <- DM_Pericyte_bedgraph[DM_Pericyte_bedgraph$strand=="-",1:4]
write.table(DM_Pericyte_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_Pericyte_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_Pericyte_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_Pericyte_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

DM_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_DM),CTSScoordinatesGR(myCAGEexp_DM))
DM_bedgraph <- cbind(DM_bedgraph[,2],DM_bedgraph[,3]-1,DM_bedgraph[,3],DM_bedgraph[,c(1,4)])
DM_bedgraph[,1] <- as.character(DM_bedgraph[,1])
DM_bedgraph <- DM_bedgraph[order(DM_bedgraph[,1]),]
DM_bedgraph_plus <-  DM_bedgraph[DM_bedgraph$strand=="+",1:4]
DM_bedgraph_minus <- DM_bedgraph[DM_bedgraph$strand=="-",1:4]
write.table(DM_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/DM_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(DM_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/DM_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")






myCAGEexp_CTRL_Rod <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_Rod.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_Cone <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_Cone.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_Rod_BC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_Rod_BC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_Cone_BC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_Cone_BC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_MG <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_MG.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_AC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_AC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_HC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_HC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_RGC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_RGC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_Micro <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_Micro.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_VEC <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_VEC.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL_Pericyte <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL_Pericyte.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")
myCAGEexp_CTRL <- CAGEexp(inputFiles = "~/zjw/20220307_DMmouse/bigwig/CTRL.bed",genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFilesType = "bed",sampleLabels = "sample")

myCAGEexp_CTRL_Rod <- getCTSS(myCAGEexp_CTRL_Rod,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_Cone <- getCTSS(myCAGEexp_CTRL_Cone,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_Rod_BC <- getCTSS(myCAGEexp_CTRL_Rod_BC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_Cone_BC <- getCTSS(myCAGEexp_CTRL_Cone_BC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_MG <- getCTSS(myCAGEexp_CTRL_MG,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_AC <- getCTSS(myCAGEexp_CTRL_AC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_HC <- getCTSS(myCAGEexp_CTRL_HC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_RGC <- getCTSS(myCAGEexp_CTRL_RGC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_Micro <- getCTSS(myCAGEexp_CTRL_Micro,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_VEC <- getCTSS(myCAGEexp_CTRL_VEC,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL_Pericyte <- getCTSS(myCAGEexp_CTRL_Pericyte,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)
myCAGEexp_CTRL <- getCTSS(myCAGEexp_CTRL,removeFirstG = FALSE,correctSystematicG = FALSE,useMulticore = T,nrCores = 8)


CTRL_Rod_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_Rod),CTSScoordinatesGR(myCAGEexp_CTRL_Rod))
CTRL_Rod_bedgraph <- cbind(CTRL_Rod_bedgraph[,2],CTRL_Rod_bedgraph[,3]-1,CTRL_Rod_bedgraph[,3],CTRL_Rod_bedgraph[,c(1,4)])
CTRL_Rod_bedgraph[,1] <- as.character(CTRL_Rod_bedgraph[,1])
CTRL_Rod_bedgraph <- CTRL_Rod_bedgraph[order(CTRL_Rod_bedgraph[,2]),]
CTRL_Rod_bedgraph <- CTRL_Rod_bedgraph[order(CTRL_Rod_bedgraph[,1]),]
CTRL_Rod_bedgraph_plus <-  CTRL_Rod_bedgraph[CTRL_Rod_bedgraph$strand=="+",1:4]
CTRL_Rod_bedgraph_minus <- CTRL_Rod_bedgraph[CTRL_Rod_bedgraph$strand=="-",1:4]
write.table(CTRL_Rod_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_Rod_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Rod_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_Rod_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_Cone_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_Cone),CTSScoordinatesGR(myCAGEexp_CTRL_Cone))
CTRL_Cone_bedgraph <- cbind(CTRL_Cone_bedgraph[,2],CTRL_Cone_bedgraph[,3]-1,CTRL_Cone_bedgraph[,3],CTRL_Cone_bedgraph[,c(1,4)])
CTRL_Cone_bedgraph[,1] <- as.character(CTRL_Cone_bedgraph[,1])
CTRL_Cone_bedgraph <- CTRL_Cone_bedgraph[order(CTRL_Cone_bedgraph[,1]),]
CTRL_Cone_bedgraph_plus <-  CTRL_Cone_bedgraph[CTRL_Cone_bedgraph$strand=="+",1:4]
CTRL_Cone_bedgraph_minus <- CTRL_Cone_bedgraph[CTRL_Cone_bedgraph$strand=="-",1:4]
write.table(CTRL_Cone_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_Cone_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Cone_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_Cone_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_Rod_BC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_Rod_BC),CTSScoordinatesGR(myCAGEexp_CTRL_Rod_BC))
CTRL_Rod_BC_bedgraph <- cbind(CTRL_Rod_BC_bedgraph[,2],CTRL_Rod_BC_bedgraph[,3]-1,CTRL_Rod_BC_bedgraph[,3],CTRL_Rod_BC_bedgraph[,c(1,4)])
CTRL_Rod_BC_bedgraph[,1] <- as.character(CTRL_Rod_BC_bedgraph[,1])
CTRL_Rod_BC_bedgraph <- CTRL_Rod_BC_bedgraph[order(CTRL_Rod_BC_bedgraph[,1]),]
CTRL_Rod_BC_bedgraph_plus <-  CTRL_Rod_BC_bedgraph[CTRL_Rod_BC_bedgraph$strand=="+",1:4]
CTRL_Rod_BC_bedgraph_minus <- CTRL_Rod_BC_bedgraph[CTRL_Rod_BC_bedgraph$strand=="-",1:4]
write.table(CTRL_Rod_BC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_Rod_BC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Rod_BC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_Rod_BC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_Cone_BC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_Cone_BC),CTSScoordinatesGR(myCAGEexp_CTRL_Cone_BC))
CTRL_Cone_BC_bedgraph <- cbind(CTRL_Cone_BC_bedgraph[,2],CTRL_Cone_BC_bedgraph[,3]-1,CTRL_Cone_BC_bedgraph[,3],CTRL_Cone_BC_bedgraph[,c(1,4)])
CTRL_Cone_BC_bedgraph[,1] <- as.character(CTRL_Cone_BC_bedgraph[,1])
CTRL_Cone_BC_bedgraph <- CTRL_Cone_BC_bedgraph[order(CTRL_Cone_BC_bedgraph[,1]),]
CTRL_Cone_BC_bedgraph_plus <-  CTRL_Cone_BC_bedgraph[CTRL_Cone_BC_bedgraph$strand=="+",1:4]
CTRL_Cone_BC_bedgraph_minus <- CTRL_Cone_BC_bedgraph[CTRL_Cone_BC_bedgraph$strand=="-",1:4]
write.table(CTRL_Cone_BC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_Cone_BC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Cone_BC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_Cone_BC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_MG_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_MG),CTSScoordinatesGR(myCAGEexp_CTRL_MG))
CTRL_MG_bedgraph <- cbind(CTRL_MG_bedgraph[,2],CTRL_MG_bedgraph[,3]-1,CTRL_MG_bedgraph[,3],CTRL_MG_bedgraph[,c(1,4)])
CTRL_MG_bedgraph[,1] <- as.character(CTRL_MG_bedgraph[,1])
CTRL_MG_bedgraph <- CTRL_MG_bedgraph[order(CTRL_MG_bedgraph[,1]),]
CTRL_MG_bedgraph_plus <-  CTRL_MG_bedgraph[CTRL_MG_bedgraph$strand=="+",1:4]
CTRL_MG_bedgraph_minus <- CTRL_MG_bedgraph[CTRL_MG_bedgraph$strand=="-",1:4]
write.table(CTRL_MG_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_MG_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_MG_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_MG_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_AC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_AC),CTSScoordinatesGR(myCAGEexp_CTRL_AC))
CTRL_AC_bedgraph <- cbind(CTRL_AC_bedgraph[,2],CTRL_AC_bedgraph[,3]-1,CTRL_AC_bedgraph[,3],CTRL_AC_bedgraph[,c(1,4)])
CTRL_AC_bedgraph[,1] <- as.character(CTRL_AC_bedgraph[,1])
CTRL_AC_bedgraph <- CTRL_AC_bedgraph[order(CTRL_AC_bedgraph[,1]),]
CTRL_AC_bedgraph_plus <-  CTRL_AC_bedgraph[CTRL_AC_bedgraph$strand=="+",1:4]
CTRL_AC_bedgraph_minus <- CTRL_AC_bedgraph[CTRL_AC_bedgraph$strand=="-",1:4]
write.table(CTRL_AC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_AC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_AC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_AC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_HC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_HC),CTSScoordinatesGR(myCAGEexp_CTRL_HC))
CTRL_HC_bedgraph <- cbind(CTRL_HC_bedgraph[,2],CTRL_HC_bedgraph[,3]-1,CTRL_HC_bedgraph[,3],CTRL_HC_bedgraph[,c(1,4)])
CTRL_HC_bedgraph[,1] <- as.character(CTRL_HC_bedgraph[,1])
CTRL_HC_bedgraph <- CTRL_HC_bedgraph[order(CTRL_HC_bedgraph[,1]),]
CTRL_HC_bedgraph_plus <-  CTRL_HC_bedgraph[CTRL_HC_bedgraph$strand=="+",1:4]
CTRL_HC_bedgraph_minus <- CTRL_HC_bedgraph[CTRL_HC_bedgraph$strand=="-",1:4]
write.table(CTRL_HC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_HC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_HC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_HC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_RGC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_RGC),CTSScoordinatesGR(myCAGEexp_CTRL_RGC))
CTRL_RGC_bedgraph <- cbind(CTRL_RGC_bedgraph[,2],CTRL_RGC_bedgraph[,3]-1,CTRL_RGC_bedgraph[,3],CTRL_RGC_bedgraph[,c(1,4)])
CTRL_RGC_bedgraph[,1] <- as.character(CTRL_RGC_bedgraph[,1])
CTRL_RGC_bedgraph <- CTRL_RGC_bedgraph[order(CTRL_RGC_bedgraph[,1]),]
CTRL_RGC_bedgraph_plus <-  CTRL_RGC_bedgraph[CTRL_RGC_bedgraph$strand=="+",1:4]
CTRL_RGC_bedgraph_minus <- CTRL_RGC_bedgraph[CTRL_RGC_bedgraph$strand=="-",1:4]
write.table(CTRL_RGC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_RGC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_RGC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_RGC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_Micro_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_Micro),CTSScoordinatesGR(myCAGEexp_CTRL_Micro))
CTRL_Micro_bedgraph <- cbind(CTRL_Micro_bedgraph[,2],CTRL_Micro_bedgraph[,3]-1,CTRL_Micro_bedgraph[,3],CTRL_Micro_bedgraph[,c(1,4)])
CTRL_Micro_bedgraph[,1] <- as.character(CTRL_Micro_bedgraph[,1])
CTRL_Micro_bedgraph <- CTRL_Micro_bedgraph[order(CTRL_Micro_bedgraph[,1]),]
CTRL_Micro_bedgraph_plus <-  CTRL_Micro_bedgraph[CTRL_Micro_bedgraph$strand=="+",1:4]
CTRL_Micro_bedgraph_minus <- CTRL_Micro_bedgraph[CTRL_Micro_bedgraph$strand=="-",1:4]
write.table(CTRL_Micro_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_Micro_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Micro_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_Micro_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_VEC_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_VEC),CTSScoordinatesGR(myCAGEexp_CTRL_VEC))
CTRL_VEC_bedgraph <- cbind(CTRL_VEC_bedgraph[,2],CTRL_VEC_bedgraph[,3]-1,CTRL_VEC_bedgraph[,3],CTRL_VEC_bedgraph[,c(1,4)])
CTRL_VEC_bedgraph[,1] <- as.character(CTRL_VEC_bedgraph[,1])
CTRL_VEC_bedgraph <- CTRL_VEC_bedgraph[order(CTRL_VEC_bedgraph[,1]),]
CTRL_VEC_bedgraph_plus <-  CTRL_VEC_bedgraph[CTRL_VEC_bedgraph$strand=="+",1:4]
CTRL_VEC_bedgraph_minus <- CTRL_VEC_bedgraph[CTRL_VEC_bedgraph$strand=="-",1:4]
write.table(CTRL_VEC_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_VEC_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_VEC_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_VEC_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_Pericyte_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL_Pericyte),CTSScoordinatesGR(myCAGEexp_CTRL_Pericyte))
CTRL_Pericyte_bedgraph <- cbind(CTRL_Pericyte_bedgraph[,2],CTRL_Pericyte_bedgraph[,3]-1,CTRL_Pericyte_bedgraph[,3],CTRL_Pericyte_bedgraph[,c(1,4)])
CTRL_Pericyte_bedgraph[,1] <- as.character(CTRL_Pericyte_bedgraph[,1])
CTRL_Pericyte_bedgraph <- CTRL_Pericyte_bedgraph[order(CTRL_Pericyte_bedgraph[,1]),]
CTRL_Pericyte_bedgraph_plus <-  CTRL_Pericyte_bedgraph[CTRL_Pericyte_bedgraph$strand=="+",1:4]
CTRL_Pericyte_bedgraph_minus <- CTRL_Pericyte_bedgraph[CTRL_Pericyte_bedgraph$strand=="-",1:4]
write.table(CTRL_Pericyte_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_Pericyte_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_Pericyte_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_Pericyte_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

CTRL_bedgraph <- data.frame(CTSStagCountDF(myCAGEexp_CTRL),CTSScoordinatesGR(myCAGEexp_CTRL))
CTRL_bedgraph <- cbind(CTRL_bedgraph[,2],CTRL_bedgraph[,3]-1,CTRL_bedgraph[,3],CTRL_bedgraph[,c(1,4)])
CTRL_bedgraph[,1] <- as.character(CTRL_bedgraph[,1])
CTRL_bedgraph <- CTRL_bedgraph[order(CTRL_bedgraph[,1]),]
CTRL_bedgraph_plus <-  CTRL_bedgraph[CTRL_bedgraph$strand=="+",1:4]
CTRL_bedgraph_minus <- CTRL_bedgraph[CTRL_bedgraph$strand=="-",1:4]
write.table(CTRL_bedgraph_plus, "~/zjw/20220307_DMmouse/bigwig/CTRL_plus.bedgraph", quote = F,row.names = F,col.names = F,sep = "\t")
write.table(CTRL_bedgraph_minus,"~/zjw/20220307_DMmouse/bigwig/CTRL_minus.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")

