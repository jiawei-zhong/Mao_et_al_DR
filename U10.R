## Signac
options(stringsAsFactors = FALSE)
options(scipen = 1000)
.libPaths("~/anaconda3/envs/R4.1/lib/R/library")


library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BuenColors)
library(Seurat)
library(SeuratData)
library(Signac)
library(TFBSTools)
library(JASPAR2020)
library(ggseqlogo)
library(ggpubr)
library(cicero)
library(monocle3)
library(SeuratWrappers)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
library(systemPipeR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)


call <- data.frame(msigdbr())
call <- call[call$gs_cat %in% c("H","C2","C3","C5"),c(3,4,15)]


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

mk <- rbind(data.frame(Rod.diff,     gene = rownames(Rod.diff),cluster = "Rod"),
            data.frame(Cone.diff,    gene = rownames(Rod.diff),cluster = "Cone"),
            data.frame(Rod_BC.diff,  gene = rownames(Rod.diff),cluster = "Rod BC"),
            data.frame(Cone_BC.diff, gene = rownames(Rod.diff),cluster = "Cone BC"),
            data.frame(MG.diff,      gene = rownames(Rod.diff),cluster = "MG"),
            data.frame(AC.diff,      gene = rownames(Rod.diff),cluster = "AC"),
            # data.frame(HC.diff,      gene = rownames(Rod.diff),cluster = "HC"),
            # data.frame(RGC.diff,     gene = rownames(Rod.diff),cluster = "RGC"),
            data.frame(Micro.diff,   gene = rownames(Rod.diff),cluster = "Micro")
            # data.frame(VEC.diff,     gene = rownames(Rod.diff),cluster = "VEC"),
            # data.frame(Pericyte.diff,gene = rownames(Rod.diff),cluster = "Pericyte")
            )

mk$cluster <- factor(mk$cluster,levels = c("Rod","Cone","Rod BC","Cone BC","MG","AC","HC","RGC","Micro","VEC","Pericyte"))

mk$peak_id <- "1"

for (i in 1:nrow(mk)) {
  mk$peak_id[i] <- paste0(CRE$V1[CRE$V4==mk$gene[i]],"-",CRE$V2[CRE$V4==mk$gene[i]],"-",CRE$V3[CRE$V4==mk$gene[i]])
  # print(i)
}

mk.sub <- mk[mk$p_val_adj<=0.05 & abs(mk$avg_log2FC)>=0.2,]



all_term_up <- data.frame()
for (i in unique(mk.sub$cluster)) {
  temp <- mk.sub[mk.sub$cluster==i,]
  temp <- temp[temp$avg_log2FC>0,]
  
  temp <- temp$gene %>% strsplit(",") %>% unlist() %>% strsplit("@") %>% lapply(function(x) x[2]) %>% unlist() %>% grep(pattern = "ADDG",value = T,invert = T) %>% unique()
  go_onto <- enrichGO(temp, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  
  # go_onto <- enricher(toupper(temp), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
  
  
  #go_onto <- simplify(go_onto)
  a <- data.frame(go_onto)
  all_term_up <- rbind(all_term_up,
                       data.frame(a,cluster=i))
}



all_term_dn <- data.frame()
for (i in unique(mk.sub$cluster)) {
  temp <- mk.sub[mk.sub$cluster==i,]
  temp <- temp[temp$avg_log2FC<0,]
  
  temp <- temp$gene %>% strsplit(",") %>% unlist() %>% strsplit("@") %>% lapply(function(x) x[2]) %>% unlist() %>% grep(pattern = "ADDG",value = T,invert = T) %>% unique()
  go_onto <- enrichGO(temp, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  
  # go_onto <- enricher(toupper(temp), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
  
  
  a <- data.frame(go_onto)
  all_term_dn <- rbind(all_term_dn,
                       data.frame(a,cluster=i))
}


# df_up <- data.frame()
# for (i in ls() %>% grep(pattern = "^go",value = T) %>% grep(pattern = "up_BP$",value = T)) {
#   a <- paste0('temp <- data.frame(',i,')')
#   eval(parse(text=a))
#   print(a)
#   df_up <- rbind(df_up,
#                  data.frame(temp,celltype=(i %>% strsplit("_") %>% unlist())[2]))
# }
# 
# df_up$celltype <- factor(df_up$celltype,levels = c("LSC","LEC","ETAC","CTAC","CEC","MC"))




df_up_new <- all_term_up[all_term_up$Description %in% c("cytoplasmic translation", "intrinsic apoptotic signaling pathway", "non-membrane-bounded organelle assembly", "positive regulation of signal transduction by p53 class mediator", "protein stabilization"),]











df_up_new$bubblesize <- df_up_new$Count/max(df_up_new$Count)

df_up_new$log <- -log10(df_up_new$p.adjust)

df_up_new$cluster <- factor(df_up_new$cluster,levels = c("Rod","Cone","Rod BC","Cone BC","MG","AC","HC","RGC","Micro","VEC","Pericyte"))

ggplot(data=df_up_new, mapping=aes(x=cluster,y=Description))+
  geom_point(stat= "identity",aes(size=log),alpha=0.7,show.legend = TRUE)+
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x=element_blank(),y=element_blank(),col="-log10 adjusted p-value") 
  # scale_color_gradientn(colors = c(colorRampPalette(as.character(jdb_palette("brewer_heat")))(250))[126:250])




df_dn_new <- all_term_dn[all_term_dn$Description %in% c("energy derivation by oxidation of organic compounds", "visual perception", "glucose metabolic process", "regulation of protein stability", "retina homeostasis"),]












df_dn_new$bubblesize <- df_dn_new$Count/max(df_dn_new$Count)

df_dn_new$log <- -log10(df_dn_new$p.adjust)

df_dn_new$cluster <- factor(df_dn_new$cluster,levels = c("Rod","Cone","Rod BC","Cone BC","MG","AC","HC","RGC","Micro","VEC","Pericyte"))

ggplot(data=df_dn_new, mapping=aes(x=cluster,y=Description))+
  geom_point(stat= "identity",aes(size=log),alpha=0.7,show.legend = TRUE)+
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x=element_blank(),y=element_blank(),col="-log10 adjusted p-value") +
  scale_color_manual(values = "red")
# scale_color_gradientn(colors = c(colorRampPalette(as.character(jdb_palette("brewer_heat")))(250))[126:250])

  



DMmouse_ATAC.new <- subset(x = DMmouse_ATAC, idents = c("HC", "RGC", "VEC","Pericyte"), invert = TRUE)

DMmouse_ATAC.new$celltype_new <- paste0(DMmouse_ATAC.new$celltype,"_",DMmouse_ATAC.new$var)

Idents(DMmouse_ATAC.new) <- "celltype_new"


open.peaks <- AccessiblePeaks(DMmouse_ATAC.new, idents = levels(Idents(DMmouse_ATAC.new)))
meta.feature <- GetAssayData(DMmouse_ATAC.new, assay = "peaks", slot = "meta.features")






# enriched.motifs_up <- data.frame()
# for (i in levels(mk.sub$cluster)[levels(mk.sub$cluster) %in% unique(mk.sub$cluster)]) {
#   
#   peaks.matched <- MatchRegionStats(
#     meta.feature = meta.feature[open.peaks, ],
#     query.feature = meta.feature[mk.sub$peak_id[mk.sub$cluster==i & mk.sub$avg_log2FC>0], ],
#     n = 50000
#   )
#   enriched.motifs_up <- rbind(enriched.motifs_up,
#                               data.frame(FindMotifs(object = DMmouse_ATAC.new,
#                                                     background = peaks.matched,
#                                                     features = mk.sub$peak_id[mk.sub$cluster==i & mk.sub$avg_log2FC>0]),
#                                          cluster=i))
#   
# }
# 
# 
# enriched.motifs_up_new <- enriched.motifs_up[enriched.motifs_up$p.adjust<0.05 & enriched.motifs_up$fold.enrichment>1.3,]
# 
# enriched.motifs_up_new$cluster <- factor(enriched.motifs_up_new$cluster,levels = unique(enriched.motifs_up_new$cluster))
# 
# 
# 
# 
# hm_up <- enriched.motifs_up_new[,c("motif","cluster")]
# 
# 
# 
# hm_up_new <- data.frame(gene=unique(hm_up$motif))
# hm_up_new[,2:8] <- 0
# hm_up_new <- hm_up_new[,-1]
# rownames(hm_up_new) <- unique(hm_up$motif)
# colnames(hm_up_new) <- c("Rod","Cone","Rod BC","Cone BC","MG","AC","Micro")
# 
#  for (i in 1:nrow(hm_up)) {
#   hm_up_new[rownames(hm_up_new)==hm_up$motif[i],colnames(hm_up_new)==hm_up$cluster[i]] <- 1
# }
# 
# hm_up_new <-rbind(hm_up_new[rowSums(hm_up_new)==7,],
#                   hm_up_new[rowSums(hm_up_new)==6,],
#                   hm_up_new[rowSums(hm_up_new)==5,],
#                   hm_up_new[rowSums(hm_up_new)==4,],
#                   hm_up_new[rowSums(hm_up_new)==3,],
#                   hm_up_new[rowSums(hm_up_new)==2,],
#                   hm_up_new[rowSums(hm_up_new)==1,])
# 
# pheatmap(hm_up_new,
#          color = colorRampPalette(c("grey", "red"))(5),cluster_cols = F,cluster_rows = F,show_rownames = F)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# enriched.motifs_dn <- data.frame()
# for (i in levels(mk.sub$cluster)[levels(mk.sub$cluster) %in% unique(mk.sub$cluster)]) {
#   
#   peaks.matched <- MatchRegionStats(
#     meta.feature = meta.feature[open.peaks, ],
#     query.feature = meta.feature[mk.sub$peak_id[mk.sub$cluster==i & mk.sub$avg_log2FC<0], ],
#     n = 50000
#   )
#   enriched.motifs_dn <- rbind(enriched.motifs_dn,
#                               data.frame(FindMotifs(object = DMmouse_ATAC.new,
#                                                     background = peaks.matched,
#                                                     features = mk.sub$peak_id[mk.sub$cluster==i & mk.sub$avg_log2FC<0]),
#                                          cluster=i))
#   
# }
# 
# 
# enriched.motifs_dn_new <- enriched.motifs_dn[enriched.motifs_dn$p.adjust<0.05 & enriched.motifs_dn$fold.enrichment>1.5,]
# 
# enriched.motifs_dn_new$cluster <- factor(enriched.motifs_dn_new$cluster,levels = unique(enriched.motifs_dn_new$cluster))
# 
# 
# 
# 
# hm_dn <- enriched.motifs_dn_new[,c("motif","cluster")]
# 
# 
# 
# hm_dn_new <- data.frame(gene=unique(hm_dn$motif))
# hm_dn_new[,2:8] <- 0
# hm_dn_new <- hm_dn_new[,-1]
# rownames(hm_dn_new) <- unique(hm_dn$motif)
# colnames(hm_dn_new) <- c("Rod","Cone","Rod BC","Cone BC","MG","AC","Micro")
# 
# 
# for (i in 1:nrow(hm_dn)) {
#   hm_dn_new[rownames(hm_dn_new)==hm_dn$motif[i],colnames(hm_dn_new)==hm_dn$cluster[i]] <- 1
# }
# 
# hm_dn_new <-rbind(hm_dn_new[rowSums(hm_dn_new)==7,],
#                   hm_dn_new[rowSums(hm_dn_new)==6,],
#                   hm_dn_new[rowSums(hm_dn_new)==5,],
#                   hm_dn_new[rowSums(hm_dn_new)==4,],
#                   hm_dn_new[rowSums(hm_dn_new)==3,],
#                   hm_dn_new[rowSums(hm_dn_new)==2,],
#                   hm_dn_new[rowSums(hm_dn_new)==1,])
# 
# pheatmap(hm_dn_new,
#          color = colorRampPalette(c("grey", "blue"))(5),cluster_cols = F,cluster_rows = F,show_rownames = F)














t.test_motif <- data.frame()
for (celltype in c("Rod", "Cone", "Rod BC" ,"Cone BC", "MG", "AC", "Micro")) {
  temp <- data.frame()
  for (i in 1:nrow(DMmouse_ATAC.new@assays$chromvar@data)) {
    
    
    
    barcode1 <- rownames(DMmouse_ATAC@meta.data)[DMmouse_ATAC@meta.data$var=="DM" & DMmouse_ATAC@meta.data$celltype==celltype]
    barcode2 <- rownames(DMmouse_ATAC@meta.data)[DMmouse_ATAC@meta.data$var=="CTRL" & DMmouse_ATAC@meta.data$celltype==celltype]
    
    temp <- rbind(temp,
                  data.frame(p.val=t.test(DMmouse_ATAC.new@assays$chromvar@data[i,barcode1],DMmouse_ATAC.new@assays$chromvar@data[i,barcode2])$p.value,
                             ave_delta=mean(DMmouse_ATAC.new@assays$chromvar@data[i,barcode1])-mean(DMmouse_ATAC.new@assays$chromvar@data[i,barcode2]),
                             mean_1=mean(DMmouse_ATAC.new@assays$chromvar@data[i,barcode1]),
                             mean_2=mean(DMmouse_ATAC.new@assays$chromvar@data[i,barcode2]),
                             cluster=celltype,
                             motif=rownames(DMmouse_ATAC.new@assays$chromvar@data)[i],
                             motif.name=enriched.motifs$motif.name[enriched.motifs$motif==rownames(DMmouse_ATAC.new@assays$chromvar@data)[i]][1]
                  )
    )
  }
  temp <- temp[order(temp$ave_delta,decreasing = T),]
  temp <- temp[order(temp$p.val,decreasing = F),]
  t.test_motif <- rbind(t.test_motif,temp)
}



t.test_motif$adj.P.Val <- p.adjust(t.test_motif$p.val,method = "BH")

t.test_motif$cluster <- factor(t.test_motif$cluster,levels = c("Rod","Cone","Rod BC","Cone BC","MG","AC","Micro"))

t.test_motif_DE <- t.test_motif[t.test_motif$p.val<0.05 & abs(t.test_motif$ave_delta)>0.2,]






hm_up <- t.test_motif_DE[t.test_motif_DE$ave_delta>0,c("motif","cluster")]



hm_up_new <- data.frame(gene=unique(hm_up$motif))
hm_up_new[,2:8] <- 0
hm_up_new <- hm_up_new[,-1]
rownames(hm_up_new) <- unique(hm_up$motif)
colnames(hm_up_new) <- c("Rod","Cone","Rod BC","Cone BC","MG","AC","Micro")


for (i in 1:nrow(hm_up)) {
  hm_up_new[rownames(hm_up_new)==hm_up$motif[i],colnames(hm_up_new)==hm_up$cluster[i]] <- 1
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











hm_dn <- t.test_motif_DE[t.test_motif_DE$ave_delta<0,c("motif","cluster")]



hm_dn_new <- data.frame(gene=unique(hm_dn$motif))
hm_dn_new[,2:8] <- 0
hm_dn_new <- hm_dn_new[,-1]
rownames(hm_dn_new) <- unique(hm_dn$motif)
colnames(hm_dn_new) <- c("Rod","Cone","Rod BC","Cone BC","MG","AC","Micro")


for (i in 1:nrow(hm_dn)) {
  hm_dn_new[rownames(hm_dn_new)==hm_dn$motif[i],colnames(hm_dn_new)==hm_dn$cluster[i]] <- 1
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




go_onto <- enrichGO(t.test_motif_DE$motif.name[t.test_motif_DE$ave_delta>0 & t.test_motif_DE$cluster=="Rod"] %>% gsub(pattern = "[(]var.1[)]",replacement = "") %>% gsub(pattern = "[(]var.2[)]",replacement = "") %>% gsub(pattern = "[(]var.3[)]",replacement = "") %>% strsplit(split="::") %>% unlist() %>% str_to_title(), 
                    OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)





all_term_motif_up <- data.frame()
for (i in unique(mk.sub$cluster)) {
  temp <- t.test_motif_DE$motif.name[t.test_motif_DE$ave_delta>0 & t.test_motif_DE$cluster==i] %>% gsub(pattern = "[(]var.1[)]",replacement = "") %>% gsub(pattern = "[(]var.2[)]",replacement = "") %>% gsub(pattern = "[(]var.3[)]",replacement = "") %>% strsplit(split="::") %>% unlist() %>% str_to_title()

  go_onto <- enrichGO(temp, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  
  # go_onto <- enricher(toupper(temp), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
  
  
  #go_onto <- simplify(go_onto)
  a <- data.frame(go_onto)
  all_term_motif_up <- rbind(all_term_motif_up,
                       data.frame(a,cluster=i))
}


all_term_motif_dn <- data.frame()
for (i in unique(mk.sub$cluster)) {
  temp <- t.test_motif_DE$motif.name[t.test_motif_DE$ave_delta<0 & t.test_motif_DE$cluster==i] %>% gsub(pattern = "[(]var.1[)]",replacement = "") %>% gsub(pattern = "[(]var.2[)]",replacement = "") %>% gsub(pattern = "[(]var.3[)]",replacement = "") %>% strsplit(split="::") %>% unlist() %>% str_to_title()
  
  go_onto <- enrichGO(temp, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  
  # go_onto <- enricher(todnper(temp), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
  
  
  #go_onto <- simplify(go_onto)
  a <- data.frame(go_onto)
  all_term_motif_dn <- rbind(all_term_motif_dn,
                             data.frame(a,cluster=i))
}




## Eno3 Zrsr1

# Zrsr1 Rod变短 Cone BC变长
# Psph  Rod变长
# Snx10 Rod_BC Cone_BC变长
# Mgat1 Cone_BC MG 变长

