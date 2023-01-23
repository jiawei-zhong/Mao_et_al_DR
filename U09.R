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



DMmouse_CRE <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/combined_cca_40.RDS")
genes_remove <- c(rownames(DMmouse_CRE@assays$RNA@counts)[rowSums(DMmouse_CRE@assays$RNA@counts)==0],"P:r1@ADDG12084996058.R")
DMmouse_CRE <- subset(DMmouse_CRE, features = rownames(DMmouse_CRE@assays$RNA@counts)[!rownames(DMmouse_CRE@assays$RNA@counts) %in% genes_remove])

CRE <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/annotate/DM/bed/DM.CRE.annot.bed.gz",sep="\t")
CRE <- CRE[!CRE$V4 %in% genes_remove,]

mk.CTRL <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/CTRL.marker.RDS")
mk.CTRL <- mk.CTRL[mk.CTRL$avg_log2FC>0.2 & mk.CTRL$p_val_adj<0.05,]
mk.CTRL$peak_id <- "1"

for (i in 1:nrow(mk.CTRL)) {
  mk.CTRL$peak_id[i] <- paste0(CRE$V1[CRE$V4==mk.CTRL$gene[i]],"-",CRE$V2[CRE$V4==mk.CTRL$gene[i]],"-",CRE$V3[CRE$V4==mk.CTRL$gene[i]])
  # print(i)
}

# counts <- Read10X_h5(filename = "Signac/atac_v1_DMmouse_ATAC_10k_filtered_peak_bc_matrix.h5")

counts <- DMmouse_CRE@assays$RNA@counts

rownames(counts) <- paste0(CRE$V1,":",CRE$V2,"-",CRE$V3)


# metadata <- read.csv(
#   file = "Signac/atac_v1_DMmouse_ATAC_10k_singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'snapATAC/output_sorted_new.bed.gz',
  min.cells = 0,
  min.features = 0
)

DMmouse_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  # meta.data = metadata
)

DMmouse_ATAC@reductions$umap <- DMmouse_CRE@reductions$umap

Idents(DMmouse_ATAC) <- Idents(DMmouse_CRE)

DMmouse_ATAC$celltype <- DMmouse_CRE$celltype
DMmouse_ATAC$var <- DMmouse_CRE$var

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(DMmouse_ATAC) <- annotations


DMmouse_ATAC <- NucleosomeSignal(object = DMmouse_ATAC)


DMmouse_ATAC$nucleosome_group <- ifelse(DMmouse_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = DMmouse_ATAC, group.by = 'nucleosome_group', region = 'chr1-1-10000000')


DMmouse_ATAC <- TSSEnrichment(DMmouse_ATAC, fast = FALSE)
DMmouse_ATAC$high.tss <- ifelse(DMmouse_ATAC$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(DMmouse_ATAC, group.by = 'high.tss') + NoLegend()


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))



# add motif information
DMmouse_ATAC <- AddMotifs(object = DMmouse_ATAC, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)


DMmouse_ATAC <- RunChromVAR(object = DMmouse_ATAC, genome = BSgenome.Mmusculus.UCSC.mm10)


DMmouse_ATAC.CTRL <- subset(DMmouse_ATAC, subset = var == "CTRL")




open.peaks <- AccessiblePeaks(DMmouse_ATAC.CTRL, idents = levels(mk.CTRL$cluster))
meta.feature <- GetAssayData(DMmouse_ATAC.CTRL, assay = "peaks", slot = "meta.features")

enriched.motifs <- data.frame()
for (i in levels(mk.CTRL$cluster)) {
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[mk.CTRL$peak_id[mk.CTRL$cluster==i], ],
    n = 50000
  )
  enriched.motifs <- rbind(enriched.motifs,
                           data.frame(FindMotifs(object = DMmouse_ATAC.CTRL,
                                                 background = peaks.matched,
                                                 features = mk.CTRL$peak_id[mk.CTRL$cluster==i]),
                                      cluster=i))
  
}


enriched.motifs_new <- enriched.motifs[enriched.motifs$p.adjust<0.05 & enriched.motifs$fold.enrichment>1.3,]

enriched.motifs_new$cluster <- factor(enriched.motifs_new$cluster,levels = unique(enriched.motifs_new$cluster))


write.table(enriched.motifs_new$motif.name[enriched.motifs_new$cluster=="Micro"],"temp.txt",quote = F,row.names = F,col.names = F)

a <- read.table("MG_string_interactions_short.tsv",header = 1,sep = "\t")

c(a[a$combined_score>0.99,1],a[a$combined_score>0.99,2]) %>% unique() %>% length()

DefaultAssay(DMmouse_ATAC.CTRL) <- 'chromvar'

ave_exp <- AverageExpression(DMmouse_ATAC.CTRL,assays = "chromvar")[[1]]

gap <- table(enriched.motifs_new$cluster)
for (i in 2:length(gap)) {
  gap[i] <- gap[i]+gap[i-1]
}



df <- data.frame()
for (i in 1:nrow(enriched.motifs_new)) {
  temp <- t(data.frame(ave_exp[rownames(ave_exp)==enriched.motifs_new$motif[i],]))
  rownames(temp) <- paste0(enriched.motifs_new$motif[i],"_",enriched.motifs_new$cluster[i])
  df <- rbind(df,temp)
}

colnames(df) <- names(ave_exp[1,])
# rownames(df) <- gene_list



library(matrixStats)
df <- (df-rowMeans(df))/(rowSds(as.matrix(df)))[row(df)]



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
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         show_rownames = F,
         # show_colnames = F,
         treeheight_row = F,
         treeheight_col = F,
         border = F,
         gaps_row = gap,
         legend = T)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))


# p1 <- FeaturePlot(object = DMmouse_ATAC.CTRL,features = "MA0712.2",min.cutoff = "q10",max.cutoff = "q90") +
#   scale_color_gradientn(colours = SpatialColors(n = 100)[51:100]) # Otx2 Rod RodBC marker

p1 <- FeaturePlot(object = DMmouse_ATAC.CTRL,features = "MA0467.1",min.cutoff = "q10",max.cutoff = "q90") +
  scale_color_gradientn(colours = SpatialColors(n = 100)[51:100]) +
  NoLegend() +
  NoAxes() # Crx Rod RodBC marker

p2 <- FeaturePlot(object = DMmouse_ATAC.CTRL,features = "MA0100.3",min.cutoff = "q10",max.cutoff = "q90") +
  scale_color_gradientn(colours = SpatialColors(n = 100)[51:100]) +
  NoLegend() +
  NoAxes() ## MYB  bc marker?

p3 <- FeaturePlot(object = DMmouse_ATAC.CTRL,features = "MA1528.1",min.cutoff = "q10",max.cutoff = "q90") +
  scale_color_gradientn(colours = SpatialColors(n = 100)[51:100]) +
  NoLegend() +
  NoAxes() ## NFIX  MG 

p4 <- FeaturePlot(object = DMmouse_ATAC.CTRL,features = "MA0081.2",min.cutoff = "q10",max.cutoff = "q90") +
  scale_color_gradientn(colours = SpatialColors(n = 100)[51:100]) +
  NoLegend() +
  NoAxes() ## SPIB  Micro


# p <- ggarrange(plotlist = list(p1,p2,p3,p4),nrow = 1)

ggsave(filename="Crx_motif.jpeg", plot=p1, width = 5, height = 5.5, units = 'in', dpi = 1000,device = "jpeg")
ggsave(filename="Myb_motif.jpeg", plot=p2, width = 5, height = 5.5, units = 'in', dpi = 1000,device = "jpeg")
ggsave(filename="Nfix_motif.jpeg", plot=p3, width = 5, height = 5.5, units = 'in', dpi = 1000,device = "jpeg")
ggsave(filename="Spib_motif.jpeg", plot=p4, width = 5, height = 5.5, units = 'in', dpi = 1000,device = "jpeg")

p <- FeaturePlot(object = DMmouse_ATAC.CTRL,features = "MA0679.2",min.cutoff = "q10",max.cutoff = "q90") +
  scale_color_gradientn(colours = SpatialColors(n = 100)[51:100]) # Onecut1





# in heatmap show MA0467.1(Crx,rod RodBC)  MA0100.3(myb,bc)  MA1528.1(NFIX,MG) MA0810.1(TFAP2A(var.2),AC HC RGC) MA0081.2(Micro,SPIB) MA0475.2(VEC,FLI1,https://onlinelibrary.wiley.com/doi/full/10.1002/adfm.202203069) 


MotifPlot(object = DMmouse_ATAC.CTRL,motifs = c("MA0467.1","MA0100.3","MA1528.1","MA0810.1","MA0081.2","MA0475.2"),assay = 'peaks')


plot(Venn(list("Rod"=enriched.motifs$motif[enriched.motifs$cluster=="Rod"],
               "RodBC"=enriched.motifs$motif[enriched.motifs$cluster=="RodBC"]))
 ,doWeight=T)

# venn.diagram(x = list("Rod"=enriched.motifs_new$motif[enriched.motifs_new$cluster=="Rod"],"VEC"=enriched.motifs_new$motif[enriched.motifs_new$cluster=="VEC"]),
#              imagetype = "png",
#              hyper.test = T,
#              total.population = length(unique(enriched.motifs$motif)),
#              lower.tail = F,
#              lwd = 4,
#              fill = c("cornflowerblue", "darkorchid1"),
#              alpha = 0.75,
#              label.col = "black",
#              cex = 2,
#              fontfamily = "serif",
#              fontface = "bold",
#              cat.col = c("cornflowerblue", "darkorchid1"),
#              cat.cex = 2,
#              cat.fontfamily = "serif",
#              cat.fontface = "bold",
#              cat.dist = c(0.03, 0.03),
#              cat.pos = c(-20, 14),
#              filename = "temp.png")


calculate.phyper.p.value <- function (list1, list2, total.size, lower.tail = TRUE, adjust = FALSE) 
{
  actual.overlap <- length(intersect(list1, list2))
  expected.overlap <- as.numeric(length(list1)) * length(list2)/total.size
  adjust.value <- 0
  if (adjust & !lower.tail) {
    adjust.value <- 1
    warning("Calculating P[X >= x]")
  }
  overlap.pvalue <- phyper(q = actual.overlap - adjust.value, 
                           m = length(list1), n = total.size - length(list1), k = length(list2), 
                           lower.tail = lower.tail)
  # return(c(actual.overlap, expected.overlap, overlap.pvalue))
  return(c(overlap.pvalue))
}

calculate.phyper.p.value(list1 = enriched.motifs_new$motif[enriched.motifs_new$cluster=="Rod"],
                         list2 = enriched.motifs_new$motif[enriched.motifs_new$cluster=="RodBC"],
                         total.size = length(unique(enriched.motifs$motif)),
                         lower.tail = F)

motif_jaccard_index_mtx <- data.frame(matrix(data = 1,nrow = length(levels(Idents(DMmouse_CRE))),ncol = length(levels(Idents(DMmouse_CRE)))))
rownames(motif_jaccard_index_mtx) <- levels(Idents(DMmouse_CRE))
colnames(motif_jaccard_index_mtx) <- levels(Idents(DMmouse_CRE))


for (i in 1:nrow(motif_jaccard_index_mtx)) {
  for (j in 1:ncol(motif_jaccard_index_mtx)) {
    temp1 <- enriched.motifs_new$motif[enriched.motifs_new$cluster==rownames(motif_jaccard_index_mtx)[i]]
    temp2 <- enriched.motifs_new$motif[enriched.motifs_new$cluster==rownames(motif_jaccard_index_mtx)[j]]
    
    motif_jaccard_index_mtx[i,j] <- length(intersect(temp1,temp2))/length(union(temp1,temp2))
    if (i==j) {
      motif_jaccard_index_mtx[i,j] <- NA
    }
  }
}



pheatmap(motif_jaccard_index_mtx,
         # scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         # col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         col = SpatialColors(n = 100)[50:1],
         # show_rownames = F,
         # show_colnames = F,
         treeheight_row = F,
         treeheight_col = F,
         border = F,
         legend = T)




motif_p_value_mtx <- data.frame(matrix(data = 1,nrow = length(levels(Idents(DMmouse_CRE))),ncol = length(levels(Idents(DMmouse_CRE)))))
rownames(motif_p_value_mtx) <- levels(Idents(DMmouse_CRE))
colnames(motif_p_value_mtx) <- levels(Idents(DMmouse_CRE))


for (i in 1:nrow(motif_p_value_mtx)) {
  for (j in 1:ncol(motif_p_value_mtx)) {
    temp1 <- enriched.motifs_new$motif[enriched.motifs_new$cluster==rownames(motif_p_value_mtx)[i]]
    temp2 <- enriched.motifs_new$motif[enriched.motifs_new$cluster==rownames(motif_p_value_mtx)[j]]
    
    motif_p_value_mtx[i,j] <- calculate.phyper.p.value(list1 = temp1,list2 = temp2,total.size = length(unique(enriched.motifs$motif)),lower.tail = F)
    if (i==j) {
      motif_p_value_mtx[i,j] <- NA
    }
  }
}


pvalue <- c()
for (i in 1:(nrow(motif_p_value_mtx)-1)) {
  pvalue <- c(pvalue,unname(motif_p_value_mtx[i,(i+1):nrow(motif_p_value_mtx)]))
}
pvalue <- data.frame(p.value=unlist(pvalue),adjust=p.adjust(unlist(pvalue),method = "BH"))

for (i in 1:nrow(motif_p_value_mtx)) {
  for (j in 1:ncol(motif_p_value_mtx)) {
    if (motif_p_value_mtx[i,j] %in% pvalue$p.value) {
      motif_p_value_mtx[i,j] <- pvalue$adjust[pvalue$p.value==motif_p_value_mtx[i,j]][1]
    }
  }
}


motif_p_value_mtx <- (-log10(motif_p_value_mtx))

pheatmap(motif_p_value_mtx,
         # scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         # col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         col = SpatialColors(n = 100)[51:100],
         # show_rownames = F,
         # show_colnames = F,
         treeheight_row = F,
         treeheight_col = F,
         border = F,
         legend = T)






mk_jaccard_index_mtx <- data.frame(matrix(data = 1,nrow = length(levels(Idents(DMmouse_CRE))),ncol = length(levels(Idents(DMmouse_CRE)))))
rownames(mk_jaccard_index_mtx) <- levels(Idents(DMmouse_CRE))
colnames(mk_jaccard_index_mtx) <- levels(Idents(DMmouse_CRE))


for (i in 1:nrow(mk_jaccard_index_mtx)) {
  for (j in 1:ncol(mk_jaccard_index_mtx)) {
    temp1 <- mk.CTRL$gene[mk.CTRL$cluster==rownames(mk_jaccard_index_mtx)[i]]
    temp2 <- mk.CTRL$gene[mk.CTRL$cluster==rownames(mk_jaccard_index_mtx)[j]]
    
    mk_jaccard_index_mtx[i,j] <- length(intersect(temp1,temp2))/length(union(temp1,temp2))
    if (i==j) {
      mk_jaccard_index_mtx[i,j] <- NA
    }
  }
}



pheatmap(mk_jaccard_index_mtx,
         # scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         # col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         col = SpatialColors(n = 100)[50:1],
         # show_rownames = F,
         # show_colnames = F,
         treeheight_row = F,
         treeheight_col = F,
         border = F,
         legend = T)



mk_p_value_mtx <- data.frame(matrix(data = 1,nrow = length(levels(Idents(DMmouse_CRE))),ncol = length(levels(Idents(DMmouse_CRE)))))
rownames(mk_p_value_mtx) <- levels(Idents(DMmouse_CRE))
colnames(mk_p_value_mtx) <- levels(Idents(DMmouse_CRE))


for (i in 1:nrow(mk_p_value_mtx)) {
  for (j in 1:ncol(mk_p_value_mtx)) {
    temp1 <- mk.CTRL$gene[mk.CTRL$cluster==rownames(mk_p_value_mtx)[i]]
    temp2 <- mk.CTRL$gene[mk.CTRL$cluster==colnames(mk_p_value_mtx)[j]]

    mk_p_value_mtx[i,j] <- calculate.phyper.p.value(list1 = temp1,list2 = temp2,total.size = length(unique(CRE$V4)),lower.tail = F)
    if (i==j) {
      mk_p_value_mtx[i,j] <- NA
    }
  }
}


pvalue <- c()
for (i in 1:(nrow(mk_p_value_mtx)-1)) {
  pvalue <- c(pvalue,unname(mk_p_value_mtx[i,(i+1):nrow(mk_p_value_mtx)]))
}
pvalue <- data.frame(p.value=unlist(pvalue),adjust=p.adjust(unlist(pvalue),method = "BH"))

for (i in 1:nrow(mk_p_value_mtx)) {
  for (j in 1:ncol(mk_p_value_mtx)) {
    if (mk_p_value_mtx[i,j] %in% pvalue$p.value) {
      mk_p_value_mtx[i,j] <- pvalue$adjust[pvalue$p.value==mk_p_value_mtx[i,j]][1]
    }
  }
}

mk_p_value_mtx <- (-log10(mk_p_value_mtx))

mk_p_value_mtx[mk_p_value_mtx==Inf] <- max(unlist(as.vector(mk_p_value_mtx))[!unlist(as.vector(mk_p_value_mtx)) %in% c(NA,Inf)])

pheatmap(mk_p_value_mtx,
         # scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         # col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         # col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         col = SpatialColors(n = 100)[51:100],
         # show_rownames = F,
         # show_colnames = F,
         treeheight_row = F,
         treeheight_col = F,
         border = F,
         legend = T)




# TF_gene <- enriched.motifs_new$motif.name %>% gsub(pattern = "[(]var.1[)]",replacement = "") %>% gsub(pattern = "[(]var.2[)]",replacement = "") %>% gsub(pattern = "[(]var.3[)]",replacement = "") %>% strsplit(split="::") %>% unlist() %>% str_to_title() %>% unique()
# 
# mk.CTRL$gene[mk.CTRL$cluster=="Rod"] %>% strsplit(",") %>% unlist() %>% strsplit("@") %>% lapply(function(x) x[2]) %>% unlist() %>% grep(pattern = "ADDG",value = T,invert = T) %>% unique() %>% intersect(y = TF_gene)
# 
# 
# 
# enriched.motifs_new$motif.name[enriched.motifs_new$cluster=="Rod"] %>% gsub(pattern = "[(]var.1[)]",replacement = "") %>% gsub(pattern = "[(]var.2[)]",replacement = "") %>% gsub(pattern = "[(]var.3[)]",replacement = "") %>% strsplit(split="::") %>% unlist() %>% str_to_title() %>% unique()
# 
# 
# venn.diagram(x = list("DEG"=mk.CTRL$gene[mk.CTRL$cluster=="Micro"] %>% strsplit(",") %>% unlist() %>% strsplit("@") %>% lapply(function(x) x[2]) %>% unlist() %>% grep(pattern = "ADDG",value = T,invert = T) %>% unique() %>% intersect(y = TF_gene),
#                       "Motif"=enriched.motifs_new$motif.name[enriched.motifs_new$cluster=="Micro"] %>% gsub(pattern = "[(]var.1[)]",replacement = "") %>% gsub(pattern = "[(]var.2[)]",replacement = "") %>% gsub(pattern = "[(]var.3[)]",replacement = "") %>% strsplit(split="::") %>% unlist() %>% str_to_title() %>% unique()
#                       ),
#              imagetype = "png",
#              hyper.test = T,
#              total.population = length(TF_gene),
#              lower.tail = F,
#              lwd = 4,
#              fill = c("cornflowerblue", "darkorchid1"),
#              alpha = 0.75,
#              label.col = "black",
#              cex = 2,
#              fontfamily = "serif",
#              fontface = "bold",
#              cat.col = c("cornflowerblue", "darkorchid1"),
#              cat.cex = 2,
#              cat.fontfamily = "serif",
#              cat.fontface = "bold",
#              cat.dist = c(0.03, 0.03),
#              cat.pos = c(-20, 14),
#              filename = "temp.png")




# DMmouse_ATAC.CTRL <- Footprint(object = DMmouse_ATAC.CTRL,motif.name ="PITX3",genome = BSgenome.Mmusculus.UCSC.mm10,in.peaks = TRUE,compute.expected = T)
# 
# Signac::Footprint(object = DMmouse_ATAC.CTRL,motif.name ="MA0467.1",genome = BSgenome.Mmusculus.UCSC.mm10)
# 
# PlotFootprint(DMmouse_ATAC.CTRL, features = c("Crx"))


# saveRDS(DMmouse_ATAC,"DMmouse_ATAC.RDS")
# saveRDS(DMmouse_ATAC.CTRL,"DMmouse_ATAC.CTRL.RDS")

# DMmouse_ATAC <- readRDS("DMmouse_ATAC.RDS")
# DMmouse_ATAC.CTRL <- readRDS("DMmouse_ATAC.CTRL.RDS")





## 全细胞建立ccans
DMmouse_ATAC.cds <- as.cell_data_set(DMmouse_ATAC)

DMmouse_ATAC.cicero <- make_cicero_cds(DMmouse_ATAC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(DMmouse_ATAC)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
# genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(DMmouse_ATAC.cicero, genomic_coords = genome.df, sample_num = 100)

ccans <- generate_ccans(conns,coaccess_cutoff_override = 0.02)


ccans_new <- merge(ccans,data.frame(Peak=paste0(CRE$V1,"-",CRE$V2,"-",CRE$V3),V2=CRE$V4, Start=CRE$V2,End=CRE$V3),by="Peak")
ccans_new <- ccans_new[order(ccans_new$Start),]
ccans_new <- ccans_new[order(ccans_new$CCAN),]

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(DMmouse_ATAC) <- links




## CTRL建立ccans
DMmouse_ATAC.CTRL.cds <- as.cell_data_set(DMmouse_ATAC.CTRL)

DMmouse_ATAC.CTRL.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(DMmouse_ATAC.CTRL)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
# genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns_CTRL <- run_cicero(DMmouse_ATAC.CTRL.cicero, genomic_coords = genome.df, sample_num = 100)

ccans_CTRL <- generate_ccans(conns_CTRL,coaccess_cutoff_override = 0.02)


ccans_CTRL_new <- merge(ccans_CTRL,data.frame(Peak=paste0(CRE$V1,"-",CRE$V2,"-",CRE$V3),V2=CRE$V4, Start=CRE$V2,End=CRE$V3),by="Peak")
ccans_CTRL_new <- ccans_CTRL_new[order(ccans_CTRL_new$Start),]
ccans_CTRL_new <- ccans_CTRL_new[order(ccans_CTRL_new$CCAN),]

# links <- ConnectionsToLinks(conns = conns_CTRL, ccans = ccans_CTRL)
links <- ConnectionsToLinks(conns = conns_CTRL, ccans = ccans)
Links(DMmouse_ATAC.CTRL) <- links

# CoveragePlot(DMmouse_ATAC.CTRL, region = "chr1-10008836-10187421",)





## DM建立ccans

DMmouse_ATAC.DM <- subset(DMmouse_ATAC, subset = var == "DM")

DMmouse_ATAC.DM.cds <- as.cell_data_set(DMmouse_ATAC.DM)

DMmouse_ATAC.DM.cicero <- make_cicero_cds(DMmouse_ATAC.DM.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(DMmouse_ATAC.DM)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
# genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns_DM <- run_cicero(DMmouse_ATAC.DM.cicero, genomic_coords = genome.df, sample_num = 100)

ccans_DM <- generate_ccans(conns_DM,coaccess_cutoff_override = 0.02)


ccans_DM_new <- merge(ccans_DM,data.frame(Peak=paste0(CRE$V1,"-",CRE$V2,"-",CRE$V3),V2=CRE$V4, Start=CRE$V2,End=CRE$V3),by="Peak")
ccans_DM_new <- ccans_DM_new[order(ccans_DM_new$Start),]
ccans_DM_new <- ccans_DM_new[order(ccans_DM_new$CCAN),]

# links <- ConnectionsToLinks(conns = conns_DM, ccans = ccans_DM)
links <- ConnectionsToLinks(conns = conns_DM, ccans = ccans)

Links(DMmouse_ATAC.DM) <- links

CoveragePlot(DMmouse_ATAC.CTRL, region = "chr1-10008836-10187421") +
  CoveragePlot(DMmouse_ATAC.DM, region = "chr1-10008836-10187421")


# FragmentHistogram()
num <- 1
df_new <- data.frame()
for (cluster in c("Rod","Cone","RodBC","ConeBC","MG","AC","Micro")) {
  df <- data.frame()
  a <- data.frame(Links(get(paste0("DMmouse_ATAC.CTRL.",cluster))))
  b <- data.frame(Links(get(paste0("DMmouse_ATAC.DM.",cluster))))
  print(length(intersect(unique(a$group),unique(b$group)))/length(union(unique(a$group),unique(b$group))))
  for (i in unique(ccans$CCAN)) {
    if (i %in% unique(a$group)) {
      temp_1 <- mean(a[a$group %in% i,]$score)
    } else {
      temp_1 <- 0
    }
    
    if (i %in% unique(b$group)) {
      temp_2 <- mean(b[b$group %in% i,]$score)
    } else {
      temp_2 <- 0
    }
    
    df <- rbind(df,data.frame(group=i,CTRL=temp_1,DM=temp_2))
  }
  rownames(df) <- df$group
  df <- df[,-1]
  
  
  my_comparison <- list(c("CTRL","DM"))
  p <- ggpaired(df, cond1 = "CTRL", cond2 = "DM",
           fill = "condition", palette = "jco",line.size = 1,point.size = 0)+
    stat_compare_means(comparisons = my_comparison,aes(label=p.sign..),label = "p-value", method = "t.test")
  assign(paste0("p",num),p) 
  df_new <- rbind(df_new,
                  data.frame(melt(df),cluster=cluster))
  num <- num + 1
}


# df_new <- melt(df)
# 
# ggdensity(df_new[df_new$value<0.2,], x = "value",
#           add = "median", rug = TRUE,
#           color = "variable", palette = c("#00AFBB", "#E7B800"))
# 
# 
# 
# ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6,p7),ncol = 7)



















ggboxplot(df_new, x = "cluster", y = "value", fill = "variable" ) +
  # stat_compare_means(aes(group=variable), label = "p.signif",method = "t.test")+
  labs(y="Score",x="") +
  scale_fill_manual(values = c("#00468B","#ED0000")) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=45,hjust = 1,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  guides(fill=FALSE)




## CTRL建立ccans Rod
DMmouse_ATAC.CTRL.Rod <- subset(DMmouse_ATAC.CTRL, idents = "Rod")
DMmouse_ATAC.CTRL.Rod.cds <- as.cell_data_set(DMmouse_ATAC.CTRL.Rod)
DMmouse_ATAC.CTRL.Rod.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.Rod.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.Rod.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.CTRL.Rod)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_CTRL.Rod <- run_cicero(DMmouse_ATAC.CTRL.Rod.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_CTRL.Rod <- generate_ccans(conns_CTRL.Rod,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_CTRL.Rod, ccans = ccans)
Links(DMmouse_ATAC.CTRL.Rod) <- links
print("finish")


## CTRL建立ccans RodBC
DMmouse_ATAC.CTRL.RodBC <- subset(DMmouse_ATAC.CTRL, idents = "RodBC")
DMmouse_ATAC.CTRL.RodBC.cds <- as.cell_data_set(DMmouse_ATAC.CTRL.RodBC)
DMmouse_ATAC.CTRL.RodBC.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.RodBC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.RodBC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.CTRL.RodBC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_CTRL.RodBC <- run_cicero(DMmouse_ATAC.CTRL.RodBC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_CTRL.RodBC <- generate_ccans(conns_CTRL.RodBC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_CTRL.RodBC, ccans = ccans)
Links(DMmouse_ATAC.CTRL.RodBC) <- links
print("finish")


## CTRL建立ccans RodBC
DMmouse_ATAC.CTRL.RodBC <- subset(DMmouse_ATAC.CTRL, idents = "Rod BC")
DMmouse_ATAC.CTRL.RodBC.cds <- as.cell_data_set(DMmouse_ATAC.CTRL.RodBC)
DMmouse_ATAC.CTRL.RodBC.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.RodBC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.RodBC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.CTRL.RodBC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_CTRL.RodBC <- run_cicero(DMmouse_ATAC.CTRL.RodBC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_CTRL.RodBC <- generate_ccans(conns_CTRL.RodBC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_CTRL.RodBC, ccans = ccans)
Links(DMmouse_ATAC.CTRL.RodBC) <- links
print("finish")


## CTRL建立ccans ConeBC
DMmouse_ATAC.CTRL.ConeBC <- subset(DMmouse_ATAC.CTRL, idents = "Cone BC")
DMmouse_ATAC.CTRL.ConeBC.cds <- as.cell_data_set(DMmouse_ATAC.CTRL.ConeBC)
DMmouse_ATAC.CTRL.ConeBC.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.ConeBC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.ConeBC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.CTRL.ConeBC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_CTRL.ConeBC <- run_cicero(DMmouse_ATAC.CTRL.ConeBC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_CTRL.ConeBC <- generate_ccans(conns_CTRL.ConeBC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_CTRL.ConeBC, ccans = ccans)
Links(DMmouse_ATAC.CTRL.ConeBC) <- links
print("finish")


## CTRL建立ccans MG
DMmouse_ATAC.CTRL.MG <- subset(DMmouse_ATAC.CTRL, idents = "MG")
DMmouse_ATAC.CTRL.MG.cds <- as.cell_data_set(DMmouse_ATAC.CTRL.MG)
DMmouse_ATAC.CTRL.MG.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.MG.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.MG.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.CTRL.MG)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_CTRL.MG <- run_cicero(DMmouse_ATAC.CTRL.MG.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_CTRL.MG <- generate_ccans(conns_CTRL.MG,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_CTRL.MG, ccans = ccans)
Links(DMmouse_ATAC.CTRL.MG) <- links
print("finish")


## CTRL建立ccans AC
DMmouse_ATAC.CTRL.AC <- subset(DMmouse_ATAC.CTRL, idents = "AC")
DMmouse_ATAC.CTRL.AC.cds <- as.cell_data_set(DMmouse_ATAC.CTRL.AC)
DMmouse_ATAC.CTRL.AC.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.AC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.AC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.CTRL.AC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_CTRL.AC <- run_cicero(DMmouse_ATAC.CTRL.AC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_CTRL.AC <- generate_ccans(conns_CTRL.AC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_CTRL.AC, ccans = ccans)
Links(DMmouse_ATAC.CTRL.AC) <- links
print("finish")


## CTRL建立ccans Micro
DMmouse_ATAC.CTRL.Micro <- subset(DMmouse_ATAC.CTRL, idents = "Micro")
DMmouse_ATAC.CTRL.Micro.cds <- as.cell_data_set(DMmouse_ATAC.CTRL.Micro)
DMmouse_ATAC.CTRL.Micro.cicero <- make_cicero_cds(DMmouse_ATAC.CTRL.Micro.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.CTRL.Micro.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.CTRL.Micro)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_CTRL.Micro <- run_cicero(DMmouse_ATAC.CTRL.Micro.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_CTRL.Micro <- generate_ccans(conns_CTRL.Micro,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_CTRL.Micro, ccans = ccans)
Links(DMmouse_ATAC.CTRL.Micro) <- links
print("finish")










## DM建立ccans Rod
DMmouse_ATAC.DM.Rod <- subset(DMmouse_ATAC.DM, idents = "Rod")
DMmouse_ATAC.DM.Rod.cds <- as.cell_data_set(DMmouse_ATAC.DM.Rod)
DMmouse_ATAC.DM.Rod.cicero <- make_cicero_cds(DMmouse_ATAC.DM.Rod.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.Rod.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.DM.Rod)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_DM.Rod <- run_cicero(DMmouse_ATAC.DM.Rod.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_DM.Rod <- generate_ccans(conns_DM.Rod,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_DM.Rod, ccans = ccans)
Links(DMmouse_ATAC.DM.Rod) <- links
print("finish")


## DM建立ccans RodBC
DMmouse_ATAC.DM.RodBC <- subset(DMmouse_ATAC.DM, idents = "RodBC")
DMmouse_ATAC.DM.RodBC.cds <- as.cell_data_set(DMmouse_ATAC.DM.RodBC)
DMmouse_ATAC.DM.RodBC.cicero <- make_cicero_cds(DMmouse_ATAC.DM.RodBC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.RodBC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.DM.RodBC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_DM.RodBC <- run_cicero(DMmouse_ATAC.DM.RodBC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_DM.RodBC <- generate_ccans(conns_DM.RodBC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_DM.RodBC, ccans = ccans)
Links(DMmouse_ATAC.DM.RodBC) <- links
print("finish")


## DM建立ccans RodBC
DMmouse_ATAC.DM.RodBC <- subset(DMmouse_ATAC.DM, idents = "Rod BC")
DMmouse_ATAC.DM.RodBC.cds <- as.cell_data_set(DMmouse_ATAC.DM.RodBC)
DMmouse_ATAC.DM.RodBC.cicero <- make_cicero_cds(DMmouse_ATAC.DM.RodBC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.RodBC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.DM.RodBC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_DM.RodBC <- run_cicero(DMmouse_ATAC.DM.RodBC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_DM.RodBC <- generate_ccans(conns_DM.RodBC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_DM.RodBC, ccans = ccans)
Links(DMmouse_ATAC.DM.RodBC) <- links
print("finish")


## DM建立ccans ConeBC
DMmouse_ATAC.DM.ConeBC <- subset(DMmouse_ATAC.DM, idents = "Cone BC")
DMmouse_ATAC.DM.ConeBC.cds <- as.cell_data_set(DMmouse_ATAC.DM.ConeBC)
DMmouse_ATAC.DM.ConeBC.cicero <- make_cicero_cds(DMmouse_ATAC.DM.ConeBC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.ConeBC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.DM.ConeBC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_DM.ConeBC <- run_cicero(DMmouse_ATAC.DM.ConeBC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_DM.ConeBC <- generate_ccans(conns_DM.ConeBC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_DM.ConeBC, ccans = ccans)
Links(DMmouse_ATAC.DM.ConeBC) <- links
print("finish")


## DM建立ccans MG
DMmouse_ATAC.DM.MG <- subset(DMmouse_ATAC.DM, idents = "MG")
DMmouse_ATAC.DM.MG.cds <- as.cell_data_set(DMmouse_ATAC.DM.MG)
DMmouse_ATAC.DM.MG.cicero <- make_cicero_cds(DMmouse_ATAC.DM.MG.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.MG.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.DM.MG)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_DM.MG <- run_cicero(DMmouse_ATAC.DM.MG.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_DM.MG <- generate_ccans(conns_DM.MG,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_DM.MG, ccans = ccans)
Links(DMmouse_ATAC.DM.MG) <- links
print("finish")


## DM建立ccans AC
DMmouse_ATAC.DM.AC <- subset(DMmouse_ATAC.DM, idents = "AC")
DMmouse_ATAC.DM.AC.cds <- as.cell_data_set(DMmouse_ATAC.DM.AC)
DMmouse_ATAC.DM.AC.cicero <- make_cicero_cds(DMmouse_ATAC.DM.AC.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.AC.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.DM.AC)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_DM.AC <- run_cicero(DMmouse_ATAC.DM.AC.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_DM.AC <- generate_ccans(conns_DM.AC,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_DM.AC, ccans = ccans)
Links(DMmouse_ATAC.DM.AC) <- links
print("finish")


## DM建立ccans Micro
DMmouse_ATAC.DM.Micro <- subset(DMmouse_ATAC.DM, idents = "Micro")
DMmouse_ATAC.DM.Micro.cds <- as.cell_data_set(DMmouse_ATAC.DM.Micro)
DMmouse_ATAC.DM.Micro.cicero <- make_cicero_cds(DMmouse_ATAC.DM.Micro.cds, reduced_coordinates = reducedDims(DMmouse_ATAC.DM.Micro.cds)$UMAP)
genome <- seqlengths(DMmouse_ATAC.DM.Micro)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns_DM.Micro <- run_cicero(DMmouse_ATAC.DM.Micro.cicero, genomic_coords = genome.df, sample_num = 100)
ccans_DM.Micro <- generate_ccans(conns_DM.Micro,coaccess_cutoff_override = 0.02)
links <- ConnectionsToLinks(conns = conns_DM.Micro, ccans = ccans)
Links(DMmouse_ATAC.DM.Micro) <- links
print("finish")

