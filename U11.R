# generate PFMatrix of RBP

RBP_Information <- read.table("~/zjw/20220307_DMmouse/CISBP_RNA/mmu/RBP_Information_all_motifs.txt",sep = "\t",header = 1)

for (i in list.files("~/zjw/20220307_DMmouse/CISBP_RNA/mmu/pwms_all_motifs/")) {

  if (!file.info(paste0("CISBP_RNA/mmu/pwms_all_motifs/",i))$size==0) {
    pwms <- read.table(paste0("CISBP_RNA/mmu/pwms_all_motifs/",i),header = 1)[,-1]
    colnames(pwms)[4] <- "T"
    pwms <- pwms %>% as.matrix() %>% t()
    pwms <- PWMatrix(ID = gsub(pattern = ".txt",replacement = "",x = i),
                     name = gsub(pattern = ".txt",replacement = "",x = i),
                     matrixClass = "Zipper-Type",strand = "+",
                     bg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                     tags = list(species="mmu",tax_group="vertebrates"),
                     profileMatrix = pwms)
    assign(gsub(pattern = ".txt",replacement = "",x = i),pwms)
    
    if (i=="M001_0.6.txt") {
      char <- paste0(gsub(pattern = ".txt",replacement = "",x = i),"=",gsub(pattern = ".txt",replacement = "",x = i))
    } else {
      char <- paste(char,paste0(gsub(pattern = ".txt",replacement = "",x = i)," = ",gsub(pattern = ".txt",replacement = "",x = i)),sep = " ,")
    }
  }
  
}

char <- paste0('pwm <- PWMatrixList(',char,')')

eval(parse(text = char))

names

counts_new <- DMmouse_CRE@assays$RNA@counts[1:317,]
rownames(counts_new) <- rownames(isoform_switch_bed)



chrom_assay_new <- CreateChromatinAssay(
  counts = counts_new,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'snapATAC/output_sorted_new.bed.gz',
  min.cells = 0,
  min.features = 0
)

DMmouse_ATAC_new <- CreateSeuratObject(
  counts = chrom_assay_new,
  assay = "peaks",
  # meta.data = metadata
)

Idents(DMmouse_ATAC_new) <- Idents(DMmouse_CRE)

DMmouse_ATAC_new$celltype <- DMmouse_CRE$celltype
DMmouse_ATAC_new$var <- DMmouse_CRE$var

# add motif information
DMmouse_ATAC_new <- AddMotifs(object = DMmouse_ATAC_new, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pwm)

open.peaks <- AccessiblePeaks(DMmouse_ATAC_new, idents = levels(mk.CTRL$cluster))
meta.feature <- GetAssayData(DMmouse_ATAC_new, assay = "peaks", slot = "meta.features")

enriched.motifs_RBP_Longer <- data.frame()
for (i in c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")) {
  gene <- get(paste0("t.test_length_CTRL_",i))$gene[get(paste0("t.test_length_CTRL_",i))$threshold=="Longer"]
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[isoform_switch_bed$V4[isoform_switch_bed$V5 %in% gene], ],
    n = 50000
  )
  
  enriched.motifs_RBP_Longer <- rbind(enriched.motifs_RBP_Longer,
                                      data.frame(FindMotifs(object = DMmouse_ATAC_new,
                                                            background = peaks.matched,
                                                            features = isoform_switch_bed$V4[isoform_switch_bed$V5 %in% gene]),
                                                 cluster=i))
  
}

enriched.motifs_RBP_Shorter <- data.frame()
for (i in c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")) {
  gene <- get(paste0("t.test_length_CTRL_",i))$gene[get(paste0("t.test_length_CTRL_",i))$threshold=="Shorter"]
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[isoform_switch_bed$V4[isoform_switch_bed$V5 %in% gene], ],
    n = 50000
  )
  
  enriched.motifs_RBP_Shorter <- rbind(enriched.motifs_RBP_Shorter,
                                       data.frame(FindMotifs(object = DMmouse_ATAC_new,
                                                             background = peaks.matched,
                                                             features = isoform_switch_bed$V4[isoform_switch_bed$V5 %in% gene]),
                                                  cluster=i))
  
}


hm_up <- enriched.motifs_RBP_Longer[enriched.motifs_RBP_Longer$pvalue<0.05,c("motif","cluster")]



hm_up_new <- data.frame(gene=unique(hm_up$motif))
hm_up_new[,2:8] <- 0
hm_up_new <- hm_up_new[,-1]
rownames(hm_up_new) <- unique(hm_up$motif)
colnames(hm_up_new) <- c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")

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
         color = colorRampPalette(c("grey", "red"))(5),cluster_cols = F,cluster_rows = F,border = F,show_rownames = F)

# M298_0.6 Rbfox1/2/3  universal
# M140_0.6 Enox1/2 universal
# M323_0.6 Nova1 Rod https://pubmed.ncbi.nlm.nih.gov/10719891/
# M231_0.6 Eif2s1 Rod
# M152_0.6 Fxr1/2 Cone
# M353_0.6 Srsf7 Cone

# M250_0.6 Ybx1/2/3 MG
# M085_0.6 Zcrb1 MG
# M102_0.6 Srsf1/9 HC
# M048_0.6 Rbm3/Cirbp HC
# M037_0.6 Mbnl1/2/3 Micro
# M040_0.6 Msi1/2 Micro

MotifPlot(
  object = DMmouse_ATAC_new,
  motifs = c("M298_0.6","M140_0.6","M323_0.6","M231_0.6","M152_0.6","M353_0.6","M250_0.6","M085_0.6","M102_0.6","M048_0.6","M037_0.6","M040_0.6")
)


temp <- enriched.motifs_RBP_Longer$motif.name[enriched.motifs_RBP_Longer$cluster=="Cone" & enriched.motifs_RBP_Longer$pvalue<0.05]

RBP_Information$RBP_Name[RBP_Information$Motif_ID %in% temp]

go_onto <- enrichGO(RBP_Information$RBP_Name[RBP_Information$Motif_ID %in% temp], OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)

a <- data.frame(go_onto)







all_term <- data.frame()
for (i in unique(enriched.motifs_RBP_Longer$cluster)) {
  temp <- enriched.motifs_RBP_Longer$motif.name[enriched.motifs_RBP_Longer$cluster==i & enriched.motifs_RBP_Longer$pvalue<0.05]
  temp <- RBP_Information$RBP_Name[RBP_Information$Motif_ID %in% temp]

  go_onto <- enrichGO(temp, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  #go_onto <- simplify(go_onto)
  a <- data.frame(go_onto)
  all_term <- rbind(all_term,
                    data.frame(a,cluster=i))
}




go_term <- rbind(all_term[all_term$Description %in% c("RNA splicing","mRNA processing") & all_term$cluster=="Rod",],
                 all_term[all_term$Description %in% c("regulation of RNA splicing","regulation of RNA binding") & all_term$cluster=="Cone",],
                 all_term[all_term$Description %in% c("regulation of translation","mRNA catabolic process") & all_term$cluster=="Rod_BC",],
                 all_term[all_term$Description %in% c("mRNA splicing, via spliceosome","regulation of cellular amide metabolic process") & all_term$cluster=="Cone_BC",],
                 all_term[all_term$Description %in% c("nucleic acid transport","mRNA destabilization") & all_term$cluster=="MG",],
                 all_term[all_term$Description %in% c("mRNA splice site selection","regulation of translation") & all_term$cluster=="AC",],
                 all_term[all_term$Description %in% c("RNA localization","mRNA modification") & all_term$cluster=="Micro",]
)

go_term$p.adjust <- (-log10(go_term$p.adjust))
# go_term$cluster <- factor(go_term$cluster,levels = c("LSC","LEC","ETAC","CTAC","CEC","ConjEC","MC","IC"))
go_term <- go_term[rev(seq(1:nrow(go_term))),]



#有一个item富集了两次
go_term$Description[duplicated(go_term$Description)] <- paste0(go_term$Description[duplicated(go_term$Description)],"--",go_term$cluster[duplicated(go_term$Description)])


ggbarplot(go_term,x = "Description",y="p.adjust",fill = "cluster",color = "cluster") +
  theme(axis.text.x = element_text(face = "plain",size = 12,angle = 270, hjust = 0,vjust = 0.5),
        axis.text.y = element_text(face = "plain",size = 12,angle = 270, hjust = 0.5,vjust = 0.5),
        axis.title.y = element_text(angle = 270),
        legend.position="none") +
  labs(x="",y="-log10 adjusted p-value")


# 
# Rod      RNA splicing                           mRNA processing
# Cone     regulation of RNA splicing             regulation of RNA binding
# Rod BC   regulation of translation              mRNA catabolic process
# Cone BC  mRNA splicing, via spliceosome         regulation of cellular amide metabolic process
# MG       nucleic acid transport                 mRNA destabilization
# AC       mRNA splice site selection             regulation of translation
# Micro    RNA localization                       mRNA modification



