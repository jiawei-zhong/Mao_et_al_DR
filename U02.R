# library(MASS)
options(bedtools.path = "~/anaconda3/bin/")
.libPaths("/disk2/user/jizhon/anaconda3/envs/R4.1/lib/R/library")

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
library(bedtoolsr)
library(systemPipeR)


gtf <- read.table("~/index/mm10/gencode.vM18.annotation.sorted.gtf",sep = "\t")
gtf <- gtf[gtf$V3=="start_codon",]

gene_id <- strsplit(gtf[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

transcript_id <- strsplit(gtf[,9],split = "transcript_id ")
transcript_id <- lapply(transcript_id, function(x) x[2])
transcript_id <- unlist(transcript_id)
transcript_id <- strsplit(transcript_id,split = ";")
transcript_id <- lapply(transcript_id, function(x) x[1])
transcript_id <- unlist(transcript_id)

gene_type <- strsplit(gtf[,9],split = "gene_type ")
gene_type <- lapply(gene_type, function(x) x[2])
gene_type <- unlist(gene_type)
gene_type <- strsplit(gene_type,split = ";")
gene_type <- lapply(gene_type, function(x) x[1])
gene_type <- unlist(gene_type)

gene_name <- strsplit(gtf[,9],split = "gene_name ")
gene_name <- lapply(gene_name, function(x) x[2])
gene_name <- unlist(gene_name)
gene_name <- strsplit(gene_name,split = ";")
gene_name <- lapply(gene_name, function(x) x[1])
gene_name <- unlist(gene_name)

transcript_name <- strsplit(gtf[,9],split = "transcript_name ")
transcript_name <- lapply(transcript_name, function(x) x[2])
transcript_name <- unlist(transcript_name)
transcript_name <- strsplit(transcript_name,split = ";")
transcript_name <- lapply(transcript_name, function(x) x[1])
transcript_name <- unlist(transcript_name)


gtf$gene_id <- gene_id
gtf$gene_name <- gene_name
gtf$gene_type <- gene_type
gtf$transcript_id <- transcript_id
gtf$transcript_name <- transcript_name
transcript_has_multi_start_codon <- (gtf$transcript_id %>% table())[(gtf$transcript_id %>% table())>1] %>% names()
transcript_has_multi_start_codon_corr_gene <- gtf$gene_id[gtf$transcript_id %in% transcript_has_multi_start_codon] %>% unique()
gtf <- gtf[!(gtf$gene_id %in% transcript_has_multi_start_codon_corr_gene),]


# generate transcript set
region <- read.table("~/zjw/annotation/gencode_mm10_all_gene_all_transcript_region.bed",sep = "\t")
region <- region[grep("UTR",region$V7),]
region$gene_id_new <- substr(region$V4,1,18)
substr(gene_id,1,18)

all_gene_all_transcript <- read.table("~/index/mm10/gencode.vM18.annotation.all.gene.all.transcript.gtf",sep = "\t")

gene_id <- strsplit(all_gene_all_transcript[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

transcript_id <- strsplit(all_gene_all_transcript[,9],split = "transcript_id ")
transcript_id <- lapply(transcript_id, function(x) x[2])
transcript_id <- unlist(transcript_id)
transcript_id <- strsplit(transcript_id,split = ";")
transcript_id <- lapply(transcript_id, function(x) x[1])
transcript_id <- unlist(transcript_id)

gene_type <- strsplit(all_gene_all_transcript[,9],split = "gene_type ")
gene_type <- lapply(gene_type, function(x) x[2])
gene_type <- unlist(gene_type)
gene_type <- strsplit(gene_type,split = ";")
gene_type <- lapply(gene_type, function(x) x[1])
gene_type <- unlist(gene_type)

gene_name <- strsplit(all_gene_all_transcript[,9],split = "gene_name ")
gene_name <- lapply(gene_name, function(x) x[2])
gene_name <- unlist(gene_name)
gene_name <- strsplit(gene_name,split = ";")
gene_name <- lapply(gene_name, function(x) x[1])
gene_name <- unlist(gene_name)

transcript_name <- strsplit(all_gene_all_transcript[,9],split = "transcript_name ")
transcript_name <- lapply(transcript_name, function(x) x[2])
transcript_name <- unlist(transcript_name)
transcript_name <- strsplit(transcript_name,split = ";")
transcript_name <- lapply(transcript_name, function(x) x[1])
transcript_name <- unlist(transcript_name)

all_gene_all_transcript <- data.frame(chr=all_gene_all_transcript$V1,
                                      stop=all_gene_all_transcript$V4-1,
                                      end=all_gene_all_transcript$V5,
                                      strand=all_gene_all_transcript$V7,
                                      gene_id,
                                      gene_name,
                                      gene_type,
                                      transcript_id,
                                      transcript_name,
                                      gene_has_start_codon=F,
                                      gene_has_5UTR=F,
                                      gene_has_unique_transcript=F,
                                      gene_has_unique_start_codon=F,
                                      all_transcript_has_same_5UTR=F,
                                      all_transcript_5UTR_no_splicing=F)

all_gene_all_transcript_temp <- data.frame(gene_id_new=substr(gene_id,1,18),
                                           gene_id) %>% unique()

region <- merge(region,all_gene_all_transcript_temp,by="gene_id_new")
region <- region[,c("V1","V2","V3","V4","V5","V6","V7","gene_id")]

gene_has_same_gene_name <- ((all_gene_all_transcript[,c("gene_id","gene_name")] %>% unique())$gene_name %>% table())[((all_gene_all_transcript[,c("gene_id","gene_name")] %>% unique())$gene_name %>% table())>1] %>% names()
gene_has_same_gene_name <- all_gene_all_transcript$gene_id[all_gene_all_transcript$gene_name %in% gene_has_same_gene_name] %>% unique()

all_gene_all_transcript <- all_gene_all_transcript[!(all_gene_all_transcript$gene_id %in% c(gene_has_same_gene_name,transcript_has_multi_start_codon_corr_gene)),]
all_gene_all_transcript$gene_has_start_codon[all_gene_all_transcript$gene_id %in% unique(gtf$gene_id)] <- T
gtf <- gtf[!(gtf$gene_id %in% gene_has_same_gene_name),]

## 判断是否为unique transcript
all_gene_all_transcript$gene_has_5UTR[(all_gene_all_transcript$gene_id %>% strsplit(split = "[.]") %>% lapply(function(x) x[1]) %>% unlist()) %in% (region$V4[region$V7=="5UTR"] %>% strsplit(split = "__") %>% lapply(function(x) x[1]) %>% unlist() %>% unique())] <- T
all_gene_all_transcript$gene_has_unique_transcript[all_gene_all_transcript$gene_id %in% names(table(all_gene_all_transcript$gene_id)[table(all_gene_all_transcript$gene_id)==1])] <- T
# all_gene_all_transcript$gene_has_unique_start_codon[ all_gene_all_transcript$gene_has_unique_transcript==T] <- T
all_gene_all_transcript$all_transcript_has_same_5UTR[all_gene_all_transcript$gene_has_unique_transcript==T & all_gene_all_transcript$gene_has_5UTR==T] <- T


## 判断是否为unique start codon
start_codon <- rbind(data.frame(transcript_id = gtf$transcript_id[gtf$V7=="+"],start_codon_position = gtf$V4[gtf$V7=="+"]+1),
                     data.frame(transcript_id = gtf$transcript_id[gtf$V7=="-"],start_codon_position = gtf$V5[gtf$V7=="-"]))
all_gene_all_transcript <- merge(all_gene_all_transcript,start_codon,by='transcript_id',all=T)
all_gene_all_transcript <- all_gene_all_transcript[order(all_gene_all_transcript$gene_id),]
all_gene_all_transcript_temp <- na.omit(all_gene_all_transcript)
# all_gene_all_transcript_temp <- all_gene_all_transcript_temp[all_gene_all_transcript_temp$gene_has_5UTR==T & all_gene_all_transcript_temp$gene_has_unique_transcript==F,]
unique_start_codon_gene <- ((all_gene_all_transcript_temp[,c("gene_id","start_codon_position")] %>% unique())[,1] %>% table())[((all_gene_all_transcript_temp[,c("gene_id","start_codon_position")] %>% unique())[,1] %>% table())==1] %>% names()
all_gene_all_transcript$gene_has_unique_start_codon[all_gene_all_transcript$gene_id %in% unique_start_codon_gene] <- T
all_gene_all_transcript$all_transcript_has_same_5UTR[all_gene_all_transcript$gene_has_unique_start_codon==F] <- F


## 判断transcript的5UTR是否相同

gene_need_to_consider <- all_gene_all_transcript$gene_id[all_gene_all_transcript$gene_has_5UTR & all_gene_all_transcript$gene_has_unique_transcript==F & all_gene_all_transcript$gene_has_unique_start_codon] %>% unique()

gtf <- read.table("~/index/mm10/gencode.vM18.annotation.sorted.gtf",sep = "\t")
gtf <- gtf[gtf$V3 %in% c("UTR","start_codon"),]

gene_id <- strsplit(gtf[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

transcript_id <- strsplit(gtf[,9],split = "transcript_id ")
transcript_id <- lapply(transcript_id, function(x) x[2])
transcript_id <- unlist(transcript_id)
transcript_id <- strsplit(transcript_id,split = ";")
transcript_id <- lapply(transcript_id, function(x) x[1])
transcript_id <- unlist(transcript_id)

gene_type <- strsplit(gtf[,9],split = "gene_type ")
gene_type <- lapply(gene_type, function(x) x[2])
gene_type <- unlist(gene_type)
gene_type <- strsplit(gene_type,split = ";")
gene_type <- lapply(gene_type, function(x) x[1])
gene_type <- unlist(gene_type)

gene_name <- strsplit(gtf[,9],split = "gene_name ")
gene_name <- lapply(gene_name, function(x) x[2])
gene_name <- unlist(gene_name)
gene_name <- strsplit(gene_name,split = ";")
gene_name <- lapply(gene_name, function(x) x[1])
gene_name <- unlist(gene_name)

transcript_name <- strsplit(gtf[,9],split = "transcript_name ")
transcript_name <- lapply(transcript_name, function(x) x[2])
transcript_name <- unlist(transcript_name)
transcript_name <- strsplit(transcript_name,split = ";")
transcript_name <- lapply(transcript_name, function(x) x[1])
transcript_name <- unlist(transcript_name)

gtf$gene_id <- gene_id
gtf$gene_name <- gene_name
gtf$gene_type <- gene_type
gtf$transcript_id <- transcript_id
gtf$transcript_name <- transcript_name
gtf <- gtf[!(gtf$gene_id %in% c(transcript_has_multi_start_codon_corr_gene,gene_has_same_gene_name)),]
gtf <- gtf[gtf$gene_id %in% gene_need_to_consider,]

same_5UTR_df <- data.frame(gene_id=gtf$gene_id %>% unique(),
                           all_transcript_has_same_5UTR=F,
                           exon_number=NA)

for (i in 1:nrow(same_5UTR_df)) {
  # print(i)
  temp <- gtf[gtf$gene_id==same_5UTR_df$gene_id[i],]
  if (temp$V7[1]=="+") {
    temp <- temp[temp$V5 < (temp$V5[temp$V3=="start_codon"][1]),]
  } else {
    temp <- temp[temp$V5 > (temp$V5[temp$V3=="start_codon"][1]),]
  }
  if (length(unique(table(temp$transcript_id)))==1) {
    exon_number <- unique(table(temp$transcript_id))
    same_5UTR_df$exon_number[i] <- exon_number
    if ((temp[,c("V4","V5")] %>% unique() %>% nrow())==exon_number) {
      same_5UTR_df$all_transcript_has_same_5UTR[i] <- T
    }
  }
}

all_gene_all_transcript$all_transcript_has_same_5UTR[all_gene_all_transcript$gene_id %in% same_5UTR_df$gene_id[same_5UTR_df$all_transcript_has_same_5UTR]] <- T

same_5UTR_df_temp <- na.omit(same_5UTR_df)

all_gene_all_transcript$all_transcript_5UTR_no_splicing[all_gene_all_transcript$gene_id %in% same_5UTR_df_temp$gene_id[same_5UTR_df_temp$exon_number==1]] <- T

# all_gene_all_transcript$gene_has_start_codon==T &
#   all_gene_all_transcript$gene_has_5UTR==T &
#   all_gene_all_transcript$gene_has_unique_transcript==T  直接使用 直接使用region文件
#
# all_gene_all_transcript$gene_has_start_codon==T &
#   all_gene_all_transcript$gene_has_5UTR==T &
#   all_gene_all_transcript$gene_has_unique_transcript==F &
#   all_gene_all_transcript$gene_has_unique_start_codon==T &
#   all_gene_all_transcript$all_transcript_has_same_5UTR==T 直接使用 直接使用region文件
#
# all_gene_all_transcript$gene_has_start_codon==T &
#   all_gene_all_transcript$gene_has_5UTR==T &
#   all_gene_all_transcript$gene_has_unique_transcript==F &
#   all_gene_all_transcript$gene_has_unique_start_codon==T &
#   all_gene_all_transcript$all_transcript_has_same_5UTR==F &
#   all_gene_all_transcript$all_transcript_5UTR_no_splicing==T 直接使用 直接使用region文件
#
# all_gene_all_transcript$gene_has_start_codon==T &
#   all_gene_all_transcript$gene_has_5UTR==T &
#   all_gene_all_transcript$gene_has_unique_transcript==F &
#   all_gene_all_transcript$gene_has_unique_start_codon==F 将各种5UTR和CRE bedtools 选取intersect最多的，如果打平 选无splicing的和长的
#
# all_gene_all_transcript$gene_has_start_codon==T &
#   all_gene_all_transcript$gene_has_5UTR==T &
#   all_gene_all_transcript$gene_has_unique_transcript==F &
#   all_gene_all_transcript$gene_has_unique_start_codon==T &
#   all_gene_all_transcript$all_transcript_has_same_5UTR==F &
#   all_gene_all_transcript$all_transcript_5UTR_no_splicing==F 将各种5UTR和CRE bedtools 选取intersect最多的，如果打平 选无splicing的和长的

gene_UTR_directly_use <- c(all_gene_all_transcript$gene_id[all_gene_all_transcript$gene_has_start_codon==T & all_gene_all_transcript$gene_has_5UTR==T & all_gene_all_transcript$gene_has_unique_transcript==T] %>% unique(),
                           all_gene_all_transcript$gene_id[all_gene_all_transcript$gene_has_start_codon==T & all_gene_all_transcript$gene_has_5UTR==T & all_gene_all_transcript$gene_has_unique_transcript==F & all_gene_all_transcript$gene_has_unique_start_codon==T & all_gene_all_transcript$all_transcript_has_same_5UTR==T] %>% unique(),
                           all_gene_all_transcript$gene_id[all_gene_all_transcript$gene_has_start_codon==T & all_gene_all_transcript$gene_has_5UTR==T & all_gene_all_transcript$gene_has_unique_transcript==F & all_gene_all_transcript$gene_has_unique_start_codon==T & all_gene_all_transcript$all_transcript_has_same_5UTR==F & all_gene_all_transcript$all_transcript_5UTR_no_splicing==T] %>% unique())

gene_need_bedtools <-    c(all_gene_all_transcript$gene_id[all_gene_all_transcript$gene_has_start_codon==T & all_gene_all_transcript$gene_has_5UTR==T & all_gene_all_transcript$gene_has_unique_transcript==F & all_gene_all_transcript$gene_has_unique_start_codon==F] %>% unique(),
                           all_gene_all_transcript$gene_id[all_gene_all_transcript$gene_has_start_codon==T & all_gene_all_transcript$gene_has_5UTR==T & all_gene_all_transcript$gene_has_unique_transcript==F & all_gene_all_transcript$gene_has_unique_start_codon==T & all_gene_all_transcript$all_transcript_has_same_5UTR==F & all_gene_all_transcript$all_transcript_5UTR_no_splicing==F] %>% unique())

# 建立transcript UTR的bed文件   gene_need_bedtools  一个一个转录本判断
# fiveUTR_bed <- region[region$V7=="5UTR",]
# threeUTR_bed <- region[region$V7=="3UTR",]
#
#
# fiveUTR_inter_threeUTR <- bedtoolsr::bt.intersect(a = fiveUTR_bed, b = threeUTR_bed, wa = T, wb = T, s= T)
# # fiveUTR_inter_threeUTR <- fiveUTR_inter_threeUTR[fiveUTR_inter_threeUTR$V4==fiveUTR_inter_threeUTR$V11,]
#
# gene_need_bedtools_overlap <- intersect(gene_need_bedtools,unique(fiveUTR_inter_threeUTR$V8))
# gene_need_bedtools_no_overlap <- setdiff(gene_need_bedtools,unique(fiveUTR_inter_threeUTR$V8))

gtf <- read.table("~/index/mm10/gencode.vM18.annotation.sorted.gtf",sep = "\t")
gtf <- gtf[gtf$V3 %in% c("UTR","start_codon"),]

gene_id <- strsplit(gtf[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

transcript_id <- strsplit(gtf[,9],split = "transcript_id ")
transcript_id <- lapply(transcript_id, function(x) x[2])
transcript_id <- unlist(transcript_id)
transcript_id <- strsplit(transcript_id,split = ";")
transcript_id <- lapply(transcript_id, function(x) x[1])
transcript_id <- unlist(transcript_id)

gtf$gene_id <- gene_id
gtf$transcript_id <- transcript_id
gtf <- gtf[!(gtf$gene_id %in% c(transcript_has_multi_start_codon_corr_gene,gene_has_same_gene_name)),]
gtf <- gtf[gtf$transcript_id %in% intersect(gtf$transcript_id[gtf$V3=="UTR"],gtf$transcript_id[gtf$V3=="start_codon"]),]
gtf_bed <- data.frame(chr=gtf$V1,
                      start=gtf$V4-1,
                      end=gtf$V5,
                      gene_id=gtf$gene_id,
                      transcript_id=gtf$transcript_id,
                      strand=gtf$V7,
                      type=gtf$V3)

# gtf_UTR <- gtf[gtf$V3=="UTR",]
# gtf_UTR_bed <- gtf_bed[gtf_bed$type=="UTR",]
# temp_1 <- bt.intersect(a = gtf_UTR_bed, b = fiveUTR_bed, wa = T, wb = T, s= T)[,c(1:7)]
# temp_1 <- temp_1[temp_1$V4 %in% gene_need_bedtools_no_overlap,]
# temp_1$V7 <- "5TUR"
gtf_bed_temp <- gtf_bed[gtf_bed$gene_id %in% gene_need_bedtools,]
# gtf_UTR_bed_temp <- gtf_UTR_bed[gtf_UTR_bed$gene_id %in% gene_need_bedtools_overlap,]

transcript_5UTR_bed <- data.frame()
for (i in unique(gtf_bed_temp$transcript_id)) {
  # print(i)
  temp <- gtf_bed_temp[gtf_bed_temp$transcript_id==i,]
  if (temp$strand[1]=="+") {
    temp <- temp[temp$start<temp$start[temp$type=="start_codon"],]
    if (!nrow(temp)==0) {
      temp$start[1] <- temp$start[1]-500
      transcript_5UTR_bed <- rbind(transcript_5UTR_bed,temp)
    }
  } else {
    temp <- temp[temp$start>temp$start[temp$type=="start_codon"],]
    if (!nrow(temp)==0) {
      temp$end[nrow(temp)] <- temp$end[nrow(temp)]+500
      transcript_5UTR_bed <- rbind(transcript_5UTR_bed,temp)
    }
  }
}

fiveUTR_information <- as.data.frame(table(transcript_5UTR_bed$transcript_id),stringsAsFactors = F)
fiveUTR_information$UTR_length <- 0
colnames(fiveUTR_information) <- c("transcript_id","exon_number","UTR_length")
fiveUTR_information <- merge(fiveUTR_information,unique(transcript_5UTR_bed[,c("gene_id","transcript_id")]),by="transcript_id")
for (i in 1:nrow(fiveUTR_information)) {
  temp <- transcript_5UTR_bed[transcript_5UTR_bed$transcript_id==fiveUTR_information$transcript_id[i],]
  fiveUTR_information$UTR_length[i] <- sum(temp$end-temp$start)
  # print(i)
}


CRE <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/annotate/DM/bed/DM.CRE.annot.bed.gz",sep="\t")
CRE_dominant_TSS <- CRE[,c("V1","V7","V8","V4","V5","V6")]
fiveUTR_CRE_inter <- bedtoolsr::bt.intersect(a = transcript_5UTR_bed, b = CRE_dominant_TSS, wa = T, wb = T, s= T)
# fiveUTR_CRE_inter <- bedtoolsr::bt.intersect(a = transcript_5UTR_bed, b = CRE, wa = T, wb = T, s= T)
temp <- fiveUTR_CRE_inter$V5 %>% table() %>% as.data.frame(stringsAsFactors = F)
colnames(temp) <- c("transcript_id","CRE_number")
fiveUTR_information <- merge(fiveUTR_information,temp,by="transcript_id",all=T)
fiveUTR_information[is.na(fiveUTR_information)] <- 0


fiveUTR_information_new <- data.frame()
for (i in unique(fiveUTR_information$gene_id)) {
  temp <- fiveUTR_information[fiveUTR_information$gene_id==i,]
  temp <- temp[temp$CRE_number==max(temp$CRE_number),]
  temp <- temp[temp$exon_number==min(temp$exon_number),]
  temp <- temp[temp$UTR_length==max(temp$UTR_length),][1,]
  fiveUTR_information_new <- rbind(fiveUTR_information_new,temp)
}


gene_need_bedtools_final <- transcript_5UTR_bed[transcript_5UTR_bed$transcript_id %in% fiveUTR_information_new$transcript_id,1:6]
temp <- region[region$gene_id %in% gene_UTR_directly_use & region$V7=="5UTR",c(1,2,3,8,5,6)]
temp <- temp[order(temp$V2,decreasing = F),]
temp <- temp[order(temp$V1,decreasing = F),]
write.table(temp,"~/zjw/20220307_DMmouse/temp.bed",quote = F,col.names = F,row.names = F,sep = "\t")

system("~/anaconda3/bin/cgat bed2bed --method=merge --merge-by-name -I ~/zjw/20220307_DMmouse/temp.bed -S ~/zjw/20220307_DMmouse/temp_temp.bed")
temp <- read.table("~/zjw/20220307_DMmouse/temp_temp.bed",sep = "\t")
colnames(temp) <- c("V1","V2","V3","gene_id","V5","V6")


gene_UTR_directly_use_final <- data.frame()
for (i in unique(temp$gene_id)) {
  temp_temp <- temp[temp$gene_id==i,]
  if (nrow(temp_temp)==1) {
    colnames(temp_temp) <- c("V1","V2","V3","gene_id","transcript_id","V6")
  } else {
    temp_temp <- data.frame(bt.merge(i = temp_temp),gene_id=i,transcript_id=".",V6=temp_temp$V6[1])
    # write.table(temp_temp,"/home/zjw/zjw/20220307_DMmouse/temp.bed",quote = F,col.names = F,row.names = F,sep = "\t")
    # system("/home/zjw/anaconda3/bin/cgat bed2bed --method=merge --merge-by-name -I /home/zjw/zjw/20220307_DMmouse/temp.bed -S /home/zjw/zjw/20220307_DMmouse/temp_temp.bed")
    # temp_temp <- read.table("/home/zjw/zjw/20220307_DMmouse/temp_temp.bed",sep = "\t")
    # temp_temp <- temp_temp[order(temp_temp$V2,decreasing = F),]
  }
  if (temp_temp$V6[1]=="+") {
    temp_temp$V2[1] <- temp_temp$V2[1]-500
  } else {
    temp_temp$V3[nrow(temp_temp)] <- temp_temp$V3[nrow(temp_temp)]+500
  }
  gene_UTR_directly_use_final <- rbind(gene_UTR_directly_use_final,temp_temp)
}

colnames(gene_UTR_directly_use_final) <- colnames(gene_need_bedtools_final)
final <- rbind(gene_need_bedtools_final,gene_UTR_directly_use_final)
final <- final[order(final$start,decreasing = F),]
final <- final[order(final$chr,decreasing = F),]
inter <- bedtoolsr::bt.intersect(a = final, b = CRE_dominant_TSS, wa = T, wb = T, s= T)
table(inter$V4)[table(inter$V4)>=2] %>% names()

CRE_need_to_calculate <- inter$V10[inter$V4 %in% names(table(inter$V4)[table(inter$V4)>=2])]
gene_need_to_consider <- names(table(inter$V4)[table(inter$V4)>=2])
fiveUTR_length <- data.frame(CRE=CRE_need_to_calculate,length=0)
inter <- inter[inter$V10 %in% CRE_need_to_calculate & inter$V4 %in% gene_need_to_consider,]
inter$UTR_length <- 0
for (i in CRE_need_to_calculate) {
  # temp <- inter[inter$V10 %in% i & inter$V4 %in% gene_need_to_consider,]
  temp <- inter[inter$V10 %in% i,]
  final_temp <- final[final$gene_id==temp$V4,]
  if (temp$V6=="+") {
    length_1 <- temp$V3-temp$V9
    if (nrow(final_temp[final_temp$end>temp$V3,])==0) {
      length_2 <- 0
    } else {
      length_2 <- sum(final_temp[final_temp$end>temp$V3,]$end-final_temp[final_temp$end>temp$V3,]$start)
    }
  } else {
    length_1 <- temp$V9-temp$V2+1
    if (nrow(final_temp[final_temp$end<temp$V3,])==0) {
      length_2 <- 0
    } else {
      length_2 <- sum(final_temp[final_temp$end<temp$V3,]$end-final_temp[final_temp$end<temp$V3,]$start)
    }
  }
  inter$UTR_length[inter$V10==i] <- c(length_1+length_2)
  # if (!length(unique(temp$V4))==1) {
  #   print("Error")
  #   print(i)
  # }
}

temp <- inter[inter$V10 %in% CRE_need_to_calculate,]

ENSG2geneID <- read.table("~/index/mm10/gencode.vM18.annotation.sorted.gtf",sep = "\t")
ENSG2geneID <- ENSG2geneID[ENSG2geneID$V3 %in% c("gene"),]

gene_id <- strsplit(ENSG2geneID[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

gene_type <- strsplit(ENSG2geneID[,9],split = "gene_type ")
gene_type <- lapply(gene_type, function(x) x[2])
gene_type <- unlist(gene_type)
gene_type <- strsplit(gene_type,split = ";")
gene_type <- lapply(gene_type, function(x) x[1])
gene_type <- unlist(gene_type)

gene_name <- strsplit(ENSG2geneID[,9],split = "gene_name ")
gene_name <- lapply(gene_name, function(x) x[2])
gene_name <- unlist(gene_name)
gene_name <- strsplit(gene_name,split = ";")
gene_name <- lapply(gene_name, function(x) x[1])
gene_name <- unlist(gene_name)

ENSG2geneID$gene_id <- gene_id
ENSG2geneID$gene_name <- gene_name
ENSG2geneID$gene_type <- gene_type

temp <- metacell@assays$RNA@counts[rownames(metacell@assays$RNA@counts) %in% ENSG2geneID$gene_name[ENSG2geneID$gene_id %in% gene_need_to_consider],]
gene_cell_5UTRlength_mtx <- as.data.frame(metacell@assays$RNA@counts[1:length(gene_need_to_consider),stringsAsFactors=F])
rownames(gene_cell_5UTRlength_mtx) <- gene_need_to_consider
CRE_cell_5UTRlength_mtx <- as.data.frame(metacell@assays$RNA@counts[rownames(metacell@assays$RNA@counts) %in% CRE_need_to_calculate,],stringsAsFactors=F)
CRE_cell_5UTRlength_mtx_new <- data.frame()



# i=gene_need_to_consider[1]
# CRE_exp <- CRE_cell_5UTRlength_mtx[inter$V10[inter$V4==i],]
# CRE_1 <- inter$V10[inter$V4==i][1]
# CRE_2 <- inter$V10[inter$V4==i][2]
# colsum <- colSums(CRE_exp)
# (inter$UTR_length[inter$V10==CRE_1])
# (inter$UTR_length[inter$V10==CRE_2])






CRE_cell_5UTRlength_mtx_new <- data.frame()
gene_cell_5UTRlength_mtx <- data.frame()
for (i in gene_need_to_consider) {
  CRE_exp <- CRE_cell_5UTRlength_mtx[inter$V10[inter$V4==i],]
  CRE_1 <- rownames(CRE_exp)[1]
  CRE_2 <- rownames(CRE_exp)[2]
  colsum <- colSums(CRE_exp)
  temp <- as.numeric(CRE_exp[1,])/colsum*(inter$UTR_length[inter$V10==CRE_1])
  temp[is.nan(temp)] <- 0
  CRE_exp[1,] <- temp
  temp <- as.numeric(CRE_exp[2,])/colsum*(inter$UTR_length[inter$V10==CRE_2])
  temp[is.nan(temp)] <- 0
  CRE_exp[2,] <- temp
  print(colSums(CRE_exp) %>% table())
  print(i)
  CRE_cell_5UTRlength_mtx_new <- rbind(CRE_cell_5UTRlength_mtx_new,CRE_exp)
  temp <- data.frame(colSums(CRE_exp))
  colnames(temp) <- i
  temp <- t(temp)
  gene_cell_5UTRlength_mtx <- rbind(gene_cell_5UTRlength_mtx,temp)
}


gene_cell_5UTRlength_deviation_mtx <- gene_cell_5UTRlength_mtx
for (i in 1:nrow(gene_cell_5UTRlength_deviation_mtx)) {
  temp_c <- gene_cell_5UTRlength_deviation_mtx[i,] %>% as.numeric()
  mean <- temp_c[!temp_c==0] %>% mean()
  temp_c[!temp_c==0] <- temp_c[!temp_c==0]-mean
  gene_cell_5UTRlength_deviation_mtx[i,] <- temp_c
  print(i)
  print(table(temp_c))
}

cell_5UTRlength_deviation_ave <- data.frame(barcode=colnames(gene_cell_5UTRlength_deviation_mtx),mean_deviation=0)
for (i in 1:nrow(cell_5UTRlength_deviation_ave)) {
  cell_5UTRlength_deviation_ave$mean_deviation[i] <- sum(gene_cell_5UTRlength_deviation_mtx[,i])/sum(!gene_cell_5UTRlength_deviation_mtx[,i]==0)
  # print(i)
}

cell_5UTRlength_deviation_ave$mean_deviation <- scale(cell_5UTRlength_deviation_ave$mean_deviation) %>% as.numeric()
cell_5UTRlength_deviation_ave$mean_deviation[cell_5UTRlength_deviation_ave$mean_deviation>10] <- 0
cell_5UTRlength_deviation_ave$mean_deviation[cell_5UTRlength_deviation_ave$mean_deviation<(-10)] <- 0


metacell$fiveUTRlength_deviation <- cell_5UTRlength_deviation_ave$mean_deviation %>% as.numeric()

FeaturePlot(metacell,features = "fiveUTRlength_deviation")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))

scale_color_manual(c(colorRampPalette(as.character(jdb_palette("aqua_brick")))(250)))

a <- data.frame(fiveUTRlength_deviation=metacell$fiveUTRlength_deviation)
a$celltype <- rownames(a) 
a <- a[grep("CTRL",a$celltype),]
a$celltype <- gsub(pattern = "CTRL.",replacement = "",x = a$celltype)
a$celltype <- substr(a$celltype,start = 1,stop = nchar(a$celltype)-4)
a <- a[a$celltype %in% c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro"),]

# order_c <- c()
# for (i in c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")) {
#   order_c <- c(order_c,mean(a$fiveUTRlength_deviation[a$celltype==i]))
# }
# names(order_c) <- c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")
# a$celltype <- factor(a$celltype,levels = names(order_c)[order(order_c,decreasing = T)])
# 
# my_comparisons <- list()
# time <- 1
# for (i in 1:6) {
#   for (j in (i+1):7) {
#     my_comparisons[[time]] <- levels(a$celltype)[c(i,j)]
#     time <- time + 1
#   }
# }

# ggboxplot(data = a,x = "celltype",y = "fiveUTRlength_deviation",xlab = "",ylab = "fiveUTRlength_deviation",fill = "celltype",palette = jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","GMP-B"))[order(order_c,decreasing = T)],legend = "right",add.params = list(binwidth = 0.02))+
ggboxplot(data = a,x = "celltype",y = "fiveUTRlength_deviation",xlab = "",ylab = "fiveUTRlength_deviation",fill = "celltype",palette = jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","GMP-B")),legend = "right",add.params = list(binwidth = 0.02))+
  # stat_compare_means(comparisons=my_comparisons, label = "p.signif",method = "t.test")+
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=90,hjust = 1,vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(10,10),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))


gene_cell_5UTRlength_mtx[1:2,1:2]
gene_cell_5UTRlength_mtx[1:2,1:8]
# gene_cell_5UTRlength_mtx[1:2,1:]

gene_cell_5UTRlength_mtx_CTRL <- gene_cell_5UTRlength_mtx[,grep("CTRL",colnames(gene_cell_5UTRlength_mtx))]
colnames(gene_cell_5UTRlength_mtx_CTRL) <- gsub(pattern = "CTRL.",replacement = "",x = colnames(gene_cell_5UTRlength_mtx_CTRL))

temp <- data.frame(ENSG=rownames(gene_cell_5UTRlength_mtx_CTRL),geneid="A")
for (i in 1:nrow(temp)) {
  temp$geneid[i] <- ENSG2geneID$gene_name[ENSG2geneID$gene_id==temp$ENSG[i]]
}

rownames(gene_cell_5UTRlength_mtx_CTRL) <- temp$geneid


isoform_switch_bed <- data.frame()
for (i in rownames(gene_cell_5UTRlength_mtx)) {
  temp <- inter[inter$V4==i,]
  isoform_switch_bed <- rbind(isoform_switch_bed,
                              data.frame(V1=temp[1,7],V2=min(temp[,8]),V3=max(temp[,8]),V4=i,V5=ENSG2geneID$gene_name[ENSG2geneID$gene_id==i],V6=temp[1,12]))
}
isoform_switch_bed$V4 <- paste0(isoform_switch_bed$V1,"-",isoform_switch_bed$V2,"-",isoform_switch_bed$V3)
rownames(isoform_switch_bed) <- paste0(isoform_switch_bed$V1,":",isoform_switch_bed$V2,"-",isoform_switch_bed$V3)

t.test_length_CTRL <- data.frame()
for (i in 1:nrow(gene_cell_5UTRlength_mtx_CTRL)) {
  for (celltype in unique(substr(colnames(gene_cell_5UTRlength_mtx_CTRL),start = 1,stop = nchar(unique(colnames(gene_cell_5UTRlength_mtx_CTRL)))-4))[1:7]) {
    col_celltype <- unique(grep(paste0(celltype,'_0'),colnames(gene_cell_5UTRlength_mtx_CTRL)),grep(paste0(celltype,'_1'),colnames(gene_cell_5UTRlength_mtx_CTRL)))
    col_rest <- setdiff(1:ncol(gene_cell_5UTRlength_mtx_CTRL),col_celltype)
    gene_cell_5UTRlength_mtx_CTRL[i,col_celltype]
    gene_cell_5UTRlength_mtx_CTRL[i,col_rest]
    
    no_zero1 <- as.numeric(gene_cell_5UTRlength_mtx_CTRL[i,col_celltype])[as.numeric(gene_cell_5UTRlength_mtx_CTRL[i,col_celltype])!=0]
    no_zero2 <- as.numeric(gene_cell_5UTRlength_mtx_CTRL[i,col_rest])[as.numeric(gene_cell_5UTRlength_mtx_CTRL[i,col_rest])!=0]
    
    
    mean1 <- mean(no_zero1)
    mean2 <- mean(no_zero2)
    
    if (length(no_zero1)>1 & length(no_zero2)>1) {
      if (length(unique(c(no_zero1,no_zero2)))>1) {
        t.test_length_CTRL <- rbind(t.test_length_CTRL,
                                    data.frame(p.val=t.test(no_zero1,no_zero2)$p.value,
                                               ave_log2FC=log2(mean1/mean2),
                                               mean_length_1=mean1,
                                               mean_length_2=mean2,
                                               cluster=celltype,
                                               gene=rownames(gene_cell_5UTRlength_mtx_CTRL)[i]
                                    )
        )
      }
      
    }
  }
}


t.test_length_CTRL_Rod <- t.test_length_CTRL[t.test_length_CTRL$cluster=="Rod",]
t.test_length_CTRL_Rod$adj.P.Val <- p.adjust(t.test_length_CTRL_Rod$p.val,method = "BH")
t.test_length_CTRL_Rod$threshold <- "No Change"
t.test_length_CTRL_Rod$threshold[t.test_length_CTRL_Rod$adj.P.Val<0.05 & t.test_length_CTRL_Rod$ave_log2FC>0] <- "Longer"
t.test_length_CTRL_Rod$threshold[t.test_length_CTRL_Rod$adj.P.Val<0.05 & t.test_length_CTRL_Rod$ave_log2FC<0] <- "Shorter"
t.test_length_CTRL_Rod$adj.P.Val <- (-log2(t.test_length_CTRL_Rod$adj.P.Val))

t.test_length_CTRL_Cone <- t.test_length_CTRL[t.test_length_CTRL$cluster=="Cone",]
t.test_length_CTRL_Cone$adj.P.Val <- p.adjust(t.test_length_CTRL_Cone$p.val,method = "BH")
t.test_length_CTRL_Cone$threshold <- "No Change"
t.test_length_CTRL_Cone$threshold[t.test_length_CTRL_Cone$adj.P.Val<0.05 & t.test_length_CTRL_Cone$ave_log2FC>0] <- "Longer"
t.test_length_CTRL_Cone$threshold[t.test_length_CTRL_Cone$adj.P.Val<0.05 & t.test_length_CTRL_Cone$ave_log2FC<0] <- "Shorter"
t.test_length_CTRL_Cone$adj.P.Val <- (-log2(t.test_length_CTRL_Cone$adj.P.Val))

t.test_length_CTRL_Rod_BC <- t.test_length_CTRL[t.test_length_CTRL$cluster=="Rod_BC",]
t.test_length_CTRL_Rod_BC$adj.P.Val <- p.adjust(t.test_length_CTRL_Rod_BC$p.val,method = "BH")
t.test_length_CTRL_Rod_BC$threshold <- "No Change"
t.test_length_CTRL_Rod_BC$threshold[t.test_length_CTRL_Rod_BC$adj.P.Val<0.05 & t.test_length_CTRL_Rod_BC$ave_log2FC>0] <- "Longer"
t.test_length_CTRL_Rod_BC$threshold[t.test_length_CTRL_Rod_BC$adj.P.Val<0.05 & t.test_length_CTRL_Rod_BC$ave_log2FC<0] <- "Shorter"
t.test_length_CTRL_Rod_BC$adj.P.Val <- (-log2(t.test_length_CTRL_Rod_BC$adj.P.Val))

t.test_length_CTRL_Cone_BC <- t.test_length_CTRL[t.test_length_CTRL$cluster=="Cone_BC",]
t.test_length_CTRL_Cone_BC$adj.P.Val <- p.adjust(t.test_length_CTRL_Cone_BC$p.val,method = "BH")
t.test_length_CTRL_Cone_BC$threshold <- "No Change"
t.test_length_CTRL_Cone_BC$threshold[t.test_length_CTRL_Cone_BC$adj.P.Val<0.05 & t.test_length_CTRL_Cone_BC$ave_log2FC>0] <- "Longer"
t.test_length_CTRL_Cone_BC$threshold[t.test_length_CTRL_Cone_BC$adj.P.Val<0.05 & t.test_length_CTRL_Cone_BC$ave_log2FC<0] <- "Shorter"
t.test_length_CTRL_Cone_BC$adj.P.Val <- (-log2(t.test_length_CTRL_Cone_BC$adj.P.Val))

t.test_length_CTRL_MG <- t.test_length_CTRL[t.test_length_CTRL$cluster=="MG",]
t.test_length_CTRL_MG$adj.P.Val <- p.adjust(t.test_length_CTRL_MG$p.val,method = "BH")
t.test_length_CTRL_MG$threshold <- "No Change"
t.test_length_CTRL_MG$threshold[t.test_length_CTRL_MG$adj.P.Val<0.05 & t.test_length_CTRL_MG$ave_log2FC>0] <- "Longer"
t.test_length_CTRL_MG$threshold[t.test_length_CTRL_MG$adj.P.Val<0.05 & t.test_length_CTRL_MG$ave_log2FC<0] <- "Shorter"
t.test_length_CTRL_MG$adj.P.Val <- (-log2(t.test_length_CTRL_MG$adj.P.Val))

t.test_length_CTRL_Micro <- t.test_length_CTRL[t.test_length_CTRL$cluster=="Micro",]
t.test_length_CTRL_Micro$adj.P.Val <- p.adjust(t.test_length_CTRL_Micro$p.val,method = "BH")
t.test_length_CTRL_Micro$threshold <- "No Change"
t.test_length_CTRL_Micro$threshold[t.test_length_CTRL_Micro$adj.P.Val<0.05 & t.test_length_CTRL_Micro$ave_log2FC>0] <- "Longer"
t.test_length_CTRL_Micro$threshold[t.test_length_CTRL_Micro$adj.P.Val<0.05 & t.test_length_CTRL_Micro$ave_log2FC<0] <- "Shorter"
t.test_length_CTRL_Micro$adj.P.Val <- (-log2(t.test_length_CTRL_Micro$adj.P.Val))

t.test_length_CTRL_AC <- t.test_length_CTRL[t.test_length_CTRL$cluster=="AC",]
t.test_length_CTRL_AC$adj.P.Val <- p.adjust(t.test_length_CTRL_AC$p.val,method = "BH")
t.test_length_CTRL_AC$threshold <- "No Change"
t.test_length_CTRL_AC$threshold[t.test_length_CTRL_AC$adj.P.Val<0.05 & t.test_length_CTRL_AC$ave_log2FC>0] <- "Longer"
t.test_length_CTRL_AC$threshold[t.test_length_CTRL_AC$adj.P.Val<0.05 & t.test_length_CTRL_AC$ave_log2FC<0] <- "Shorter"
t.test_length_CTRL_AC$adj.P.Val <- (-log2(t.test_length_CTRL_AC$adj.P.Val))


p1 <- ggplot(data = t.test_length_CTRL_Rod, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 80))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Rod") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))



p2 <- ggplot(data = t.test_length_CTRL_Cone, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 25))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Cone") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))


p3 <- ggplot(data = t.test_length_CTRL_Rod_BC, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 40))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Rod BC") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))


p4 <- ggplot(data = t.test_length_CTRL_Cone_BC, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 65))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Cone BC") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))



p5 <- ggplot(data = t.test_length_CTRL_MG, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 51))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "MG") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))


p6 <- ggplot(data = t.test_length_CTRL_AC, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 74))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "AC") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))




p7 <- ggplot(data = t.test_length_CTRL_Micro, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 46))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Micro") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))



ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6,p7),ncol = 7)


# Mk_gene_CTRL <- readRDS("~/zjw/20220307_DMmouse/Seurat_merge_CRE/results/CTRL.marker.old.RDS")
# Mk_gene_CTRL <- Mk_gene_CTRL[Mk_gene_CTRL$p_val_adj<0.05 & abs(Mk_gene_CTRL$avg_log2FC)>=0.2,]
# 
# AC.diff[AC.diff$p_val_adj<=0.05 & AC.diff$avg_log2FC<=(-0.2)
        
vennPlot(overLapper(list("High expressed genes"=Mk_gene_CTRL$gene[Mk_gene_CTRL$cluster=="Cone BC" & Mk_gene_CTRL$avg_log2FC>0],
                         "Low expressed genes"=Mk_gene_CTRL$gene[Mk_gene_CTRL$cluster=="Cone BC" & Mk_gene_CTRL$avg_log2FC<0],
                         "Longer genes"=t.test_length_CTRL_Cone_BC$gene[t.test_length_CTRL_Cone_BC$threshold=="Longer"],
                         "Shorter genes"=t.test_length_CTRL_Cone_BC$gene[t.test_length_CTRL_Cone_BC$threshold=="Shorter"]), type="vennsets"
)
)



go_Rod_Longer <- enrichGO(t.test_length_CTRL_Rod$gene[t.test_length_CTRL_Rod$threshold=="Longer"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Rod_Shorter <- enrichGO(t.test_length_CTRL_Rod$gene[t.test_length_CTRL_Rod$threshold=="Shorter"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Rod_all <- enrichGO(t.test_length_CTRL_Rod$gene[!t.test_length_CTRL_Rod$threshold=="No change"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")

go_Cone_Longer <- enrichGO(t.test_length_CTRL_Cone$gene[t.test_length_CTRL_Cone$threshold=="Longer"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Cone_Shorter <- enrichGO(t.test_length_CTRL_Cone$gene[t.test_length_CTRL_Cone$threshold=="Shorter"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Cone_all <- enrichGO(t.test_length_CTRL_Cone$gene[!t.test_length_CTRL_Cone$threshold=="No change"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")

go_Rod_BC_Longer <- enrichGO(t.test_length_CTRL_Rod_BC$gene[t.test_length_CTRL_Rod_BC$threshold=="Longer"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Rod_BC_Shorter <- enrichGO(t.test_length_CTRL_Rod_BC$gene[t.test_length_CTRL_Rod_BC$threshold=="Shorter"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Rod_BC_all <- enrichGO(t.test_length_CTRL_Rod_BC$gene[!t.test_length_CTRL_Rod_BC$threshold=="No change"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")

go_Cone_BC_Longer <- enrichGO(t.test_length_CTRL_Cone_BC$gene[t.test_length_CTRL_Cone_BC$threshold=="Longer"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Cone_BC_Shorter <- enrichGO(t.test_length_CTRL_Cone_BC$gene[t.test_length_CTRL_Cone_BC$threshold=="Shorter"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Cone_BC_all <- enrichGO(t.test_length_CTRL_Cone_BC$gene[!t.test_length_CTRL_Cone_BC$threshold=="No change"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")

go_MG_Longer <- enrichGO(t.test_length_CTRL_MG$gene[t.test_length_CTRL_MG$threshold=="Longer"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_MG_Shorter <- enrichGO(t.test_length_CTRL_MG$gene[t.test_length_CTRL_MG$threshold=="Shorter"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_MG_all <- enrichGO(t.test_length_CTRL_MG$gene[!t.test_length_CTRL_MG$threshold=="No change"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")

go_AC_Longer <- enrichGO(t.test_length_CTRL_AC$gene[t.test_length_CTRL_AC$threshold=="Longer"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_AC_Shorter <- enrichGO(t.test_length_CTRL_AC$gene[t.test_length_CTRL_AC$threshold=="Shorter"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_AC_all <- enrichGO(t.test_length_CTRL_AC$gene[!t.test_length_CTRL_AC$threshold=="No change"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")

go_Micro_Longer <- enrichGO(t.test_length_CTRL_Micro$gene[t.test_length_CTRL_Micro$threshold=="Longer"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Micro_Shorter <- enrichGO(t.test_length_CTRL_Micro$gene[t.test_length_CTRL_Micro$threshold=="Shorter"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
go_Micro_all <- enrichGO(t.test_length_CTRL_Micro$gene[!t.test_length_CTRL_Micro$threshold=="No change"],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")


dotplot(go_AC_all,
        title='Top5 GO terms of Gene cluster 1',
        showCategory=5,split='ONTOLOGY')+
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))







t.test_length_CTRL_temp <- rbind(t.test_length_CTRL_Rod,
                                 t.test_length_CTRL_Cone,
                                 t.test_length_CTRL_Rod_BC,
                                 t.test_length_CTRL_Cone_BC,
                                 t.test_length_CTRL_MG,
                                 t.test_length_CTRL_AC,
                                 t.test_length_CTRL_Micro)




hm_up <- t.test_length_CTRL_temp[t.test_length_CTRL_temp$threshold=="Longer",c("gene","cluster")]



hm_up_new <- data.frame(gene=unique(hm_up$gene))
hm_up_new[,2:8] <- 0
hm_up_new <- hm_up_new[,-1]
rownames(hm_up_new) <- unique(hm_up$gene)
colnames(hm_up_new) <- c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")

for (i in 1:nrow(hm_up)) {
  hm_up_new[rownames(hm_up_new)==hm_up$gene[i],colnames(hm_up_new)==hm_up$cluster[i]] <- 1
}

hm_up_new <-rbind(hm_up_new[rowSums(hm_up_new)==7,],
                  hm_up_new[rowSums(hm_up_new)==6,],
                  hm_up_new[rowSums(hm_up_new)==5,],
                  hm_up_new[rowSums(hm_up_new)==4,],
                  hm_up_new[rowSums(hm_up_new)==3,],
                  hm_up_new[rowSums(hm_up_new)==2,],
                  hm_up_new[rowSums(hm_up_new)==1,])

pheatmap(hm_up_new,
         color = colorRampPalette(c("grey", "red"))(5),cluster_cols = F,cluster_rows = F,border = F,show_rownames = T)









gene_cell_5UTRlength_mtx_all <- gene_cell_5UTRlength_mtx

temp <- data.frame(ENSG=rownames(gene_cell_5UTRlength_mtx_all),geneid="A")
for (i in 1:nrow(temp)) {
  temp$geneid[i] <- ENSG2geneID$gene_name[ENSG2geneID$gene_id==temp$ENSG[i]]
}

rownames(gene_cell_5UTRlength_mtx_all) <- temp$geneid


t.test_length <- data.frame()
for (i in 1:nrow(gene_cell_5UTRlength_mtx_all)) {
  for (celltype in c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")) {
    col_celltype <- unique(grep(paste0("DM.",celltype,'_0'),colnames(gene_cell_5UTRlength_mtx_all)),grep("DM.",paste0(celltype,'_1'),colnames(gene_cell_5UTRlength_mtx_all)))
    col_rest <- unique(grep(paste0("CTRL.",celltype,'_0'),colnames(gene_cell_5UTRlength_mtx_all)),grep("CTRL.",paste0(celltype,'_1'),colnames(gene_cell_5UTRlength_mtx_all)))
    gene_cell_5UTRlength_mtx_all[i,col_celltype]
    gene_cell_5UTRlength_mtx_all[i,col_rest]
    
    no_zero1 <- as.numeric(gene_cell_5UTRlength_mtx_all[i,col_celltype])[as.numeric(gene_cell_5UTRlength_mtx_all[i,col_celltype])!=0]
    no_zero2 <- as.numeric(gene_cell_5UTRlength_mtx_all[i,col_rest])[as.numeric(gene_cell_5UTRlength_mtx_all[i,col_rest])!=0]
    
    
    mean1 <- mean(no_zero1)
    mean2 <- mean(no_zero2)
    
    if (length(no_zero1)>1 & length(no_zero2)>1) {
      if (length(unique(c(no_zero1,no_zero2)))>1) {
        t.test_length <- rbind(t.test_length,
                               data.frame(p.val=t.test(no_zero1,no_zero2)$p.value,
                                          ave_log2FC=log2(mean1/mean2),
                                          mean_length_1=mean1,
                                          mean_length_2=mean2,
                                          cluster=celltype,
                                          gene=rownames(gene_cell_5UTRlength_mtx_all)[i]
                               )
        )
      }
      
    }
  }
}




t.test_length_Rod <- t.test_length[t.test_length$cluster=="Rod",]
t.test_length_Rod$adj.P.Val <- p.adjust(t.test_length_Rod$p.val,method = "hochberg")
t.test_length_Rod$threshold <- "No Change"
t.test_length_Rod$threshold[t.test_length_Rod$adj.P.Val<0.05 & t.test_length_Rod$ave_log2FC>0] <- "Longer"
t.test_length_Rod$threshold[t.test_length_Rod$adj.P.Val<0.05 & t.test_length_Rod$ave_log2FC<0] <- "Shorter"
t.test_length_Rod$adj.P.Val <- (-log2(t.test_length_Rod$adj.P.Val))
t.test_length_Rod$p.val <- (-log2(t.test_length_Rod$p.val))

t.test_length_Cone <- t.test_length[t.test_length$cluster=="Cone",]
t.test_length_Cone$adj.P.Val <- p.adjust(t.test_length_Cone$p.val,method = "hochberg")
t.test_length_Cone$threshold <- "No Change"
t.test_length_Cone$threshold[t.test_length_Cone$adj.P.Val<0.05 & t.test_length_Cone$ave_log2FC>0] <- "Longer"
t.test_length_Cone$threshold[t.test_length_Cone$adj.P.Val<0.05 & t.test_length_Cone$ave_log2FC<0] <- "Shorter"
t.test_length_Cone$adj.P.Val <- (-log2(t.test_length_Cone$adj.P.Val))
t.test_length_Cone$p.val <- (-log2(t.test_length_Cone$p.val))


t.test_length_Rod_BC <- t.test_length[t.test_length$cluster=="Rod_BC",]
t.test_length_Rod_BC$adj.P.Val <- p.adjust(t.test_length_Rod_BC$p.val,method = "hochberg")
t.test_length_Rod_BC$threshold <- "No Change"
t.test_length_Rod_BC$threshold[t.test_length_Rod_BC$adj.P.Val<0.05 & t.test_length_Rod_BC$ave_log2FC>0] <- "Longer"
t.test_length_Rod_BC$threshold[t.test_length_Rod_BC$adj.P.Val<0.05 & t.test_length_Rod_BC$ave_log2FC<0] <- "Shorter"
t.test_length_Rod_BC$adj.P.Val <- (-log2(t.test_length_Rod_BC$adj.P.Val))
t.test_length_Rod_BC$p.val <- (-log2(t.test_length_Rod_BC$p.val))

t.test_length_Cone_BC <- t.test_length[t.test_length$cluster=="Cone_BC",]
t.test_length_Cone_BC$adj.P.Val <- p.adjust(t.test_length_Cone_BC$p.val,method = "hochberg")
t.test_length_Cone_BC$threshold <- "No Change"
t.test_length_Cone_BC$threshold[t.test_length_Cone_BC$adj.P.Val<0.05 & t.test_length_Cone_BC$ave_log2FC>0] <- "Longer"
t.test_length_Cone_BC$threshold[t.test_length_Cone_BC$adj.P.Val<0.05 & t.test_length_Cone_BC$ave_log2FC<0] <- "Shorter"
t.test_length_Cone_BC$adj.P.Val <- (-log2(t.test_length_Cone_BC$adj.P.Val))
t.test_length_Cone_BC$p.val <- (-log2(t.test_length_Cone_BC$p.val))

t.test_length_MG <- t.test_length[t.test_length$cluster=="MG",]
t.test_length_MG$adj.P.Val <- p.adjust(t.test_length_MG$p.val,method = "hochberg")
t.test_length_MG$threshold <- "No Change"
t.test_length_MG$threshold[t.test_length_MG$adj.P.Val<0.05 & t.test_length_MG$ave_log2FC>0] <- "Longer"
t.test_length_MG$threshold[t.test_length_MG$adj.P.Val<0.05 & t.test_length_MG$ave_log2FC<0] <- "Shorter"
t.test_length_MG$adj.P.Val <- (-log2(t.test_length_MG$adj.P.Val))
t.test_length_MG$p.val <- (-log2(t.test_length_MG$p.val))

t.test_length_Micro <- t.test_length[t.test_length$cluster=="Micro",]
t.test_length_Micro$adj.P.Val <- p.adjust(t.test_length_Micro$p.val,method = "hochberg")
t.test_length_Micro$threshold <- "No Change"
t.test_length_Micro$threshold[t.test_length_Micro$adj.P.Val<0.05 & t.test_length_Micro$ave_log2FC>0] <- "Longer"
t.test_length_Micro$threshold[t.test_length_Micro$adj.P.Val<0.05 & t.test_length_Micro$ave_log2FC<0] <- "Shorter"
t.test_length_Micro$adj.P.Val <- (-log2(t.test_length_Micro$adj.P.Val))
t.test_length_Micro$p.val <- (-log2(t.test_length_Micro$p.val))

t.test_length_AC <- t.test_length[t.test_length$cluster=="AC",]
t.test_length_AC$adj.P.Val <- p.adjust(t.test_length_AC$p.val,method = "hochberg")
t.test_length_AC$threshold <- "No Change"
t.test_length_AC$threshold[t.test_length_AC$adj.P.Val<0.05 & t.test_length_AC$ave_log2FC>0] <- "Longer"
t.test_length_AC$threshold[t.test_length_AC$adj.P.Val<0.05 & t.test_length_AC$ave_log2FC<0] <- "Shorter"
t.test_length_AC$adj.P.Val <- (-log2(t.test_length_AC$adj.P.Val))
t.test_length_AC$p.val <- (-log2(t.test_length_AC$p.val))



t.test_length_temp <- rbind(t.test_length_Rod,
                            t.test_length_Cone,
                            t.test_length_Rod_BC,
                            t.test_length_Cone_BC,
                            t.test_length_MG,
                            t.test_length_AC,
                            t.test_length_Micro)

t.test_length_temp$cluster <- factor(t.test_length_temp$cluster,levels = c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro"))

ggplot(data = t.test_length_temp,aes(x = cluster, y = ave_log2FC, color = threshold))+
  geom_jitter(size = 0.85,
              width =0.4)+
  scale_color_manual(values=c("#ED0000","grey","#00468B")) +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))
  

# 
# hm_up <- t.test_length_temp[t.test_length_temp$threshold=="Shorter",c("gene","cluster")]
# 
# 
# 
# hm_up_new <- data.frame(gene=unique(hm_up$gene))
# hm_up_new[,2:8] <- 0
# hm_up_new <- hm_up_new[,-1]
# rownames(hm_up_new) <- unique(hm_up$gene)
# colnames(hm_up_new) <- c("Rod","Cone","Rod_BC","Cone_BC","MG","AC","Micro")
# 
# for (i in 1:nrow(hm_up)) {
#   hm_up_new[rownames(hm_up_new)==hm_up$gene[i],colnames(hm_up_new)==hm_up$cluster[i]] <- 1
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
#          color = colorRampPalette(c("grey", "red"))(5),cluster_cols = F,cluster_rows = F,border = F,show_rownames = F)



p1 <- ggplot(data = t.test_length_Rod, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 30))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Rod") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))



p2 <- ggplot(data = t.test_length_Cone, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 25))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Cone") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))


p3 <- ggplot(data = t.test_length_Rod_BC, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 40))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Rod BC") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))


p4 <- ggplot(data = t.test_length_Cone_BC, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 65))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Cone BC") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))



p5 <- ggplot(data = t.test_length_MG, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 51))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "MG") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))


p6 <- ggplot(data = t.test_length_AC, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 74))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "AC") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))




p7 <- ggplot(data = t.test_length_Micro, aes(x = ave_log2FC, y = adj.P.Val, color = threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","grey","#ED0000")) +
  xlim(c(-2, 2)) +
  ylim(c(0, 46))+
  geom_vline(xintercept = c(0),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=2,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)",title = "Micro") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position=c(10,10),
        legend.title = element_blank(),
        legend.text = element_text(size = 11))



ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6,p7),ncol = 7)



a <- read.table("~/zjw/20220307_DMmouse/gwas/list_gwas_summary_statistics.tsv",header = 1,sep = "\t")


gwas <- read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0001365-associations-2022-12-2.csv") # AMD
gwas <- read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0003770-associations-2022-12-2.csv") # diabetic retinopathy
gwas <- read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0004190-associations-2022-12-2.csv") # primary open-angle glaucoma
gwas <- read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0009322-associations-2022-12-2.csv") # proliferative diabetic retinopathy
gwas <- read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_1001506-associations-2022-12-2.csv") # primary-angle closure glaucoma
gwas <- read.csv("~/zjw/20220307_DMmouse/efotraits_Orphanet_791-associations-2022-12-2.csv") # Retinitis pigmentosa

gwas <- rbind(data.frame(gene=read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0001365-associations-2022-12-2.csv")$Mapped.gene %>% strsplit(", ") %>% unlist(),disease="AMD"),
              data.frame(gene=read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0003770-associations-2022-12-2.csv")$Mapped.gene %>% strsplit(", ") %>% unlist(),disease="DR"),
              data.frame(gene=read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0004190-associations-2022-12-2.csv")$Mapped.gene %>% strsplit(", ") %>% unlist(),disease="POAG"),
              data.frame(gene=read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_0009322-associations-2022-12-2.csv")$Mapped.gene %>% strsplit(", ") %>% unlist(),disease="PDR"),
              data.frame(gene=read.csv("~/zjw/20220307_DMmouse/efotraits_EFO_1001506-associations-2022-12-2.csv")$Mapped.gene %>% strsplit(", ") %>% unlist(),disease="PACG"),
              data.frame(gene=read.csv("~/zjw/20220307_DMmouse/efotraits_Orphanet_791-associations-2022-12-2.csv")$Mapped.gene %>% strsplit(", ") %>% unlist(),disease="RP")
              )

# gwas <- gwas[,c("Mapped.gene","disease")]
gwas <- gwas[!gwas$gene=="'-",]


write.table(CRE[,c(1,2,3)],"all_peaks.bed",quote = F,row.names = F,col.names = F,sep = "\t")

samplesheet <- data.frame()
for (i in unique(mk.CTRL$cluster)) {
  write.table(CRE[CRE$V4 %in% mk.CTRL$gene[mk.CTRL$cluster==i],c(1,2,3)],paste0(gsub(pattern = " ",replacement = "_",x = i),"_peaks.bed"),quote = F,row.names = F,col.names = F,sep = "\t")
  
  samplesheet <- rbind(samplesheet,
                       data.frame(sample_id=i,sites=paste0(gsub(pattern = " ",replacement = "_",x = i),"_peaks.bed")))
}

write.table(CRE[,c(1,2,3)],"all_peaks.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(samplesheet,"samplesheet.txt",quote = F,row.names = F,col.names = T,sep = "\t")


write.table(data.frame(sample_id=mk.CTRL$cluster,sites=mk.CTRL$gene),"samplesheet.txt",quote = F,row.names = F,col.names = T,sep = "\t")

gwas <- read.table("~/zjw/20220307_DMmouse/gwas/GCST90013904_buildGRCh38.tsv.gz",sep = "\t",header = 1)  # 403833
gwas$chromosome <- paste0("chr",gwas$chromosome)
gwas$ID <- paste(gwas$chromosome,gwas$base_pair_location,gwas$effect_allele,gwas$other_allele)
write.table(data.frame(gwas$chromosome,gwas$base_pair_location,gwas$base_pair_location,gwas$effect_allele,gwas$other_allele),"ex1.avinput",quote = F,row.names = F,col.names = F,sep = "\t")

system("~/biosoft/annovar/annotate_variation.pl ex1.avinput ~/biosoft/annovar/humandb/ -filter -build hg38 -dbtype avsnp150")

gwas_new <- read.table("~/zjw/20220307_DMmouse/ex1.avinput.hg38_avsnp150_dropped",sep = "\t")
gwas_new$ID <- paste(gwas_new$V3,gwas_new$V4,gwas_new$V6,gwas_new$V7)
gwas_new <- gwas_new[order(gwas_new$ID),]

gwas <- gwas[gwas$ID %in% gwas_new$ID,]
gwas <- gwas[order(gwas$ID),]
gwas$ID <- gwas_new$V2



# write.table(data.frame(gwas$chromosome,gwas$base_pair_location-1,gwas$base_pair_location,gwas$ID),"temp.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# 
# system("/disk2/user/jizhon/anaconda3/bin/liftOver ~/zjw/20220307_DMmouse/temp.bed ~/index/liftover/hg38ToMm10.over.chain.gz ~/zjw/20220307_DMmouse/temp_mm10.bed ~/zjw/20220307_DMmouse/temp_mm10_umapped.bed")
# 
# 
# gwas_mm10 <- read.table("~/zjw/20220307_DMmouse/temp_mm10.bed")
# gwas_mm10 <- gwas_mm10[order(gwas_mm10$V4),]
# 
# 
# gwas <- gwas[gwas$ID %in% gwas_mm10$V4,]
# gwas <- gwas[order(gwas$ID),]
# gwas$chromosome <- gwas_mm10$V1
# gwas$base_pair_location <- gwas_mm10$V3
colnames(gwas) <- c("hg18chr","bp","pval","a1","a2","beta","se","snpid")
gwas <- gwas[,c("snpid","hg18chr","bp","a1","a2","beta","se","pval")]
write.table(gwas,"temp.txt",quote = F,row.names = F,col.names = T,sep = "\t")

# run in sh
# ~/biosoft/ldsc/munge_sumstats.py --sumstats temp.txt --N 403833 --out temp_new" 

