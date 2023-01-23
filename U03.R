options(stringsAsFactors = FALSE)
options(scipen = 100)
options(bedtools.path = "~/anaconda3/bin/")
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/disk2/user/jizhon/biosoft/homer/bin/", sep=":"))


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
library(bedtoolsr)

region <- read.table("~/zjw/annotation/gencode_mm10_all_gene_all_transcript_region.bed",sep = "\t")

enhancer <- read.table("~/zjw/20220307_DMmouse/F5.mm10.enhancers.bed",sep = "\t")[,1:5]
enhancer <- rbind(data.frame(enhancer,V6="+",V7="Enhancer"),data.frame(enhancer,V6="-",V7="Enhancer"))
region <- rbind(region,enhancer)




### DM1
temp_tss <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443300/annotate/SRR12443300/bed/SRR12443300.CRE.coord.bed.gz",sep = "\t")
temp_tss <- temp_tss[,c(1,7,8,4,5,6)]

temp_tss_region <- bedtoolsr::bt.intersect(a = temp_tss, b = region, wa = T, wb = T, s= T)



for(i in 1:nrow(temp_tss)) {
  a <- temp_tss_region[temp_tss_region[,1]==temp_tss[i,1] & temp_tss_region[,3]==temp_tss[i,3] & temp_tss_region[,6]==temp_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000","Enhancer"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  temp_tss[i,7] <- b
}
temp_tss[is.na(temp_tss[,7]),7] <- "intergenic"
temp_tss[temp_tss[,7]=="TES1000",7] <- "intergenic"
temp_tss[temp_tss[,7]=="last_exon",7] <- "other_exon"
temp_tss[temp_tss[,7]=="last_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="one_exon",7] <- "first_exon"
temp_tss[temp_tss[,7]=="one_intron",7] <- "first_intron"
temp_tss[temp_tss[,7]=="first_exon",7] <- "first exon"
temp_tss[temp_tss[,7]=="first_intron",7] <- "first intron"
temp_tss[temp_tss[,7]=="other_exon",7] <- "other exon"
temp_tss[temp_tss[,7]=="other_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="5UTR",7] <- "5'UTR"
temp_tss[temp_tss[,7]=="3UTR",7] <- "3'UTR"
temp_tss[temp_tss[,7]=="TSS1000",7] <- "TSS±1000"
temp_tss_df <- as.data.frame(table(temp_tss[,7]))
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c("TSS±1000", "Enhancer", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]
temp_tss_df[,1] <- paste(temp_tss_df[,1]," (",round(temp_tss_df[,2]/sum(temp_tss_df[,2])*100,2),"%)",sep = "")
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c(grep("TSS",temp_tss_df[,1],value = T),
                                                      grep("Enhancer",temp_tss_df[,1],value = T),
                                                      grep("5'UTR",temp_tss_df[,1],value = T),
                                                      grep("first exon",temp_tss_df[,1],value = T),
                                                      grep("first intron",temp_tss_df[,1],value = T),
                                                      grep("other exon",temp_tss_df[,1],value = T),
                                                      grep("other intron",temp_tss_df[,1],value = T),
                                                      grep("3'UTR",temp_tss_df[,1],value = T),
                                                      grep("intergenic",temp_tss_df[,1],value = T)), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]




p1 <- ggplot(temp_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL,title = "DM1") +
  scale_x_continuous(breaks = NULL)


write.table(temp_tss[!temp_tss$V7.1=="TSS±1000",1:6],"temp.bed",quote = F,row.names = F,col.names = F,sep = "\t")

system("findMotifsGenome.pl temp.bed mm10 homer_DM1")

### DM2
temp_tss <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443301/annotate/SRR12443301/bed/SRR12443301.CRE.coord.bed.gz",sep = "\t")
temp_tss <- temp_tss[,c(1,7,8,4,5,6)]

temp_tss_region <- bedtoolsr::bt.intersect(a = temp_tss, b = region, wa = T, wb = T, s= T)



for(i in 1:nrow(temp_tss)) {
  a <- temp_tss_region[temp_tss_region[,1]==temp_tss[i,1] & temp_tss_region[,3]==temp_tss[i,3] & temp_tss_region[,6]==temp_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000","Enhancer"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  temp_tss[i,7] <- b
}
temp_tss[is.na(temp_tss[,7]),7] <- "intergenic"
temp_tss[temp_tss[,7]=="TES1000",7] <- "intergenic"
temp_tss[temp_tss[,7]=="last_exon",7] <- "other_exon"
temp_tss[temp_tss[,7]=="last_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="one_exon",7] <- "first_exon"
temp_tss[temp_tss[,7]=="one_intron",7] <- "first_intron"
temp_tss[temp_tss[,7]=="first_exon",7] <- "first exon"
temp_tss[temp_tss[,7]=="first_intron",7] <- "first intron"
temp_tss[temp_tss[,7]=="other_exon",7] <- "other exon"
temp_tss[temp_tss[,7]=="other_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="5UTR",7] <- "5'UTR"
temp_tss[temp_tss[,7]=="3UTR",7] <- "3'UTR"
temp_tss[temp_tss[,7]=="TSS1000",7] <- "TSS±1000"
temp_tss_df <- as.data.frame(table(temp_tss[,7]))
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c("TSS±1000", "Enhancer", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]
temp_tss_df[,1] <- paste(temp_tss_df[,1]," (",round(temp_tss_df[,2]/sum(temp_tss_df[,2])*100,2),"%)",sep = "")
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c(grep("TSS",temp_tss_df[,1],value = T),
                                                      grep("Enhancer",temp_tss_df[,1],value = T),
                                                      grep("5'UTR",temp_tss_df[,1],value = T),
                                                      grep("first exon",temp_tss_df[,1],value = T),
                                                      grep("first intron",temp_tss_df[,1],value = T),
                                                      grep("other exon",temp_tss_df[,1],value = T),
                                                      grep("other intron",temp_tss_df[,1],value = T),
                                                      grep("3'UTR",temp_tss_df[,1],value = T),
                                                      grep("intergenic",temp_tss_df[,1],value = T)), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]




p2 <- ggplot(temp_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL,title = "DM2") +
  scale_x_continuous(breaks = NULL)

write.table(temp_tss[!temp_tss$V7.1=="TSS±1000",1:6],"temp.bed",quote = F,row.names = F,col.names = F,sep = "\t")

system("findMotifsGenome.pl temp.bed mm10 homer_DM2")


### DM3
temp_tss <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443302/annotate/SRR12443302/bed/SRR12443302.CRE.coord.bed.gz",sep = "\t")
temp_tss <- temp_tss[,c(1,7,8,4,5,6)]

temp_tss_region <- bedtoolsr::bt.intersect(a = temp_tss, b = region, wa = T, wb = T, s= T)



for(i in 1:nrow(temp_tss)) {
  a <- temp_tss_region[temp_tss_region[,1]==temp_tss[i,1] & temp_tss_region[,3]==temp_tss[i,3] & temp_tss_region[,6]==temp_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000","Enhancer"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  temp_tss[i,7] <- b
}
temp_tss[is.na(temp_tss[,7]),7] <- "intergenic"
temp_tss[temp_tss[,7]=="TES1000",7] <- "intergenic"
temp_tss[temp_tss[,7]=="last_exon",7] <- "other_exon"
temp_tss[temp_tss[,7]=="last_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="one_exon",7] <- "first_exon"
temp_tss[temp_tss[,7]=="one_intron",7] <- "first_intron"
temp_tss[temp_tss[,7]=="first_exon",7] <- "first exon"
temp_tss[temp_tss[,7]=="first_intron",7] <- "first intron"
temp_tss[temp_tss[,7]=="other_exon",7] <- "other exon"
temp_tss[temp_tss[,7]=="other_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="5UTR",7] <- "5'UTR"
temp_tss[temp_tss[,7]=="3UTR",7] <- "3'UTR"
temp_tss[temp_tss[,7]=="TSS1000",7] <- "TSS±1000"
temp_tss_df <- as.data.frame(table(temp_tss[,7]))
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c("TSS±1000", "Enhancer", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]
temp_tss_df[,1] <- paste(temp_tss_df[,1]," (",round(temp_tss_df[,2]/sum(temp_tss_df[,2])*100,2),"%)",sep = "")
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c(grep("TSS",temp_tss_df[,1],value = T),
                                                      grep("Enhancer",temp_tss_df[,1],value = T),
                                                      grep("5'UTR",temp_tss_df[,1],value = T),
                                                      grep("first exon",temp_tss_df[,1],value = T),
                                                      grep("first intron",temp_tss_df[,1],value = T),
                                                      grep("other exon",temp_tss_df[,1],value = T),
                                                      grep("other intron",temp_tss_df[,1],value = T),
                                                      grep("3'UTR",temp_tss_df[,1],value = T),
                                                      grep("intergenic",temp_tss_df[,1],value = T)), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]




p3 <- ggplot(temp_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL,title = "DM3") +
  scale_x_continuous(breaks = NULL)

write.table(temp_tss[!temp_tss$V7.1=="TSS±1000",1:6],"temp.bed",quote = F,row.names = F,col.names = F,sep = "\t")

system("findMotifsGenome.pl temp.bed mm10 homer_DM3")


### CTRL1
temp_tss <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443303/annotate/SRR12443303/bed/SRR12443303.CRE.coord.bed.gz",sep = "\t")
temp_tss <- temp_tss[,c(1,7,8,4,5,6)]

temp_tss_region <- bedtoolsr::bt.intersect(a = temp_tss, b = region, wa = T, wb = T, s= T)



for(i in 1:nrow(temp_tss)) {
  a <- temp_tss_region[temp_tss_region[,1]==temp_tss[i,1] & temp_tss_region[,3]==temp_tss[i,3] & temp_tss_region[,6]==temp_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000","Enhancer"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  temp_tss[i,7] <- b
}
temp_tss[is.na(temp_tss[,7]),7] <- "intergenic"
temp_tss[temp_tss[,7]=="TES1000",7] <- "intergenic"
temp_tss[temp_tss[,7]=="last_exon",7] <- "other_exon"
temp_tss[temp_tss[,7]=="last_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="one_exon",7] <- "first_exon"
temp_tss[temp_tss[,7]=="one_intron",7] <- "first_intron"
temp_tss[temp_tss[,7]=="first_exon",7] <- "first exon"
temp_tss[temp_tss[,7]=="first_intron",7] <- "first intron"
temp_tss[temp_tss[,7]=="other_exon",7] <- "other exon"
temp_tss[temp_tss[,7]=="other_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="5UTR",7] <- "5'UTR"
temp_tss[temp_tss[,7]=="3UTR",7] <- "3'UTR"
temp_tss[temp_tss[,7]=="TSS1000",7] <- "TSS±1000"
temp_tss_df <- as.data.frame(table(temp_tss[,7]))
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c("TSS±1000", "Enhancer", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]
temp_tss_df[,1] <- paste(temp_tss_df[,1]," (",round(temp_tss_df[,2]/sum(temp_tss_df[,2])*100,2),"%)",sep = "")
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c(grep("TSS",temp_tss_df[,1],value = T),
                                                      grep("Enhancer",temp_tss_df[,1],value = T),
                                                      grep("5'UTR",temp_tss_df[,1],value = T),
                                                      grep("first exon",temp_tss_df[,1],value = T),
                                                      grep("first intron",temp_tss_df[,1],value = T),
                                                      grep("other exon",temp_tss_df[,1],value = T),
                                                      grep("other intron",temp_tss_df[,1],value = T),
                                                      grep("3'UTR",temp_tss_df[,1],value = T),
                                                      grep("intergenic",temp_tss_df[,1],value = T)), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]




p4 <- ggplot(temp_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL,title = "CTRL1") +
  scale_x_continuous(breaks = NULL)

write.table(temp_tss[!temp_tss$V7.1=="TSS±1000",1:6],"temp.bed",quote = F,row.names = F,col.names = F,sep = "\t")

system("findMotifsGenome.pl temp.bed mm10 homer_CTRL1")


### CTRL2
temp_tss <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443304/annotate/SRR12443304/bed/SRR12443304.CRE.coord.bed.gz",sep = "\t")
temp_tss <- temp_tss[,c(1,7,8,4,5,6)]

temp_tss_region <- bedtoolsr::bt.intersect(a = temp_tss, b = region, wa = T, wb = T, s= T)



for(i in 1:nrow(temp_tss)) {
  a <- temp_tss_region[temp_tss_region[,1]==temp_tss[i,1] & temp_tss_region[,3]==temp_tss[i,3] & temp_tss_region[,6]==temp_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000","Enhancer"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  temp_tss[i,7] <- b
}
temp_tss[is.na(temp_tss[,7]),7] <- "intergenic"
temp_tss[temp_tss[,7]=="TES1000",7] <- "intergenic"
temp_tss[temp_tss[,7]=="last_exon",7] <- "other_exon"
temp_tss[temp_tss[,7]=="last_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="one_exon",7] <- "first_exon"
temp_tss[temp_tss[,7]=="one_intron",7] <- "first_intron"
temp_tss[temp_tss[,7]=="first_exon",7] <- "first exon"
temp_tss[temp_tss[,7]=="first_intron",7] <- "first intron"
temp_tss[temp_tss[,7]=="other_exon",7] <- "other exon"
temp_tss[temp_tss[,7]=="other_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="5UTR",7] <- "5'UTR"
temp_tss[temp_tss[,7]=="3UTR",7] <- "3'UTR"
temp_tss[temp_tss[,7]=="TSS1000",7] <- "TSS±1000"
temp_tss_df <- as.data.frame(table(temp_tss[,7]))
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c("TSS±1000", "Enhancer", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]
temp_tss_df[,1] <- paste(temp_tss_df[,1]," (",round(temp_tss_df[,2]/sum(temp_tss_df[,2])*100,2),"%)",sep = "")
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c(grep("TSS",temp_tss_df[,1],value = T),
                                                      grep("Enhancer",temp_tss_df[,1],value = T),
                                                      grep("5'UTR",temp_tss_df[,1],value = T),
                                                      grep("first exon",temp_tss_df[,1],value = T),
                                                      grep("first intron",temp_tss_df[,1],value = T),
                                                      grep("other exon",temp_tss_df[,1],value = T),
                                                      grep("other intron",temp_tss_df[,1],value = T),
                                                      grep("3'UTR",temp_tss_df[,1],value = T),
                                                      grep("intergenic",temp_tss_df[,1],value = T)), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]




p5 <- ggplot(temp_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL,title = "CTRL2") +
  scale_x_continuous(breaks = NULL)

write.table(temp_tss[!temp_tss$V7.1=="TSS±1000",1:6],"temp.bed",quote = F,row.names = F,col.names = F,sep = "\t")

system("findMotifsGenome.pl temp.bed mm10 homer_CTRL2")


### CTRL3
temp_tss <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/SRR12443305/annotate/SRR12443305/bed/SRR12443305.CRE.coord.bed.gz",sep = "\t")
temp_tss <- temp_tss[,c(1,7,8,4,5,6)]

temp_tss_region <- bedtoolsr::bt.intersect(a = temp_tss, b = region, wa = T, wb = T, s= T)



for(i in 1:nrow(temp_tss)) {
  a <- temp_tss_region[temp_tss_region[,1]==temp_tss[i,1] & temp_tss_region[,3]==temp_tss[i,3] & temp_tss_region[,6]==temp_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000","Enhancer"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  temp_tss[i,7] <- b
}
temp_tss[is.na(temp_tss[,7]),7] <- "intergenic"
temp_tss[temp_tss[,7]=="TES1000",7] <- "intergenic"
temp_tss[temp_tss[,7]=="last_exon",7] <- "other_exon"
temp_tss[temp_tss[,7]=="last_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="one_exon",7] <- "first_exon"
temp_tss[temp_tss[,7]=="one_intron",7] <- "first_intron"
temp_tss[temp_tss[,7]=="first_exon",7] <- "first exon"
temp_tss[temp_tss[,7]=="first_intron",7] <- "first intron"
temp_tss[temp_tss[,7]=="other_exon",7] <- "other exon"
temp_tss[temp_tss[,7]=="other_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="5UTR",7] <- "5'UTR"
temp_tss[temp_tss[,7]=="3UTR",7] <- "3'UTR"
temp_tss[temp_tss[,7]=="TSS1000",7] <- "TSS±1000"
temp_tss_df <- as.data.frame(table(temp_tss[,7]))
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c("TSS±1000", "Enhancer", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]
temp_tss_df[,1] <- paste(temp_tss_df[,1]," (",round(temp_tss_df[,2]/sum(temp_tss_df[,2])*100,2),"%)",sep = "")
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c(grep("TSS",temp_tss_df[,1],value = T),
                                                      grep("Enhancer",temp_tss_df[,1],value = T),
                                                      grep("5'UTR",temp_tss_df[,1],value = T),
                                                      grep("first exon",temp_tss_df[,1],value = T),
                                                      grep("first intron",temp_tss_df[,1],value = T),
                                                      grep("other exon",temp_tss_df[,1],value = T),
                                                      grep("other intron",temp_tss_df[,1],value = T),
                                                      grep("3'UTR",temp_tss_df[,1],value = T),
                                                      grep("intergenic",temp_tss_df[,1],value = T)), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]




p6 <- ggplot(temp_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL,title = "CTRL3") +
  scale_x_continuous(breaks = NULL)

write.table(temp_tss[!temp_tss$V7.1=="TSS±1000",1:6],"temp.bed",quote = F,row.names = F,col.names = F,sep = "\t")

system("findMotifsGenome.pl temp.bed mm10 homer_CTRL3")


a <- rbind(data.frame(fread("~/zjw/20220307_DMmouse/homer_DM1/knownResults.txt",sep = "\t"),type="DM1"),
           data.frame(fread("~/zjw/20220307_DMmouse/homer_DM2/knownResults.txt",sep = "\t"),type="DM2"),
           data.frame(fread("~/zjw/20220307_DMmouse/homer_DM3/knownResults.txt",sep = "\t"),type="DM3"),
           data.frame(fread("~/zjw/20220307_DMmouse/homer_CTRL1/knownResults.txt",sep = "\t"),type="CTRL1"),
           data.frame(fread("~/zjw/20220307_DMmouse/homer_CTRL2/knownResults.txt",sep = "\t"),type="CTRL2"),
           data.frame(fread("~/zjw/20220307_DMmouse/homer_CTRL3/knownResults.txt",sep = "\t"),type="CTRL3")
           )


a$P.value <- as.numeric(a$P.value)

a <- a[grep("Promoter",a$Motif.Name),]
a <- a[grep("ox",a$Motif.Name),]

write.csv(a,"homer.csv",quote = F,row.names = F,col.names = T)

### Sum
temp_tss <- read.table("~/zjw/20220307_DMmouse/SCAFE_outs/cm.aggregate/annotate/DM/bed/DM.CRE.annot.bed.gz",sep = "\t")
temp_tss <- temp_tss[,c(1,7,8,4,5,6)]
temp_tss <- temp_tss[temp_tss$V5>10,]

temp_tss_region <- bedtoolsr::bt.intersect(a = temp_tss, b = region, wa = T, wb = T, s= T)



for(i in 1:nrow(temp_tss)) {
  a <- temp_tss_region[temp_tss_region[,1]==temp_tss[i,1] & temp_tss_region[,3]==temp_tss[i,3] & temp_tss_region[,6]==temp_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000","Enhancer"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  temp_tss[i,7] <- b
}
temp_tss[is.na(temp_tss[,7]),7] <- "intergenic"
temp_tss[temp_tss[,7]=="TES1000",7] <- "intergenic"
temp_tss[temp_tss[,7]=="last_exon",7] <- "other_exon"
temp_tss[temp_tss[,7]=="last_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="one_exon",7] <- "first_exon"
temp_tss[temp_tss[,7]=="one_intron",7] <- "first_intron"
temp_tss[temp_tss[,7]=="first_exon",7] <- "first exon"
temp_tss[temp_tss[,7]=="first_intron",7] <- "first intron"
temp_tss[temp_tss[,7]=="other_exon",7] <- "other exon"
temp_tss[temp_tss[,7]=="other_intron",7] <- "other intron"
temp_tss[temp_tss[,7]=="5UTR",7] <- "5'UTR"
temp_tss[temp_tss[,7]=="3UTR",7] <- "3'UTR"
temp_tss[temp_tss[,7]=="TSS1000",7] <- "TSS±1000"
temp_tss_df <- as.data.frame(table(temp_tss[,7]))
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c("TSS±1000", "Enhancer", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]
temp_tss_df[,1] <- paste(temp_tss_df[,1]," (",round(temp_tss_df[,2]/sum(temp_tss_df[,2])*100,2),"%)",sep = "")
temp_tss_df$Var1 <- factor(temp_tss_df$Var1, levels=c(grep("TSS",temp_tss_df[,1],value = T),
                                                      grep("Enhancer",temp_tss_df[,1],value = T),
                                                      grep("5'UTR",temp_tss_df[,1],value = T),
                                                      grep("first exon",temp_tss_df[,1],value = T),
                                                      grep("first intron",temp_tss_df[,1],value = T),
                                                      grep("other exon",temp_tss_df[,1],value = T),
                                                      grep("other intron",temp_tss_df[,1],value = T),
                                                      grep("3'UTR",temp_tss_df[,1],value = T),
                                                      grep("intergenic",temp_tss_df[,1],value = T)), ordered=TRUE)
temp_tss_df <- temp_tss_df[order(temp_tss_df[,1]),]


temp_tss_df <- temp_tss_df[grep("intergenic",temp_tss_df$Var1,invert = T),]

p7 <- ggplot(temp_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL,title = "Sum") +
  scale_x_continuous(breaks = NULL)




p1+p2+p3+p4+p5+p6