# Sys.setenv(PKG_CONFIG_PATH = "/GPUFS/gzzoc_yjhu_3/anaconda3/lib/pkgconfig")
# Sys.setenv("PROJ_LIB"="/GPUFS/gzzoc_yjhu_3/anaconda3/include")

options(stringsAsFactors = FALSE)
# options(scipen = 100)
options(bedtools.path = "~/anaconda3/bin/")



.libPaths()

# library(MASS)
# library(BuenColors)
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
library(proj4,lib.loc = "/GPUFS/gzzoc_yjhu_3/anaconda3/envs/Benchmarking/lib/R/library")
library(ggalt,lib.loc = "/GPUFS/gzzoc_yjhu_3/anaconda3/envs/Benchmarking/lib/R/library")

gencode <- read.table("~/zjw/annotation/gencode_mm9_all_gene_all_transcript.bed",header = FALSE,sep = "\t")
gencode_plus <- gencode[gencode[,6]=="+",]
gencode_minus <- gencode[gencode[,6]=="-",]


tss_refer <- rbind(data.frame(gencode_plus[,c(1,2)],V3=gencode_plus[,2]+1,gencode_plus[,c(4,5,6)]),
                   data.frame(V1=gencode_minus[,1],V2=gencode_minus[,3]-1,gencode_minus[c(3,4,5,6)]))



## DM1

DM1_tss_gene <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/SRR12443300/annotate/SRR12443300/bed/SRR12443300.CRE.coord.bed.gz",sep = "\t")
DM1_tss_gene <- DM1_tss_gene[,c(1,7,8,4,5,6)]

DM1_tss_gene_ref <- c()
for(i in 1:nrow(DM1_tss_gene)) {
  b <- tss_refer[tss_refer[,1]==DM1_tss_gene[i,1] & tss_refer[,6]==DM1_tss_gene[i,6],]
  if(DM1_tss_gene[i,6]=="+") {
    DM1_tss_gene_ref <- c(DM1_tss_gene_ref,(DM1_tss_gene[i,3]-b[,3])[order(abs(DM1_tss_gene[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    DM1_tss_gene_ref <- c(DM1_tss_gene_ref,(b[,3]-DM1_tss_gene[i,3])[order(abs(b[,3]-DM1_tss_gene[i,3]),decreasing = FALSE)][1])
  }
}


DM1_tss_df <- data.frame(V1=as.data.frame(table(DM1_tss_gene_ref))[,1],V2=as.data.frame(table(DM1_tss_gene_ref))[,2],V3=rep("test",nrow(as.data.frame(table(DM1_tss_gene_ref)))))

DM1_tss_df[,1] <- as.numeric(as.character(DM1_tss_df[,1]))
DM1_tss_df <- DM1_tss_df[DM1_tss_df[,1]<301 & DM1_tss_df[,1]>-301,] 
# DM1_tss_df <- rbind(DM1_tss_df,
#                     data.frame(V1=setdiff(-1000:1000,DM1_tss_df$V1),V2=0,V3="test")
#                     )



## DM2

DM2_tss_gene <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/SRR12443301/annotate/SRR12443301/bed/SRR12443301.CRE.coord.bed.gz",sep = "\t")
DM2_tss_gene <- DM2_tss_gene[,c(1,7,8,4,5,6)]

DM2_tss_gene_ref <- c()
for(i in 1:nrow(DM2_tss_gene)) {
  b <- tss_refer[tss_refer[,1]==DM2_tss_gene[i,1] & tss_refer[,6]==DM2_tss_gene[i,6],]
  if(DM2_tss_gene[i,6]=="+") {
    DM2_tss_gene_ref <- c(DM2_tss_gene_ref,(DM2_tss_gene[i,3]-b[,3])[order(abs(DM2_tss_gene[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    DM2_tss_gene_ref <- c(DM2_tss_gene_ref,(b[,3]-DM2_tss_gene[i,3])[order(abs(b[,3]-DM2_tss_gene[i,3]),decreasing = FALSE)][1])
  }
}


DM2_tss_df <- data.frame(V1=as.data.frame(table(DM2_tss_gene_ref))[,1],V2=as.data.frame(table(DM2_tss_gene_ref))[,2],V3=rep("test",nrow(as.data.frame(table(DM2_tss_gene_ref)))))

DM2_tss_df[,1] <- as.numeric(as.character(DM2_tss_df[,1]))
DM2_tss_df <- DM2_tss_df[DM2_tss_df[,1]<301 & DM2_tss_df[,1]>-301,] 




## DM3

DM3_tss_gene <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/SRR12443302/annotate/SRR12443302/bed/SRR12443302.CRE.coord.bed.gz",sep = "\t")
DM3_tss_gene <- DM3_tss_gene[,c(1,7,8,4,5,6)]

DM3_tss_gene_ref <- c()
for(i in 1:nrow(DM3_tss_gene)) {
  b <- tss_refer[tss_refer[,1]==DM3_tss_gene[i,1] & tss_refer[,6]==DM3_tss_gene[i,6],]
  if(DM3_tss_gene[i,6]=="+") {
    DM3_tss_gene_ref <- c(DM3_tss_gene_ref,(DM3_tss_gene[i,3]-b[,3])[order(abs(DM3_tss_gene[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    DM3_tss_gene_ref <- c(DM3_tss_gene_ref,(b[,3]-DM3_tss_gene[i,3])[order(abs(b[,3]-DM3_tss_gene[i,3]),decreasing = FALSE)][1])
  }
}


DM3_tss_df <- data.frame(V1=as.data.frame(table(DM3_tss_gene_ref))[,1],V2=as.data.frame(table(DM3_tss_gene_ref))[,2],V3=rep("test",nrow(as.data.frame(table(DM3_tss_gene_ref)))))

DM3_tss_df[,1] <- as.numeric(as.character(DM3_tss_df[,1]))
DM3_tss_df <- DM3_tss_df[DM3_tss_df[,1]<301 & DM3_tss_df[,1]>-301,] 




## CTRL1

CTRL1_tss_gene <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/SRR12443303/annotate/SRR12443303/bed/SRR12443303.CRE.coord.bed.gz",sep = "\t")
CTRL1_tss_gene <- CTRL1_tss_gene[,c(1,7,8,4,5,6)]

CTRL1_tss_gene_ref <- c()
for(i in 1:nrow(CTRL1_tss_gene)) {
  b <- tss_refer[tss_refer[,1]==CTRL1_tss_gene[i,1] & tss_refer[,6]==CTRL1_tss_gene[i,6],]
  if(CTRL1_tss_gene[i,6]=="+") {
    CTRL1_tss_gene_ref <- c(CTRL1_tss_gene_ref,(CTRL1_tss_gene[i,3]-b[,3])[order(abs(CTRL1_tss_gene[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    CTRL1_tss_gene_ref <- c(CTRL1_tss_gene_ref,(b[,3]-CTRL1_tss_gene[i,3])[order(abs(b[,3]-CTRL1_tss_gene[i,3]),decreasing = FALSE)][1])
  }
}


CTRL1_tss_df <- data.frame(V1=as.data.frame(table(CTRL1_tss_gene_ref))[,1],V2=as.data.frame(table(CTRL1_tss_gene_ref))[,2],V3=rep("test",nrow(as.data.frame(table(CTRL1_tss_gene_ref)))))

CTRL1_tss_df[,1] <- as.numeric(as.character(CTRL1_tss_df[,1]))
CTRL1_tss_df <- CTRL1_tss_df[CTRL1_tss_df[,1]<301 & CTRL1_tss_df[,1]>-301,] 




## CTRL2

CTRL2_tss_gene <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/SRR12443304/annotate/SRR12443304/bed/SRR12443304.CRE.coord.bed.gz",sep = "\t")
CTRL2_tss_gene <- CTRL2_tss_gene[,c(1,7,8,4,5,6)]

CTRL2_tss_gene_ref <- c()
for(i in 1:nrow(CTRL2_tss_gene)) {
  b <- tss_refer[tss_refer[,1]==CTRL2_tss_gene[i,1] & tss_refer[,6]==CTRL2_tss_gene[i,6],]
  if(CTRL2_tss_gene[i,6]=="+") {
    CTRL2_tss_gene_ref <- c(CTRL2_tss_gene_ref,(CTRL2_tss_gene[i,3]-b[,3])[order(abs(CTRL2_tss_gene[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    CTRL2_tss_gene_ref <- c(CTRL2_tss_gene_ref,(b[,3]-CTRL2_tss_gene[i,3])[order(abs(b[,3]-CTRL2_tss_gene[i,3]),decreasing = FALSE)][1])
  }
}


CTRL2_tss_df <- data.frame(V1=as.data.frame(table(CTRL2_tss_gene_ref))[,1],V2=as.data.frame(table(CTRL2_tss_gene_ref))[,2],V3=rep("test",nrow(as.data.frame(table(CTRL2_tss_gene_ref)))))

CTRL2_tss_df[,1] <- as.numeric(as.character(CTRL2_tss_df[,1]))
CTRL2_tss_df <- CTRL2_tss_df[CTRL2_tss_df[,1]<301 & CTRL2_tss_df[,1]>-301,] 



## CTRL3

CTRL3_tss_gene <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/SRR12443305/annotate/SRR12443305/bed/SRR12443305.CRE.coord.bed.gz",sep = "\t")
CTRL3_tss_gene <- CTRL3_tss_gene[,c(1,7,8,4,5,6)]

CTRL3_tss_gene_ref <- c()
for(i in 1:nrow(CTRL3_tss_gene)) {
  b <- tss_refer[tss_refer[,1]==CTRL3_tss_gene[i,1] & tss_refer[,6]==CTRL3_tss_gene[i,6],]
  if(CTRL3_tss_gene[i,6]=="+") {
    CTRL3_tss_gene_ref <- c(CTRL3_tss_gene_ref,(CTRL3_tss_gene[i,3]-b[,3])[order(abs(CTRL3_tss_gene[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    CTRL3_tss_gene_ref <- c(CTRL3_tss_gene_ref,(b[,3]-CTRL3_tss_gene[i,3])[order(abs(b[,3]-CTRL3_tss_gene[i,3]),decreasing = FALSE)][1])
  }
}


CTRL3_tss_df <- data.frame(V1=as.data.frame(table(CTRL3_tss_gene_ref))[,1],V2=as.data.frame(table(CTRL3_tss_gene_ref))[,2],V3=rep("test",nrow(as.data.frame(table(CTRL3_tss_gene_ref)))))

CTRL3_tss_df[,1] <- as.numeric(as.character(CTRL3_tss_df[,1]))
CTRL3_tss_df <- CTRL3_tss_df[CTRL3_tss_df[,1]<301 & CTRL3_tss_df[,1]>-301,] 




## Sum

Sum_tss_gene <- read.table("~/zjw/20220307DMmouse/SCAFE_outs/cm.aggregate/annotate/DM/bed/DM.CRE.annot.bed.gz",sep = "\t")
Sum_tss_gene <- Sum_tss_gene[,c(1,7,8,4,5,6)]

Sum_tss_gene_ref <- c()
for(i in 1:nrow(Sum_tss_gene)) {
  b <- tss_refer[tss_refer[,1]==Sum_tss_gene[i,1] & tss_refer[,6]==Sum_tss_gene[i,6],]
  if(Sum_tss_gene[i,6]=="+") {
    Sum_tss_gene_ref <- c(Sum_tss_gene_ref,(Sum_tss_gene[i,3]-b[,3])[order(abs(Sum_tss_gene[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    Sum_tss_gene_ref <- c(Sum_tss_gene_ref,(b[,3]-Sum_tss_gene[i,3])[order(abs(b[,3]-Sum_tss_gene[i,3]),decreasing = FALSE)][1])
  }
}


Sum_tss_df <- data.frame(V1=as.data.frame(table(Sum_tss_gene_ref))[,1],V2=as.data.frame(table(Sum_tss_gene_ref))[,2],V3=rep("test",nrow(as.data.frame(table(Sum_tss_gene_ref)))))

Sum_tss_df[,1] <- as.numeric(as.character(Sum_tss_df[,1]))
Sum_tss_df <- Sum_tss_df[Sum_tss_df[,1]<301 & Sum_tss_df[,1]>-301,] 





tss_df <- rbind(data.frame(DM1_tss_df[,c(1,2)],V3="DM1"),
                data.frame(DM2_tss_df[,c(1,2)],V3="DM2"),
                data.frame(DM3_tss_df[,c(1,2)],V3="DM3"),
                data.frame(CTRL1_tss_df[,c(1,2)],V3="CTRL1"),
                data.frame(CTRL2_tss_df[,c(1,2)],V3="CTRL2"),
                data.frame(CTRL3_tss_df[,c(1,2)],V3="CTRL3"))

#tss_df[tss_df[,3]=="tss",2] <- tss_df[tss_df[,3]=="tss",2]/13249
#tss_df[tss_df[,3]=="tes",2] <- tss_df[tss_df[,3]=="tes",2]/20376






p1 <- ggplot(DM1_tss_df,aes(V1,V2))+ 
  scale_x_continuous(breaks=seq(-300,300,100))+
  # scale_colour_manual(values=c("red", "blue"),
  #                     breaks=c("tes","tss"),
  #                     labels=c("tes","tss"))+ 
  geom_xspline(data = DM1_tss_df,spline_shape = 1,size=1)+
  theme_minimal()+theme(text=element_text(size=17))+
  labs(x="Distance to the annotated TSS",y="Number of CRE",col="cell type",title = "DM1")+
  #guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))




p2 <- ggplot(DM2_tss_df,aes(V1,V2))+ 
  scale_x_continuous(breaks=seq(-300,300,100))+
  # scale_colour_manual(values=c("red", "blue"),
  #                     breaks=c("tes","tss"),
  #                     labels=c("tes","tss"))+ 
  geom_xspline(data = DM2_tss_df,spline_shape = 1,size=1)+
  theme_minimal()+theme(text=element_text(size=17))+
  labs(x="Distance to the annotated TSS",y="Number of CRE",col="cell type",title = "DM2")+
  #guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))




p3 <- ggplot(DM3_tss_df,aes(V1,V2))+ 
  scale_x_continuous(breaks=seq(-300,300,100))+
  # scale_colour_manual(values=c("red", "blue"),
  #                     breaks=c("tes","tss"),
  #                     labels=c("tes","tss"))+ 
  geom_xspline(data = DM3_tss_df,spline_shape = 1,size=1)+
  theme_minimal()+theme(text=element_text(size=17))+
  labs(x="Distance to the annotated TSS",y="Number of CRE",col="cell type",title = "DM3")+
  #guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))




p4 <- ggplot(CTRL1_tss_df,aes(V1,V2))+ 
  scale_x_continuous(breaks=seq(-300,300,100))+
  # scale_colour_manual(values=c("red", "blue"),
  #                     breaks=c("tes","tss"),
  #                     labels=c("tes","tss"))+ 
  geom_xspline(data = CTRL1_tss_df,spline_shape = 1,size=1)+
  theme_minimal()+theme(text=element_text(size=17))+
  labs(x="Distance to the annotated TSS",y="Number of CRE",col="cell type",,title = "CTRL1")+
  #guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))




p5 <- ggplot(CTRL2_tss_df,aes(V1,V2))+ 
  scale_x_continuous(breaks=seq(-300,300,100))+
  # scale_colour_manual(values=c("red", "blue"),
  #                     breaks=c("tes","tss"),
  #                     labels=c("tes","tss"))+ 
  geom_xspline(data = CTRL2_tss_df,spline_shape = 1,size=1)+
  theme_minimal()+theme(text=element_text(size=17))+
  labs(x="Distance to the annotated TSS",y="Number of CRE",col="cell type",title = "CTRL2")+
  #guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))




p6 <- ggplot(CTRL3_tss_df,aes(V1,V2))+ 
  scale_x_continuous(breaks=seq(-300,300,100))+
  # scale_colour_manual(values=c("red", "blue"),
  #                     breaks=c("tes","tss"),
  #                     labels=c("tes","tss"))+ 
  geom_xspline(data = CTRL3_tss_df,spline_shape = 1,size=1)+
  theme_minimal()+theme(text=element_text(size=17))+
  labs(x="Distance to the annotated TSS",y="Number of CRE",col="cell type",title = "CTRL3")+
  #guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))



p7 <- ggplot(Sum_tss_df,aes(V1,V2))+ 
  scale_x_continuous(breaks=seq(-250,250,50))+
  # scale_colour_manual(values=c("red", "blue"),
  #                     breaks=c("tes","tss"),
  #                     labels=c("tes","tss"))+ 
  geom_xspline(data = Sum_tss_df,spline_shape = 1,size=1)+
  theme_minimal()+theme(text=element_text(size=17))+
  labs(x="Distance to the annotated TSS",y="Number of CRE",col="cell type",title = "Sum")+
  #guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))

p1+p2+p3+p4+p5+p6
