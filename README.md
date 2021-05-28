Supporting Material for Mathias2021
================
Carolina Mathias, João Carlos Degraf Muzzi, Bruna Borba Antunes, Daniela F. Gradia, Mauro A. A. Castro and Jaqueline C. de Oliveira. <br>

04 May 2021.

Context
----
The immune subtype classification presented by Thorsson and colleagues (2018) using mRNA represented a mark in immuno-oncology. In Mathias et al., 2021, we add the complexity of lncRNAs to this characterization. We show that lncRNAs may add valuable information on the immune and molecular subtypes characterization in breast cancer and we present a signature composed by 10 to 11 lncRNAs for each molecular subtype associated with the immune subtypes and immune features. These signatures highlight little-known and some better-studied lncRNAs that may help guide future studies on possible biomarkers for therapeutic design and patients’ selection for clinical trials. This script reproduces all results published in Mathias et al. (2021) and serves as complementary material.

In Preprocessing.Rmd file, we show how to obtain and preprocess the gene expression matrix and Thorsson et al. (2018) master table used in this study. At the end of the section the objects are saved as RData and can be imported for further analysis. Here we use the preprocessed data and show the steps for data analysis and how to obtain the results and plots shown in Mathias et al. (2021). 

Loading packages and pre-process data
----
```r
# Load Libraries
library(tidyverse)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(cowplot)
library(ggpubr)
library(fgsea)
library(msigdbr)

```

```r
# Import preprocessed Thorsson et al (2018) data set
repo_link <- "https://github.com/sysbiolab/Sup_Material_Mathias2021/"
name.file <- "blob/main/thorsson_BRCA.RData"
download.file(
   url= paste0(repo_link,name.file),
   destfile = "thorsson_BRCA.RData")
 
load("thorsson_BRCA.RData")
 
# Loading gene expression matrix saved in Preprocess_BRCA.Rmd
load("gexp_BRCA.RData")

longos.names <- gexp$lnc.Names
gene.names <- gexp$gene.names
rna <- gexp$mRNA.Tumor
rna_healthy <-gexp$mRNA.Healthy
rna_lnc <- gexp$lncRNA.Tumor
rna_lnc_healthy <- gexp$lncRNA.Healthy

rm(gexp)
```



Code snippet for Sup. figure 1
----
```r
### Heatmap for Immuno modulator genes
#### Immune Modulator genes described by Thorsson et al 2018
#### C10orf14 not found in gene expression matrix
co_stimulator <- sort(c("CD80","CD28","ICOSLG"))
co_inhibitor <- sort(c("PDCD1LG2","CD274","VTCN1","SLAMF7","BTN3A2","BTN3A1","CD276")) 
ligand <- sort(c("TNFSF9","TNF","TNFSF4","IL1B","CXCL9","CXCL10","CCL5","VEGFB","CX3CL1",
                 "TGFB1",
                 "VEGFA","CD70","CD40LG","IL10", "IFNG", "IL1A","IL12A", "IFNA2",
                 "IFNA1","IL4","IL2","IL13"))
receptor <- sort(c("TNFRSF18","TIGIT","PDCD1","CTLA4","IL2RA","TNFRSF4","CD27",
                   "LAG3","TNFRSF9","ICOS","BTLA","KIR2DL3","KIR2DL1",
                   "TNFRSF14","EDNRB","CD40","ADORA2A","TLR4","HAVCR2")) 
cell_adhesion <- sort(c("ITGB2","ICAM1","SELP"))
antigen_presentation <-
  sort(c("HLA-DRB5","HLA-DQA1","HLA-DQB1","MICA","MICB","HLA-DQA2","HLA-DQB2","HLA-B",
         "HLA-A","HLA-C","HLA-DRA","HLA-DRB1","HLA-DPB1","HLA-DPA1"))
other <- sort(c("IDO1","GZMA","PRF1","ARG1","HMGB1","ENTPD1"))

IM <- c(co_stimulator,co_inhibitor,ligand,receptor,cell_adhesion,antigen_presentation,other)

rm(antigen_presentation,cell_adhesion,co_inhibitor,
   co_stimulator,ligand,other,receptor)

#### Extracting IM genes from gene expression matrix
dat_imun <- gene.names[gene.names$`HGNC symbol` %in% IM,]
dat_imun <- rna[rownames(rna) %in% dat_imun$`Gene stable ID`,]

# preparing data.frame for heatmap
rows <- gene.names[gene.names$`Gene stable ID` %in% rownames(dat_imun),]
all(rows$`Gene stable ID` %in% rownames(dat_imun)) #TRUE
all(IM %in% rows$`HGNC symbol`) #TRUE
rownames(rows) <- rows$`HGNC symbol`
rows <- rows[IM,]
dat_imun <- dat_imun[rows$`Gene stable ID`,]
rownames(dat_imun) <- rows$`HGNC symbol`

unique(factor(thorsson$TCGA.Subtype))
thorsson <- thorsson[thorsson$TCGA.Subtype!="NA",]
dat_imun <- dat_imun[,thorsson$TCGA.Participant.Barcode]
identical(thorsson$TCGA.Participant.Barcode, colnames(dat_imun)) #TRUE

# Column-wise Z-score and setting maximum and minimum to +2 and -2 
dat_imun <- data.frame(t(dat_imun))
dat_imun <- data.frame(scale(dat_imun))
range(dat_imun)
summary(dat_imun)
dat_imun[dat_imun>2]<-2
dat_imun[dat_imun<(-2)]<-(-2)
range(dat_imun)

#############
# Heatmap row split
split <- as.vector(c(rep("Co-\nStimulator",3),
                     rep("Co-\ninhibitor",7),
                     rep("Ligand",22),
                     rep("Receptor",19),
                     rep("Cell\nadhesion",3),
                     rep("Antigen\npresentation",14),
                     rep("Other",6) ))

split <- factor(split,
                levels = c("Co-\nStimulator","Co-\ninhibitor","Ligand","Receptor",
                           "Cell\nadhesion", "Antigen\npresentation","Other"))

### generating top annotation data.frame
annot_hm <- thorsson[,c("TCGA.Participant.Barcode", "TCGA.Subtype", "OS", "Immune.Subtype",
                        "Leukocyte.Fraction",
                        "Lymphocyte.Infiltration.Signature.Score","Macrophage.Regulation",
                        "Wound.Healing","IFN.gamma.Response", "TGF.beta.Response")]
rownames(annot_hm) <- annot_hm$TCGA.Participant.Barcode
colnames(annot_hm)[c(3,6)] <- c("Vital.Status","Lymphocyte.Inf.Score")
annot_hm$Vital.Status <- ifelse(annot_hm$Vital.Status==1, "Dead", "Alive")
annot_hm <- annot_hm[,-1]
annot_hm[,1:3] <- lapply(annot_hm[,1:3], as.factor)
annot_hm[,4:9] <- lapply(annot_hm[,4:9], scales::rescale)
annot_hm$TCGA.Subtype <- factor(annot_hm$TCGA.Subtype,
                                levels=c("BRCA.Basal", "BRCA.Normal",
                                         "BRCA.Her2", "BRCA.LumA", "BRCA.LumB" ))

identical(rownames(annot_hm), rownames(dat_imun)) #TRUE 

#Plotting Heatmap
h1<- Heatmap(matrix = t(dat_imun),
             name = "Log2 Gene\nexpression\nz-score",
             show_column_names = F, 
             show_column_dend = F,
             column_names_side="bottom",
             column_names_gp= gpar(col="white", fontsize=1),
             column_title_gp = gpar(fontsize=8), 
             column_title = c("Basal\n(n=169)", "Normal\n(n=136)",
                              "LumA\n(n=499)", "LumB\n(n=184)","Her2\n(n=72)"),
             column_split = factor(annot_hm$TCGA.Subtype,
                                   levels=c("BRCA.Basal", "BRCA.Normal",
                                            "BRCA.LumA", "BRCA.LumB","BRCA.Her2" )),
             cluster_columns = T, 
             cluster_column_slices = F,
             cluster_rows = F, 
             row_split = split, 
             row_title_rot = 0, 
             row_title_gp = gpar(fontsize=8),
             row_names_gp = gpar(fontsize=7),
             row_title_side ="left",
             heatmap_legend_param = list(
               labels_gp=gpar(fontsize=8),
               at=c(-2,0,2),
               title_gp=gpar(fontsize=8, fontface="bold"),
               direction = "horizontal",
               border=T,
               legend_width= unit(1.7,"cm")),
             top_annotation = columnAnnotation(
               df=annot_hm[,-c(2,5:9)],
               simple_anno_size= unit(3, "mm"),
               height=unit(1,"cm"),
               annotation_name_gp = gpar(fontsize=8),
               annotation_name_side = "right",
               na_col="white", 
               annotation_legend_param = list(
                 TCGA.Subtype=list(
                   labels_gp=gpar(fontsize=8),
                   title_gp=gpar(fontsize=8, fontface="bold"),
                   nrow=3, 
                   title="TCGA Subtype"),
                 Immune.Subtype=list(
                   labels_gp=gpar(fontsize=8),
                   title_gp=gpar(fontsize=8, fontface="bold"),
                   nrow=3, 
                   title="Immune\nSubtype"),
                 Leukocyte.Fraction=list(
                   title="Leukocyte\nFraction",
                   border=T, legend_width= unit(1.7,"cm"),
                   labels_gp=gpar(fontsize=8),
                   title_gp=gpar(fontsize=8, fontface="bold"),
                   direction="horizontal", at=c(0,0.5,1))),
               col=list(
                 Immune.Subtype=c("C1"="red","C2"="yellow","C3"="green3",
                                  "C4"="cyan","C5"="blue","C6"="pink"),
                 TCGA.Subtype=c("BRCA.Basal"="sienna1", "BRCA.Normal"="purple3",
                                "BRCA.Her2" ="hotpink", "BRCA.LumA"="cornflowerblue",
                                "BRCA.LumB"="springgreen2"),
                 Leukocyte.Fraction=circlize::colorRamp2(c(0,0.5), c("white","darkgreen"))
               ))    )
g<- grid.grabExpr(
  draw(h1,
       heatmap_legend_side="bottom",
       annotation_legend_side="bottom",
       merge_legends=T),
  height = 9,
  width = 7)
plot_grid(g)

rm(dat_imun, g, h1, rna_IM,rows, genes, IM, patients,split)

```

Signal to Noise Ratio Calculation
----
```r
#####################################################################################
# Filtering lncRNA matrix for 75%+ expression
#####################################################################################
rna <- rna[,thorsson$TCGA.Participant.Barcode]
rna_lnc <- rna_lnc[,thorsson$TCGA.Participant.Barcode]
rows <- rowSums(rna_lnc, na.rm = T)
rna_lnc <- rna_lnc[rowSums(rna_lnc, na.rm = T) > quantile(rows, probs = 0.25),]
# lncRNAs with sum greater than 8.04 in 1,061 patients 
#####################################################################################
# Signal to Noise ratio for BC subtypes
#####################################################################################
subtype<- c("BRCA.Basal"  ,"BRCA.LumA"  , "BRCA.LumB"  , "BRCA.Normal", "BRCA.Her2"  )
mean_na <- function(a){mean(a,na.rm=T)}
sd_na <- function(a){sd(a,na.rm=T)}
stnr <- NULL
for( i in 1:length(subtype)){
  gexp <- rna_lnc
  print(paste0(subtype[i]))
  idx <- thorsson$TCGA.Participant.Barcode[thorsson$TCGA.Subtype==subtype[i]]
  gexp1 <- (apply(gexp[,names(gexp) %in% idx],1,mean_na)-apply(gexp[,!(names(gexp) %in% idx)],1,mean_na)) / 
        (apply(gexp[,names(gexp) %in% idx],1,sd_na)+apply(gexp[,!(names(gexp) %in% idx)],1,sd_na)  )
  gexp1 <- c(subtype[i],gexp1)
  stnr <- rbind(stnr,gexp1, deparse.level = 0)
}


stnr <- data.frame(stnr)
rownames(stnr) <- stnr[,1]
rows <- rownames(stnr)
stnr <- stnr[,-1]
stnr <- data.frame(lapply(stnr,as.numeric))
rownames(stnr) <- rows

# Supplementary table with SNR results and quantile distribution
sup.tab <- as.data.frame(t(stnr))
lnc.names <- longos.names
lnc.names$`Gene stable ID`[duplicated(lnc.names$`Gene stable ID`)] 
#[1] "ENSG00000230417" "ENSG00000254876"
lnc.names$`Gene stable ID`[duplicated(lnc.names$`Gene stable ID`)] <-
  paste0(lnc.names$`Gene stable ID`[duplicated(lnc.names$`Gene stable ID`)],"_")
lnc.names$`Gene stable ID`[duplicated(lnc.names$`Gene stable ID`)] #0
rownames(lnc.names) <- lnc.names$`Gene stable ID`
sup.tab <- cbind(lnc.names[rownames(sup.tab),],sup.tab)
sup.tab <- sup.tab[,-c(2:4)]
sup.tab2 <- sup.tab
sup.tab2[,3:7] <- (apply(sup.tab2[,3:7], 2, rank))/nrow(sup.tab2)
sup.tab <- cbind(sup.tab[,1:3],
                 BRCA.Basal.Quantile=sup.tab2$BRCA.Basal,
                 BRCA.LumA=sup.tab$BRCA.LumA,
                 BRCA.LumA.Quantile=sup.tab2$BRCA.LumA,
                 BRCA.LumB=sup.tab$BRCA.LumB,
                 BRCA.LumB.Quantile=sup.tab2$BRCA.LumB,
                 BRCA.Normal=sup.tab$BRCA.Normal,
                 BRCA.Normal.Quantile=sup.tab2$BRCA.Normal,
                 BRCA.Her2=sup.tab$BRCA.Her2,
                 BRCA.Her2.Quantile=sup.tab2$BRCA.Her2)
sup.tab <- sup.tab[order(sup.tab$BRCA.Basal.Quantile,
                         decreasing = T),]
writexl::write_xlsx(sup.tab,path="BRCA.1STNR.xlsx")
# Define quantile (up 95%)
quantil <- 0.95 # For generating lncRNAs list as presented in Sup. Figure 2, alter here for 0.90
# Separating and saving lncRNAs for each molecular subtype  
basal <- t(stnr[1,])
basal <- names(basal[basal > quantile(basal,probs = quantil),])
write_tsv(as.data.frame(basal), file="BRCA.Basal_1snr.tsv", col_names = F)
lumA <- t(stnr[2,])
lumA <- names(lumA[lumA > quantile(lumA,probs = quantil),])
write_tsv(as.data.frame(lumA), file="BRCA.LumA_1snr.tsv", col_names = F)
lumB <- t(stnr[3,])
lumB <- names(lumB[lumB > quantile(lumB,probs = quantil),])
write_tsv(as.data.frame(lumB), file="BRCA.LumB_1snr.tsv", col_names = F)
normal <- t(stnr[4,])
normal <- names(normal[normal > quantile(normal,probs = quantil),])
write_tsv(as.data.frame(normal), file="BRCA.Normal_1snr.tsv", col_names = F)
her2 <- t(stnr[5,])
her2 <- names(her2[her2 > quantile(her2,probs = quantil),])
write_tsv(as.data.frame(her2), file="BRCA.Her2_1snr.tsv", col_names = F)
```


```r
#####################################################################################
# Signal to Noise ratio for Immune Subtypes
#####################################################################################
subtype<- c("BRCA.Basal"  ,"BRCA.LumA"  , "BRCA.LumB"  , "BRCA.Normal", "BRCA.Her2"  )
immune_subtype <- c("C1","C2","C3","C4","C5","C6")
lnc.names <- lnc.names[,c("Gene stable ID","symbol.ens")]
colnames(lnc.names) <- c("Gene stable ID", "Gene name" )
rownames(lnc.names) <- lnc.names$`Gene stable ID`

for( i in 1:length(subtype)){
  thorsson_sub <- thorsson[thorsson$TCGA.Subtype==subtype[i],]
  rna_lnc_sub <- rna_lnc[,(thorsson_sub$TCGA.Participant.Barcode)]
  if(subtype[i]=="BRCA.Basal"){lnc_sub <- basal}
  if(subtype[i]=="BRCA.LumA"){lnc_sub <- lumA}
  if(subtype[i]=="BRCA.LumB"){lnc_sub <- lumB}  
  if(subtype[i]=="BRCA.Normal"){lnc_sub <- normal}
  if(subtype[i]=="BRCA.Her2"){lnc_sub <- her2}
  
  gexp <- rna_lnc_sub[lnc_sub, ] 
  stnr <- NULL
  x=1

  while(x <= 6){
    print(paste0(subtype[i]," ", immune_subtype[x]))
    idx <- thorsson_sub$TCGA.Participant.Barcode[thorsson_sub$Immune.Subtype==immune_subtype[x]]
    
    if(sum(names(gexp) %in% idx) >5){
      gexp1 <- (apply(gexp[,names(gexp) %in% idx],1,mean_na)-apply(gexp[,!(names(gexp) %in% idx)],1,mean_na)) / 
        (apply(gexp[,names(gexp) %in% idx],1,sd_na)+apply(gexp[,!(names(gexp) %in% idx)],1,sd_na)  )
      gexp1 <- c(immune_subtype[x],gexp1)
      stnr <- rbind(stnr,gexp1, deparse.level = 0)
    }
    x<-x+1
  }
  
  stnr <- data.frame(stnr)
  rownames(stnr) <- stnr[,1]
  rows <- rownames(stnr)
  stnr <- stnr[,-1]
  stnr <- data.frame(lapply(stnr,as.numeric))
  rownames(stnr) <- rows
  soma <- data.frame(colSums(abs(stnr)))
  colnames(soma) <- "STNR.Sum"
  soma2<- soma
  colnames(soma2) <- paste0(subtype[i],".STNR.Sum")
  soma2$Quantile <- (rank(soma2[,1]))/nrow(soma2)
  soma2$Z.score <- (scale(soma2[,1]))
  mess <- paste0("Order quantile = z-score? ",
               identical(rownames(soma2[order(soma2$Quantile),]),
                         rownames(soma2[order(soma2$Z.score),])))
  print(mess)

  soma2 <- cbind(lnc.names[rownames(soma2),], soma2)
  soma2 <- soma2[order(soma2$Quantile, decreasing = T),]
  writexl::write_xlsx(soma2, path=paste0(subtype[i],"_2SNR_sum.xlsx"))
  # Define quantile 98%
  quantil <- 0.98
  # Selecting lncRNAs for each molecular subtype
  if(subtype[i]=="BRCA.Her2"){
    soma_her2 <- soma
    # select lncRNAs above quantile definition
    lnc_98stn <- quantile(soma[,1], quantil)
    lnc_98stn <- as.data.frame(rownames(soma[soma$STNR.Sum> lnc_98stn ,1, drop=F]))
    lnc.h <- lnc_98stn
    # Save selected lncRNAs
    #write_tsv(lnc_98stn, file= "BRCA.Her2_2stn.tsv", col_names = F)
  }
  if(subtype[i]=="BRCA.Basal"){
    soma_basal <- soma
    # select lncRNAs above quantile definition
    lnc_98stn <- quantile(soma[,1], quantil)
    lnc_98stn <- as.data.frame(rownames(soma[soma$STNR.Sum> lnc_98stn ,1, drop=F]))
    lnc.b <- lnc_98stn
    # Save selected lncRNAs
    #write_tsv(lnc_98stn, file= "BRCA.Basal_2stn.tsv", col_names = F)
  }
  if(subtype[i]=="BRCA.LumA"){
    soma_lumA <- soma
    # select lncRNAs above quantile definition
    lnc_98stn <- quantile(soma[,1], quantil)
    lnc_98stn <- as.data.frame(rownames(soma[soma$STNR.Sum> lnc_98stn ,1, drop=F]))
    lnc.la <- lnc_98stn
    # Save selected lncRNAs
    #write_tsv(lnc_98stn, file= "BRCA.LumA_2stn.tsv", col_names = F)
  }
  if(subtype[i]=="BRCA.LumB"){
    soma_lumB <- soma
    # select lncRNAs above quantile definition
    lnc_98stn <- quantile(soma[,1], quantil)
    lnc_98stn <- as.data.frame(rownames(soma[soma$STNR.Sum> lnc_98stn ,1, drop=F]))
    lnc.lb <- lnc_98stn
    # Save selected lncRNAs
    #write_tsv(lnc_98stn, file= "BRCA.LumB_2stn.tsv", col_names = F)
  }
  if(subtype[i]=="BRCA.Normal"){
    soma_normal <- soma
    # select lncRNAs above quantile definition
    lnc_98stn <- quantile(soma[,1], quantil)
    lnc_98stn <- as.data.frame(rownames(soma[soma$STNR.Sum> lnc_98stn ,1, drop=F]))
    lnc.n <- lnc_98stn
    # Save selected lncRNAs
    #write_tsv(lnc_98stn, file= "BRCA.Normal_2stn.tsv", col_names = F)
  }
}
# All 53 lncRNAs
lnc <- rbind(lnc.b, lnc.n, lnc.la, lnc.lb, lnc.h)

rm(gexp, gexp1,lnc_98stn, rna_lnc_sub,soma,thorsson_sub, i, idx,rows, x)

```

Code snippet for Figure 2B
----
```r
### Plotting histogram to show selection criteria and distribution of SNR 
### in each molecular subtype
#Scale sum (range 0 to 1)
soma_basal$STNR.Sum <- scale(soma_basal$STNR.Sum) 
soma_normal$STNR.Sum <- scale(soma_normal$STNR.Sum) 
soma_lumA$STNR.Sum <- scale(soma_lumA$STNR.Sum) 
soma_lumB$STNR.Sum <- scale(soma_lumB$STNR.Sum)   
soma_her2$STNR.Sum <- scale(soma_her2$STNR.Sum) 

g <- ggplot(soma_basal,aes(x=STNR.Sum))+
  geom_histogram(color="black", fill="sienna1", bins = 50)+ 
  xlab("Scaled \U03A3 |SNR|")+
  ylab("Number of lncRNAs")+
  labs(title="Basal")+ 
  geom_vline(xintercept = quantile(soma_basal[,1], 0.98), linetype=2, color="#800026")+
  annotate(
     geom="text",
     y=25, x=quantile(soma_basal[,1], 0.98)*1.05,
     label="quantile > 0.98", 
     color="#800026", 
     size=3, 
     hjust = 0, 
     fontface="bold")+
  coord_cartesian(ylim=c(0,75), xlim=c(-3,8))
g.basal<- plot_grid(g)

g <- ggplot(soma_normal,aes(x=STNR.Sum))+ 
  geom_histogram(color="black", fill="purple3", bins = 50)+ 
  xlab("Scaled \U03A3 |SNR|")+
  ylab("Number of lncRNAs")+
  labs(title="Normal")+ 
  geom_vline(xintercept = quantile(soma_normal[,1], 0.98), linetype=2, color="#800026")+
  annotate(
     geom="text",
     y=25, quantile(soma_normal[,1], 0.98)*1.05,
     label="quantile > 0.98", 
     color="#800026", 
     size=3, 
     hjust = 0, 
     fontface="bold")+
  coord_cartesian(ylim=c(0,75), xlim=c(-3,8))
g.normal<- plot_grid(g) 

g <- ggplot(soma_lumA,aes(x=STNR.Sum))+
  geom_histogram(color="black", fill="cornflowerblue", bins = 50)+
  xlab("Scaled \U03A3 |SNR|")+
  ylab("Number of lncRNAs")+
  labs(title="LumA")+ 
  geom_vline(xintercept = quantile(soma_lumA[,1], 0.98), linetype=2, color="#800026")+
  annotate(
     geom="text",
     y=25, x=quantile(soma_lumA[,1], 0.98)*1.05,
     label="quantile > 0.98", 
     color="#800026", 
     size=3, 
     hjust = 0, 
     fontface="bold")+
  coord_cartesian(ylim=c(0,75), xlim=c(-3,8))
g.lumA <- plot_grid(g) 

g <- ggplot(soma_lumB,aes(x=STNR.Sum))+
  geom_histogram(color="black", fill="springgreen2", bins = 50)+ 
  xlab("Scaled \U03A3 |SNR|")+
  ylab("Number of lncRNAs")+
  labs(title="LumB")+ 
  geom_vline(xintercept = quantile(soma_lumB[,1], 0.98), linetype=2, color="#800026")+
  annotate(
     geom="text",
     y=25, quantile(soma_lumB[,1], 0.98)*1.05,
     label="quantile > 0.98", 
     color="#800026", 
     size=3, 
     hjust = 0, 
     fontface="bold")+
  coord_cartesian(ylim=c(0,75), xlim=c(-3,8))
g.lumB <- plot_grid(g) 

g <- ggplot(soma_her2,aes(x=STNR.Sum))+ 
  geom_histogram(color="black", fill="hotpink", bins = 50)+ 
  xlab("Scaled \U03A3 |SNR|")+
  ylab("Number of lncRNAs")+
  labs(title="Her2")+ 
  geom_vline(xintercept = quantile(soma_her2[,1], 0.98), linetype=2, color="#800026")+
  annotate(
     geom="text",
     y=25, x=quantile(soma_her2[,1], 0.98)*1.05,
     label="quantile > 0.98", 
     color="#800026", 
     size=3, 
     hjust = 0, 
     fontface="bold")+
  coord_cartesian(ylim=c(0,75), xlim=c(-3,8))
g.her2 <- plot_grid(g) 

g<- plot_grid(g.basal,g.normal,g.lumA, g.lumB, g.her2, ncol=1) #pdf 4x8in
g
ggsave(plot=g, device = "svg", width = 4, height = 8, filename = "histogramas_98.svg", dpi=300)

#rm(g, g.basal, g.her2, g.lumA, g.lumB, g.normal, soma_basal, soma_her2, soma_lumA,
#   soma_lumB, soma_normal)
```


Specific lncRNAs for each molecular subtype
----
```r
#####################################################################################
# Specific lncRNAs for each molecular lncRNAs
#####################################################################################
# Excluding shared lncRNAs between the molecular subtypes
for(i in 1:length(subtype)){
  if(subtype[i]=="BRCA.Normal"){
    print(subtype[i])
    lnc_nao <- rbind(lnc.la, lnc.lb, lnc.h, lnc.b)
    spec <- lnc.n[!lnc.n[,1] %in% lnc_nao[,1],]
    spec_normal <- spec
    #write_tsv(spec, "spec_normal.tsv")
  }
  
  if(subtype[i]=="BRCA.Her2"){
    print(subtype[i])
    lnc_nao <- rbind(lnc.b, lnc.la, lnc.lb, lnc.n)
    spec <- lnc.h[!lnc.h[,1] %in% lnc_nao[,1],]
    spec_her2 <- spec
    #write_tsv(spec, "spec_her2.tsv")
  }
  
  if(subtype[i]=="BRCA.Basal"){
    print(subtype[i])
    lnc_nao <- rbind(lnc.h, lnc.la, lnc.lb, lnc.n)
    spec <- lnc.b[!lnc.b[,1] %in% lnc_nao[,1],]
    spec_basal <- spec
    #write_tsv(spec, "spec_Basal.tsv")
  }
  
  if(subtype[i]=="BRCA.LumA"){
    print(subtype[i])
    lnc_nao <- rbind(lnc.h, lnc.b, lnc.lb, lnc.n)
    spec <- lnc.la[!lnc.la[,1] %in% lnc_nao[,1],]
    spec_la <- spec
    #write_tsv(spec, "spec_lumA.tsv")
  }
  
  if(subtype[i]=="BRCA.LumB"){
    print(subtype[i])
    lnc_nao <- rbind(lnc.h, lnc.b, lnc.la, lnc.n)
    spec <- lnc.lb[!lnc.lb[,1] %in% lnc_nao[,1],]
    spec_lb <- spec
    #write_tsv(spec, "spec_lumB.tsv")
  }
}  
```
Code snippet for Supplementary Table 4 - Cox univariate results and Leukocyte fraction correlation
----
```r
### Cox univariate analysis
##### For Overall Survival (OS)
for( i in 1:length(subtype)){
  print(subtype[i])
  annot <- thorsson[thorsson$TCGA.Subtype==subtype[i],c(1:5,33:36)]
  rownames(annot) <- annot$TCGA.Participant.Barcode

  if(subtype[i]=="BRCA.Her2"){spec <- spec_her2}
  if(subtype[i]=="BRCA.Basal"){spec <- spec_basal}
  if(subtype[i]=="BRCA.LumA"){spec <- spec_la}
  if(subtype[i]=="BRCA.LumB"){spec <- spec_lb}
  if(subtype[i]=="BRCA.Normal"){spec <- spec_normal}

  gexp <- rna_lnc[spec, rownames(annot)]
  gexp<- (as.data.frame(t(gexp)))
  print(paste0(subtype[i], 
               " Identical= ",
               identical(rownames(gexp), rownames(annot))))
  
  d <- as.vector(colnames(gexp))
  gexp <- cbind(annot,rownames(gexp), gexp)
  # Function for Cox univariate analysis
  coxUnivariateAnalyses <- function(annot, vars4test, time, event, save = TRUE,
                                    fpath = subtype[i]) {
    #-- Rename Survvars
    annot_rn <- annot %>%
      dplyr::rename(time = !! time, event = !! event) %>%
      mutate(across(contains(" mutation"), ~ factor(., levels = c("WT", "Mut"))))
    
    models <- lapply(vars4test, function(varn) {
      message(varn)
      ts_form <- paste0("Surv(time, event) ~ `", varn, "`") %>% as.formula()
      model_res <- coxph(ts_form, data = annot_rn) %>% summary()
    })
    names(models) <- vars4test
    
    var_pvals <- sapply(models, function(model) { model$waldtest[3] })
    
    var_hrs <- lapply(models, function(model) { 
      model$conf.int %>%
        as.data.frame() %>%
        rownames_to_column("Category") %>%
        dplyr::rename(HR = "exp(coef)") %>%
        dplyr::select(-"exp(-coef)")
    }) %>% bind_rows(.id = "Variable")
    
    univar_cox <- data.frame(Variable = names(var_pvals), OS.Cox_pvalue = var_pvals) %>%
      mutate(OS.Cox.p.adj = p.adjust(OS.Cox_pvalue, method = "fdr") ) %>%
      mutate(Variable = str_remove(Variable, ".pvalue"),
             OS.Significant = ifelse(OS.Cox.p.adj < 0.1, "*", "ns")) %>%
      arrange(OS.Cox.p.adj) %>%
      full_join(var_hrs, by = "Variable") %>%
      dplyr::select(Variable, Category, everything())
    colnames(univar_cox)[6:8] <- c("OS.HR","OS.lower.95" ,"OS.upper.95")
    rownames(univar_cox) <- univar_cox$Variable
   
    if(save) openxlsx::write.xlsx(univar_cox, file = paste0(fpath, "_Cox_results.xlsx"))
    
    return(list(cox_models = models, table = univar_cox))
  }
# Running analysis
final <- coxUnivariateAnalyses(gexp,d,time="OS.Time", event = "OS", save = F)
final <- as.data.frame(final$table)
rownames(final) <- final[,1]
final <- final[,-1]
if(subtype[i]=="BRCA.Her2"){lnc_cox_her2 <- final}
if(subtype[i]=="BRCA.Basal"){lnc_cox_basal <- final}
if(subtype[i]=="BRCA.LumA"){lnc_cox_lumA <- final}
if(subtype[i]=="BRCA.LumB"){lnc_cox_lumB <- final}
if(subtype[i]=="BRCA.Normal"){lnc_cox_normal <- final}
}

###################
### For PFI
for( i in 1:length(subtype)){
  print(subtype[i])
  annot <- thorsson[thorsson$TCGA.Subtype==subtype[i],c(1:5,33:36)]
  rownames(annot) <- annot$TCGA.Participant.Barcode
  
  if(subtype[i]=="BRCA.Her2"){spec <- spec_her2}
  if(subtype[i]=="BRCA.Basal"){spec <- spec_basal}
  if(subtype[i]=="BRCA.LumA"){spec <- spec_la}
  if(subtype[i]=="BRCA.LumB"){spec <- spec_lb}
  if(subtype[i]=="BRCA.Normal"){spec <- spec_normal}
  
  gexp <- rna_lnc[spec, rownames(annot)]
  gexp<- (as.data.frame(t(gexp)))
  print(paste0(subtype[i], 
               " Identical= ",
               identical(rownames(gexp), rownames(annot))))
  
  d <- as.vector(colnames(gexp))
  gexp <- cbind(annot,rownames(gexp), gexp)
  # Function for Cox univariate analysis
  coxUnivariateAnalyses <- function(annot, vars4test, time, event, save = TRUE,
                                    fpath = subtype[i]) {
    #-- Rename Survvars
    annot_rn <- annot %>%
      dplyr::rename(time = !! time, event = !! event) %>%
      mutate(across(contains(" mutation"), ~ factor(., levels = c("WT", "Mut"))))
    
    models <- lapply(vars4test, function(varn) {
      message(varn)
      ts_form <- paste0("Surv(time, event) ~ `", varn, "`") %>% as.formula()
      model_res <- coxph(ts_form, data = annot_rn) %>% summary()
    })
    names(models) <- vars4test
    
    var_pvals <- sapply(models, function(model) { model$waldtest[3] })
    
    var_hrs <- lapply(models, function(model) { 
      model$conf.int %>%
        as.data.frame() %>%
        rownames_to_column("Category") %>%
        dplyr::rename(HR = "exp(coef)") %>%
        dplyr::select(-"exp(-coef)")
    }) %>% bind_rows(.id = "Variable")
    
    univar_cox <- data.frame(Variable = names(var_pvals), PFI.Cox_pvalue = var_pvals) %>%
      mutate(PFI.Cox.p.adj = p.adjust(PFI.Cox_pvalue, method = "fdr") ) %>%
      mutate(Variable = str_remove(Variable, ".pvalue"),
             PFI.Significant = ifelse(PFI.Cox.p.adj < 0.1, "*", "ns")) %>%
      arrange(PFI.Cox.p.adj) %>%
      full_join(var_hrs, by = "Variable") %>%
      dplyr::select(Variable, Category, everything())
    colnames(univar_cox)[6:8] <- c("PFI.HR","PFI.lower.95" ,"PFI.upper.95")
    
    if(save) openxlsx::write.xlsx(univar_cox, file = paste0(fpath, "_Cox_results.xlsx"))
    
    return(list(cox_models = models, table = univar_cox))
  }
  # Running analysis
  final <- coxUnivariateAnalyses(gexp,d,time="PFI.Time", event = "PFI", save=F)
  final <- as.data.frame(final$table)
  rownames(final) <- final[,1]
  final <- final[,-1]
  if(subtype[i]=="BRCA.Her2"){
    lnc_cox_her2 <- cbind(lnc_cox_her2,final[rownames(lnc_cox_her2),]) }
  if(subtype[i]=="BRCA.Basal"){
    lnc_cox_basal <- cbind(lnc_cox_basal,final[rownames(lnc_cox_basal),]) }
  if(subtype[i]=="BRCA.LumA"){
    lnc_cox_lumA <- cbind(lnc_cox_lumA,final[rownames(lnc_cox_lumA),]) }
  if(subtype[i]=="BRCA.LumB"){
    lnc_cox_lumB <- cbind(lnc_cox_lumB,final[rownames(lnc_cox_lumB),]) }
  if(subtype[i]=="BRCA.Normal"){
    lnc_cox_normal <- cbind(lnc_cox_normal,final[rownames(lnc_cox_normal),]) }
}

```

```r
## Creating Data.frame for 53 lncRNAs names
longos <- c(lnc_cox_basal[,1], lnc_cox_her2[,1], lnc_cox_lumA[,1],
            lnc_cox_lumB[,1], lnc_cox_normal[,1])
longos.names <- as.data.frame(longos.names[longos.names$`Gene stable ID` %in% longos, ])
rownames(longos.names) <- longos.names$`Gene stable ID`
openxlsx::write.xlsx(longos.names, file="Sup_table_2_lncRNAs_names.xlsx")

### Leukocyte fraction correlation
for(x in 1:length(subtype)){
  if(subtype[x]=="BRCA.Basal"){lnc <- as.data.frame(lnc_cox_basal)}
  if(subtype[x]=="BRCA.LumA"){lnc <- as.data.frame(lnc_cox_lumA)}
  if(subtype[x]=="BRCA.LumB"){lnc <- as.data.frame(lnc_cox_lumB)}
  if(subtype[x]=="BRCA.Normal"){lnc <- as.data.frame(lnc_cox_normal)}
  if(subtype[x]=="BRCA.Her2"){lnc <- as.data.frame(lnc_cox_her2)}
  print(subtype[x])
  
  # Data.frame with Cox results will be used to save leukocyte fraction correlation
  annot <- thorsson[thorsson$TCGA.Subtype==subtype[x],c(1:5,10:14,33:36)]
  rownames(annot) <- annot$TCGA.Participant.Barcode
  regact <- as.data.frame(t(rna_lnc[lnc[,1], rownames(annot)]))
  identical(rownames(regact), rownames(annot))
  regact$Leukocyte.Fraction <- annot$Leukocyte.Fraction
  regact_na <- regact[!is.na(regact$Leukocyte.Fraction),, drop=F]
  ## Correlation matrix
  cor_mat <- regact_na[1,1:(ncol(regact_na)-1)]
  rownames(cor_mat) <- "Leuk.Fraction.Spearman.Corr"
  cor_mat[1,1:ncol(cor_mat)] <- rep(NA,ncol(cor_mat))
  ### Spearman rho inferrence
  for(i in 1:(ncol(regact_na)-1)){
    cor_mat[1,i] <- cor(regact_na[,i], regact_na$Leukocyte.Fraction, method = "spearman")
  }
  #### P-value
  for(i in 1:(ncol(regact_na)-1)){
    cor_mat[2,i] <- cor.test(regact_na[,i], regact_na$Leukocyte.Fraction, method = "spearman")[3]
  }
  #################
  cor_mat <- data.frame(t(cor_mat))
  rownames(lnc) <- lnc[,1]
  cor_mat <- cor_mat[rownames(lnc), , drop=F]
  identical(rownames(cor_mat), rownames(lnc)) #TRUE
  lnc <- cbind(lnc, cor_mat)
  ### Adjusting p-value with FDR method
  lnc$Spear.Corr.P.adj <- p.adjust(lnc[,16], method="fdr")
  colnames(lnc)[c(15,16)] <- c("Leuk.Fraction.corr", "Spear.Corr.P.value")
  #################
  # Adding friendly names to data.frame
  lnc <- cbind(lncRNAs = longos.names[rownames(lnc),"symbol.ens"], lnc)
  
  if(subtype[x]=="BRCA.Her2"){lnc_cox_her2 <- lnc}
  if(subtype[x]=="BRCA.Basal"){lnc_cox_basal <- lnc}
  if(subtype[x]=="BRCA.LumA"){lnc_cox_lumA <- lnc}
  if(subtype[x]=="BRCA.LumB"){lnc_cox_lumB <- lnc}
  if(subtype[x]=="BRCA.Normal"){lnc_cox_normal <- lnc}
}

# Saving data.frame - Supplementary Table 4
sheet <- list(BRCA.Basal = lnc_cox_basal,
              BRCA.LumA = lnc_cox_lumA,
              BRCA.LumB = lnc_cox_lumB,
              BRCA.Normal = lnc_cox_normal,
              BRCA.Her2 = lnc_cox_her2)

writexl::write_xlsx(sheet, path="Sup_Table_4_BRCA_Cox_results.xlsx")

```

Code snippet for Figure 3
----
```r
### Heatmap for Specific lncRNAs
immune_subtype2 <- immune_subtype[-5]
for(x in 1:5){
print(subtype[x])
if(subtype[x]=="BRCA.Her2"){lnc <- as.data.frame(lnc_cox_her2)}
if(subtype[x]=="BRCA.Basal"){lnc <- as.data.frame(lnc_cox_basal)}
if(subtype[x]=="BRCA.LumA"){lnc <- as.data.frame(lnc_cox_lumA)}
if(subtype[x]=="BRCA.LumB"){lnc <- as.data.frame(lnc_cox_lumB)}
if(subtype[x]=="BRCA.Normal"){lnc <- as.data.frame(lnc_cox_normal)}

  # Data.frame with Cox results and Leukocyte Fraction correlation
annot_row <- as.data.frame(lnc)
# 53 lncRNAs expression
dat_imun <- rna_lnc[lnc[,2],]
rownames(annot_row) <- annot_row[,2]
dat_imun <- data.frame(t(dat_imun))
# Data.frame for top annotation
annot_sub <- annot_hm
# Transforming Thorsson et al immune features into categorical data (quartiles)
annot_sub$Leukocyte.Fraction <- as.factor(
  ifelse(annot_sub$Leukocyte.Fraction < quantile(annot_sub$Leukocyte.Fraction, probs = 0.25, na.rm = T),
         "1st quartile",
         ifelse(annot_sub$Leukocyte.Fraction < quantile(annot_sub$Leukocyte.Fraction, probs = 0.50, na.rm = T),
                "2nd quartile",
                ifelse(annot_sub$Leukocyte.Fraction < quantile(annot_sub$Leukocyte.Fraction, probs = 0.75, na.rm = T),
                "3rd quartile",
                "4th quartile"))))
annot_sub$Lymphocyte.Inf.Score <- as.factor(
  ifelse(annot_sub$Lymphocyte.Inf.Score < quantile(annot_sub$Lymphocyte.Inf.Score, probs = 0.25, na.rm = T),
         "1st quartile",
         ifelse(annot_sub$Lymphocyte.Inf.Score < quantile(annot_sub$Lymphocyte.Inf.Score, probs = 0.50, na.rm = T),
                "2nd quartile",
                ifelse(annot_sub$Lymphocyte.Inf.Score < quantile(annot_sub$Lymphocyte.Inf.Score, probs = 0.75, na.rm = T),
                       "3rd quartile",
                       "4th quartile"))))                
annot_sub$Macrophage.Regulation <- as.factor(
  ifelse(annot_sub$Macrophage.Regulation < quantile(annot_sub$Macrophage.Regulation, probs = 0.25, na.rm = T),
         "1st quartile",
         ifelse(annot_sub$Macrophage.Regulation < quantile(annot_sub$Macrophage.Regulation, probs = 0.50, na.rm = T),
                "2nd quartile",
                ifelse(annot_sub$Macrophage.Regulation < quantile(annot_sub$Macrophage.Regulation, probs = 0.75, na.rm = T),
                       "3rd quartile",
                       "4th quartile"))))
annot_sub$Wound.Healing <- as.factor(
  ifelse(annot_sub$Wound.Healing < quantile(annot_sub$Wound.Healing, probs = 0.25, na.rm = T),
         "1st quartile",
         ifelse(annot_sub$Wound.Healing < quantile(annot_sub$Wound.Healing, probs = 0.50, na.rm = T),
                "2nd quartile",
                ifelse(annot_sub$Wound.Healing < quantile(annot_sub$Wound.Healing, probs = 0.75, na.rm = T),
                       "3rd quartile",
                       "4th quartile"))))
annot_sub$IFN.gamma.Response <- as.factor(
  ifelse(annot_sub$IFN.gamma.Response < quantile(annot_sub$IFN.gamma.Response, probs = 0.25, na.rm = T),
         "1st quartile",
         ifelse(annot_sub$IFN.gamma.Response < quantile(annot_sub$IFN.gamma.Response, probs = 0.50, na.rm = T),
                "2nd quartile",
                ifelse(annot_sub$IFN.gamma.Response < quantile(annot_sub$IFN.gamma.Response, probs = 0.75, na.rm = T),
                       "3rd quartile",
                       "4th quartile"))))
annot_sub$TGF.beta.Response <- as.factor(
  ifelse(annot_sub$TGF.beta.Response < quantile(annot_sub$TGF.beta.Response, probs = 0.25, na.rm = T),
         "1st quartile",
         ifelse(annot_sub$TGF.beta.Response < quantile(annot_sub$TGF.beta.Response, probs = 0.50, na.rm = T),
                "2nd quartile",
                ifelse(annot_sub$TGF.beta.Response < quantile(annot_sub$TGF.beta.Response, probs = 0.75, na.rm = T),
                       "3rd quartile",
                       "4th quartile"))))

# Defining Immune subtypes with more than 5 patients for each molecular subtype
annot_sub <- annot_sub[complete.cases(annot_sub$Immune.Subtype),]
annot_sub <- annot_sub[annot_sub$TCGA.Subtype==subtype[x],]
if(subtype[x]=="BRCA.Her2"){lev <- c("C1","C2")}
if(subtype[x]=="BRCA.Basal"){lev <- c("C1","C2")}
if(subtype[x]=="BRCA.LumA"){lev <- c("C1","C2", "C3","C4","C6")}
if(subtype[x]=="BRCA.LumB"){lev <- c("C1","C2", "C3","C4")}
if(subtype[x]=="BRCA.Normal"){lev <- c("C1","C2", "C3","C4","C6")}

summary(annot_sub$Immune.Subtype)
annot_sub <- annot_sub[annot_sub$Immune.Subtype %in% lev,]
dat_imun <- dat_imun[rownames(annot_sub),]
identical(rownames(dat_imun), rownames(annot_sub)) #TRUE
summary(annot_sub$TCGA.Subtype)
#####################
### Column-wise Z-score and setting maximum and minimum standard deviations to plot
dat_imun <- data.frame(scale(dat_imun))
range(dat_imun)

summary(dat_imun)
dat_imun[dat_imun>2]<-2
dat_imun[dat_imun<(-2)]<-(-2)
range(dat_imun)

### Adjusting Survival data for significant results only
identical(rownames(annot_sub), rownames(dat_imun)) #TRUE 
dat_imun <- dat_imun[,rownames(annot_row)]
identical(rownames(annot_row), colnames(dat_imun)) # TRUE
annot_row$OS.Hazard.Ratio <- ifelse(annot_row$OS.HR>1 & annot_row$OS.Cox.p.adj<0.1,
                                 1,
                                 ifelse(annot_row$OS.HR < 1 & annot_row$OS.Cox.p.adj<0.1,-1,0))
annot_row$PFI.Hazard.Ratio <- ifelse(annot_row$PFI.HR>1 & annot_row$PFI.Cox.p.adj<0.1,
                                    1,
                                    ifelse(annot_row$PFI.HR < 1 & annot_row$PFI.Cox.p.adj<0.1,-1,0))

colnames(annot_row)[16] <- "LF\ncorrelation"

#### Changing lncRNAs names to symbols
longos <- as.data.frame(longos.names[longos.names$`Gene stable ID` %in% colnames(dat_imun),])
rownames(longos) <- longos$`Gene stable ID`
longos <- longos[colnames(dat_imun),]
identical(rownames(longos), colnames(dat_imun)) #true
colnames(dat_imun) <- longos$symbol.ens
identical(rownames(longos), rownames(annot_row)) #true
rownames(annot_row) <- longos$symbol.ens

identical(rownames(annot_row), colnames(dat_imun)) # TRUE
#############

# Heatmap

h1<- Heatmap(matrix = t(dat_imun),
             name= "lncRNA expression\nz-score",
             show_heatmap_legend = F,
             show_column_names = F, 
             show_column_dend = F,
             column_names_side="bottom",
             column_names_gp= gpar(col="white", fontsize=1),
             column_title_gp = gpar(fontsize=8, fontface="bold"), 
             # column_title = c("Basal\n(n=169)", "Normal\n(n=136)",
             #                  "LumA\n(n=499)", "LumB\n(n=184)","Her2\n(n=72)"),
             column_split = factor(annot_sub$Immune.Subtype,
                                   levels=c("C1","C2","C3","C4","C5","C6" )),
             cluster_columns = T, 
             cluster_column_slices = F,
             cluster_rows = T, 
             show_row_names = T,
             show_row_dend = T,
             row_dend_width = unit(5, "mm"),
             #row_split = split, 
             row_title_rot = 0, 
             row_title_gp = gpar(fontsize=8),
             row_names_gp = gpar(fontsize=7),
             row_title_side ="left",
             heatmap_legend_param = list(
               at=c(-2,0,2),
               legend_width= unit(1.7,"cm"),
               labels_gp=gpar(fontsize=8), 
               title_gp=gpar(fontsize=8, fontface="bold"),
               border=T, 
               direction = "horizontal"),
            top_annotation = columnAnnotation(
              df=annot_sub[,-c(2)],
              simple_anno_size= unit(2.5, "mm"),
              height=unit(1,"cm"),
              annotation_name_gp = gpar(fontsize=8),
              annotation_name_side = "right",
              #show_legend=c(T,T,T,F,F,F,F,F),
              show_legend=F,
              na_col="white", 
              annotation_legend_param = list(
                TCGA.Subtype=list(
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=3, 
                  title="TCGA Subtype"),
                Immune.Subtype=list(
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=3, 
                  title="Immune\nSubtype"),
                Leukocyte.Fraction=list(
                  title="Leukocyte\nFraction",
                  border=T, #legend_width= unit(1.7,"cm"),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=4,
                  direction="horizontal"),
                Lymphocyte.Inf.Score=list(
                  title="Leukocyte\nFraction",
                  border=T, #legend_width= unit(1.7,"cm"),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=2,
                  direction="horizontal"),
                Macrophage.Regulation=list(
                  title="Leukocyte\nFraction",
                  border=T, #legend_width= unit(1.7,"cm"),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=2,
                  direction="horizontal"),
                Wound.Healing=list(
                  title="Leukocyte\nFraction",
                  border=T, #legend_width= unit(1.7,"cm"),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=2,
                  direction="horizontal"),
                IFN.gamma.Response=list(
                  title="Leukocyte\nFraction",
                  border=T, #legend_width= unit(1.7,"cm"),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=2,
                  direction="horizontal"),
                TGF.beta.Response=list(
                  title="Leukocyte\nFraction",
                  border=T, #legend_width= unit(1.7,"cm"),
                  labels_gp=gpar(fontsize=8),
                  title_gp=gpar(fontsize=8, fontface="bold"),
                  nrow=2,
                  direction="horizontal")
              ),
              col=list(
                Immune.Subtype=c("C1"="red","C2"="yellow","C3"="green3",
                                 "C4"="cyan","C5"="blue","C6"="pink"),
                TCGA.Subtype=c("BRCA.Basal"="sienna1", "BRCA.Normal"="purple3",
                               "BRCA.Her2" ="hotpink", "BRCA.LumA"="cornflowerblue",
                               "BRCA.LumB"="springgreen2"),
                Leukocyte.Fraction=c("1st quartile"="gray90", "2nd quartile"="gray60",
                                     "3rd quartile"= "gray40", "4th quartile"= "black") ,
                Lymphocyte.Inf.Score=c("1st quartile"="gray90", "2nd quartile"="gray60",
                                       "3rd quartile"= "gray40", "4th quartile"= "black") ,
                Macrophage.Regulation =c("1st quartile"="gray90", "2nd quartile"="gray60",
                                         "3rd quartile"= "gray40", "4th quartile"= "black") ,
                Wound.Healing=c("1st quartile"="gray90", "2nd quartile"="gray60",
                                "3rd quartile"= "gray40", "4th quartile"= "black") ,
                IFN.gamma.Response=c("1st quartile"="gray90", "2nd quartile"="gray60",
                                     "3rd quartile"= "gray40", "4th quartile"= "black") ,
                TGF.beta.Response=c("1st quartile"="gray90", "2nd quartile"="gray60",
                                    "3rd quartile"= "gray40", "4th quartile"= "black"))
              ),
            right_annotation = rowAnnotation(
              df=data.frame(annot_row[,c("OS.Hazard.Ratio","PFI.Hazard.Ratio"), drop=F]),
              annotation_name_rot=90,
              #show_legend=c(T,F),
              show_legend=c(F),
              annotation_name_side="bottom",
              annotation_legend_param = list(
                title= "Hazard Ratio",
                at=c(-1,0,1),
                labels=c("<1","NS",">1"),
                labels_gp=gpar(fontsize=8),
                title_gp=gpar(fontsize=8, fontface="bold"),
                border=T,
                #legend_width= unit(1.7,"cm"),
                direction="horizontal",
                nrow=3),
              annotation_label = c("OS","PFI"),
              annotation_name_gp = gpar(fontsize=8) ,
              simple_anno_size = unit(0.2, "cm"),
              #gap=unit(3.5,"mm"),
              col=list(OS.Hazard.Ratio = c("-1"= "green4","0" ="grey", "1"= "magenta4"),
                       PFI.Hazard.Ratio = c("-1"= "green4","0" ="grey", "1"= "magenta4"))),
            left_annotation = rowAnnotation(
              df=data.frame(annot_row[,"LF\ncorrelation", drop=F]),
              show_legend=F,
              annotation_name_rot=90,
              annotation_name_side="bottom",
              annotation_legend_param = list(
                title="LF Spearman\nCorrelation",
                at=c(-0.8,0,0.8),
                #labels=c("<1","NS",">1"),
                labels_gp=gpar(fontsize=8),
                title_gp=gpar(fontsize=8, fontface="bold"),
                border=T,
                #legend_width= unit(1.7,"cm"),
                direction="horizontal"),
              annotation_label = "LF\ncorr.",
              annotation_name_gp = gpar(fontsize=8) ,
              simple_anno_size = unit(0.2, "cm"),
              col=list("LF.correlation" = circlize::colorRamp2(c(-1,0,1), 
                                                    c("darkgreen","white", "orangered4")))))
        
g<- grid.grabExpr(
  draw(h1, 
       heatmap_legend_side="bottom",
       merge_legends=F))
plot_grid(g)

if(subtype[x]=="BRCA.Her2"){g.her2 <- plot_grid(g)}
if(subtype[x]=="BRCA.Basal"){g.basal <- plot_grid(g)}
if(subtype[x]=="BRCA.LumA"){g.lumA <- plot_grid(g)}
if(subtype[x]=="BRCA.LumB"){g.lumB <- plot_grid(g)}
if(subtype[x]=="BRCA.Normal"){g.normal <- plot_grid(g)}
}
row1 <- plot_grid(g.lumA)
row2 <- plot_grid(g.lumB, g.her2, rel_widths = c(1,0.7))
row3 <- plot_grid(g.normal, g.basal, rel_widths = c(1,0.7))
g <-plot_grid(row1,row2,row3, nrow = 3)
g
ggsave(plot=g,filename = "Heatmap_lncs_v5.svg",device = "svg",
       dpi=300, height = 200, width=180, units="mm")
# Cleaning environment
rm(annot_row, annot_sub, cor_mat, 
   dat_imun, final, g.basal, g.her2, g.lumA, g.lumB, g.normal, gexp, h1, ha, lnc, lnc_cox,
   lnc_nao, lnc.b, lnc.h, lnc.la, lnc.lb, lnc.n, row1, row2, row3, spec,
   spec_basal, spec_her2, spec_la, spec_lb, spec_normal, basal, d, her2, i, immune_subtype2,
   lev, lnc_sub, lumA, lumB, normal, quantil, x)

```

GSEA analysis
----
```r
# Function to plot GSEA results from prof Mauro Castro
library(purrr)
library(dplyr)
library(ggplot2)

my_dot_plot_gsea <- function(gsea_obj, top=20, abrv=60, fdr_level=0.05, 
                             title="", ylab="Hallmark gene sets") {
  
  # get fdr threshold
  fdr_threshold <- p.threshold(gsea_obj$pval, alpha = fdr_level, method = "BH")
  
  # adjust pval
  gsea_obj$padj <- p.adjust(gsea_obj$pval, method = "BH")
  
  # Subset results
  df <- gsea_obj %>%
    purrr::when(!is.null(top) ~ head(., top), ~ .)
  df <- as.data.frame(df)
  df <- df[,c("pathway","size", "ES","NES", "pval", "padj")]
  
  # Handle empty dataframes
  if (nrow(df) == 0) return(ggempty())
  
  # Plotting variables
  df$significance <- df$pval
  
  # Order by significance value
  df <- df[order(-df$pval),]
  
  # Abbreviate labels
  label.abrv <- substr(df$pathway, 1, abrv)
  if (any(duplicated(label.abrv))) {
    stop("Non-unique labels after abbreviating")
  } else {
    df$label.abrv <- factor(label.abrv, levels=label.abrv)   
  }
  .reverselog_trans <- function(base=exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
                      scales::log_breaks(base=base), 
                      domain=c(1e-100, Inf))
  }
  df$effect <- c("<0","=0",">0")[sign(df$NES)+2]
  df$effect[df$padj > fdr_level] <- "=0"
  df$effect <- factor(df$effect, levels = c("<0","=0",">0"))
  cols <- c("<0" = "blue", "=0"="grey",">0" = "red")
  #---at
  at <- -log10(df$pval)
  at <- sort(c(range(at),mean(range(at))))
  at <- 10^-at
  if(all(at>fdr_level))at <- sort(c(max(at),fdr_level))
  #---size
  bksz <- sort(unique(df$size))
  bksz <- pretty(bksz)
  if(length(bksz)>5){
    bksz <- sort(unique(as.numeric(sapply(bksz, pretty))))
    bksz <- unique(round(quantile(bksz)))
  }
  #---plot paste0
  gp <- ggplot(df, aes(y=label.abrv, x=significance, size=size, fill = effect, )) +
    geom_vline(xintercept=fdr_threshold, linetype=2, color="#800026") +
    geom_point(show.legend = c(size=T, colour=T), alpha=0.7, shape = 21) +
    scale_fill_manual(values = cols, aesthetics="fill", drop=F) +
    # scale_colour_manual(values = cols, aesthetics="fill") +
    scale_size(range = c(1, 6), breaks = bksz) +
    labs(size="Gene set size",fill="Enrichment (NES)", title=title) +
    scale_x_continuous(trans=.reverselog_trans(10), breaks = at,
                       expand = expansion(mult = c(0.07,0.1)), 
                       labels = format(-log10(at), digits = 1)) +
    scale_y_discrete(expand = expansion(mult = c(0.05,0.07)), position = "right") +
    ylab(paste0(ylab," (n=",nrow(df),")")) + 
    xlab("-Log10(P-value)") + 
    annotate(geom="text", y=(top), x=fdr_threshold-(fdr_threshold/10), label=paste0("FDR<",fdr_level), 
             color="#800026", size=3, hjust = 0, fontface="bold") +
    theme(
      plot.title=element_text(hjust=0),
      legend.position="bottom",
      legend.margin = margin(-5, -5, 5, -5),
      legend.spacing.y = unit(0, 'cm'),
      title = element_text(size=12), 
      axis.title.y = element_text(size=12), 
      axis.text.y=element_text(size=8), 
      axis.ticks = element_line(linetype=1, color="grey"),
      legend.text=element_text(size=7),
      legend.title=element_text(size=8),
      # aspect.ratio=2/1,
      # axis.line = element_blank(),
      # panel.grid.major.y = element_blank(),
      # panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype=111, color="grey80", size = 0.4),
      panel.grid.major = element_line(linetype=3, color="grey", size = 0.2),
      panel.background = element_rect(fill = "grey98", colour = "grey50"),
      panel.border = element_rect(colour = "grey", fill=NA, size=1)) + 
    guides(size = guide_legend(label.position = "bottom", title.position="bottom", title.hjust=0.5)) +
    guides(fill = guide_legend(label.position = "bottom", title.position="bottom", title.hjust=0.5))
  return(gp)
}

##------------------------------------------------------------------------
##returns rejection threshold for methods in 'p.adjust'
p.threshold <- function(pvals, alpha=0.05, method="BH") {
  pvals <- sort(pvals)
  padj <- p.adjust(pvals, method = method)
  thset <- which(padj <= alpha)
  if(length(thset)>0){
    mx1 <- mx2 <- which.max(thset)
    if(mx2<length(padj)) mx2 <- mx2 + 1
    th <- (pvals[mx1] + min(pvals[mx2],alpha) ) / 2
  } else {
    th <- min(c(alpha,pvals))
  }
  return(th)
}
```


```r
# Cleaning environment
rm(annot, g, longos, regact,regact_na, sheet, soma_basal, soma_her2, soma_lumA, soma_lumB, soma_normal, soma2, stnr, sup.tab, sup.tab2, mess)

# Enrichment analysis (GSEA)

###### Choosing Hallmarks modules for enrichment
h_gene_sets = msigdbr(species = "Homo sapiens", 
                      category = "H")
head(h_gene_sets)
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

#### Preparing matrix with all subtypes including healthy
gexp_lnc <- as.data.frame(t(rna_lnc))
identical(thorsson$TCGA.Participant.Barcode, rownames(gexp_lnc)) #TRUE
gexp_lnc <- cbind(TCGA.Subtype= thorsson[,c("TCGA.Subtype")],
                  gexp_lnc)
gexp_healthy <- as.data.frame(t(rna_lnc_healthy))
gexp_healthy <- gexp_healthy[,rownames(rna_lnc)]
gexp_healthy <- cbind(TCGA.Subtype= rep("Non-tumoral",nrow(gexp_healthy)),
                      gexp_healthy)
identical(colnames(gexp_healthy), colnames(gexp_lnc)) #TRUE
gexp_lnc <- rbind(gexp_lnc, gexp_healthy)
rm(gexp_healthy)
gexp_lnc <- gexp_lnc[complete.cases(gexp_lnc$TCGA.Subtype),]
gexp_lnc[,2:ncol(gexp_lnc)] <- lapply(gexp_lnc[,2:ncol(gexp_lnc)], as.numeric)
gexp_lnc$TCGA.Subtype <- factor(gexp_lnc$TCGA.Subtype,
                                   levels = c("BRCA.Basal","BRCA.LumA" ,
                                              "BRCA.LumB" ,"BRCA.Normal",
                                              "BRCA.Her2", "Non-tumoral"))

#####################################################################################
# Looping for enrichment
#####################################################################################
for(x in 1:5){
  lnc <- NULL
  if(subtype[x]=="BRCA.Basal"){lnc <- lnc_cox_basal[,2, drop=T]}
  if(subtype[x]=="BRCA.LumA"){lnc <- lnc_cox_lumA[,2, drop=T]}
  if(subtype[x]=="BRCA.LumB"){lnc <- lnc_cox_lumB[,2, drop=T]}
  if(subtype[x]=="BRCA.Normal"){lnc <-lnc_cox_normal[,2, drop=T]}
  if(subtype[x]=="BRCA.Her2"){lnc <- lnc_cox_her2[,2, drop=T]}
  print(subtype[x])
  sheet <- list(NULL)
  for(i in 1:length(lnc)){
    varn <- lnc[i]
    print(subtype[x]," ", varn)

#####################################################################################
# Preparing data for enrichment
#####################################################################################
    annot_sub <- thorsson[thorsson$TCGA.Subtype==subtype[x],]
    rownames(annot_sub) <- annot_sub$TCGA.Participant.Barcode
    longo <- as.data.frame(t(rna_lnc[varn,rownames(annot_sub), drop=F]))
    gexp <- rna[,rownames(annot_sub)]
    gexp <- gexp[rowSums(gexp, na.rm = T)>0,]
    
    longo$Actv <- factor(ifelse(longo[,varn]>median(longo[,varn]), "High", "Low"),
                            levels=c("High", "Low"))
    
    identical(rownames(longo), colnames(gexp)) #TRUE
    
    ###########
    # Differential expression by Signal to Noise Ratio (SNR)
    high <- rownames(longo[longo$Actv=="High",])
    print("SNR calculus")
    snr <- (apply(gexp[,names(gexp) %in% high],1,mean_na)-apply(gexp[,!(names(gexp)%in%high)],1,mean_na))/
      (apply(gexp[,names(gexp)%in%high],1,sd_na)+apply(gexp[,!(names(gexp)%in%high)],1,sd_na))
    snr <- as.data.frame(snr[order(snr, decreasing = T)])  
    colnames(snr) <- "SNR.High.vs.Low"   
    
    #####################################################################################
    # Enrichment analysis
    #####################################################################################
    
    ### Ordering genes
    genes.order <- gene.names
    genes.order <- genes.order[!duplicated(genes.order$`Gene stable ID`),]
    rownames(genes.order) <- genes.order$`Gene stable ID`
    genes.order <- genes.order[rownames(snr),]
    genes.order <- cbind(genes.order[,2:3], rownames(snr),snr)
    genes.order <- genes.order[genes.order$`HGNC symbol`!="",]
    dup <- genes.order$`HGNC symbol`[duplicated(genes.order$`HGNC symbol`)]
    genes.order$`HGNC symbol`[match(dup, genes.order$`HGNC symbol`)] <- paste0(dup,"_")
    genes.order <- genes.order[complete.cases(genes.order$`HGNC symbol`),]
    rownames(genes.order) <- genes.order$`HGNC symbol`
    gen <- as.vector(genes.order$SNR.High.vs.Low)
    names(gen) <- genes.order$`HGNC symbol`
    
    ##### enrichment analysis with fgsea
    print("GSEA")
    df.gsea <- fgsea(pathways = msigdbr_list, stats= gen, eps=0, nPermSimple=10000)
    print("GSEA concluded")
    gsea_obj <- df.gsea
    gsea_obj <- gsea_obj[order(gsea_obj$pval),]
    gsea_obj$leadingEdge <- as.character(gsea_obj$leadingEdge)
    
    sheet[[i]] <- gsea_obj
    
    ######
    # Plotting results for selected lncRNAs
    gsea_obj <- gsea_obj[complete.cases(gsea_obj$pval),]
    gsea_obj$pathway <- substr(gsea_obj$pathway,10, nchar(gsea_obj$pathway))
    if(varn %in% c("ENSG00000235576","ENSG00000214548","ENSG00000281649",
                   "ENSG00000231367", "ENSG00000230266" )){
     gsea_plot <-my_dot_plot_gsea(gsea_obj, 
                       top=20,
                       title = paste0(subtype[x],"\n",varn," High vs Low"))
      ggsave(plot= gsea_plot, 
             filename = paste0("GSEA_",subtype[x],"_",varn,"_Hallmark.svg"), 
             device = "svg",
             dpi=300, height = 120 , width = 120, units="mm")
    }
  }
  
  lnc_name <- as.data.frame(longos.names[longos.names$`Gene stable ID` %in% lnc,])
  rownames(lnc_name) <- lnc_name$`Gene stable ID`
  lnc_name <- lnc_name[lnc,]
  names(sheet) <- lnc_name$symbol.ens
  ### Saving enrichment results for each molecular subtype
  writexl::write_xlsx(sheet, path = paste0(subtype[x],"_GSEA_all_H.xlsx"))
  print("Done")
}
```

Code snippet for Supplementary Figure 3
----
```r
### Plotting GSEA tables
hallmark <- data.frame(
        stringsAsFactors = FALSE,
                               Number = c(1L,2L,3L,4L,5L,6L,7L,8L,
                                          9L,10L,11L,12L,13L,14L,15L,16L,
                                          17L,18L,19L,20L,21L,22L,23L,
                                          24L,25L,26L,27L,28L,29L,30L,31L,
                                          32L,33L,34L,35L,36L,37L,38L,
                                          39L,40L,41L,42L,43L,44L,45L,46L,
                                          47L,48L,49L,50L),
                       Hallmark.Name = c("APICAL_JUNCTION",
                                          "APICAL_SURFACE","PEROXISOME","ADIPOGENESIS",
                                          "ANGIOGENESIS",
                                          "EPITHELIAL_MESENCHYMAL_TRANSITION","MYOGENESIS",
                                          "SPERMATOGENESIS","PANCREAS_BETA_CELLS",
                                          "DNA_REPAIR","UV_RESPONSE_DN",
                                          "UV_RESPONSE_UP","ALLOGRAFT_REJECTION",
                                          "COAGULATION","COMPLEMENT",
                                          "INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE",
                                          "IL6_JAK_STAT3_SIGNALING",
                                          "INFLAMMATORY_RESPONSE","BILE_ACID_METABOLISM",
                                          "CHOLESTEROL_HOMEOSTASIS",
                                          "FATTY_ACID_METABOLISM","GLYCOLYSIS",
                                          "HEME_METABOLISM","OXIDATIVE_PHOSPHORYLATION",
                                          "XENOBIOTIC_METABOLISM","APOPTOSIS",
                                          "HYPOXIA","PROTEIN_SECRETION",
                                          "UNFOLDED_PROTEIN_RESPONSE",
                                          "REACTIVE_OXYGEN_SPECIES_PATHWAY","E2F_TARGETS",
                                          "G2M_CHECKPOINT","MYC_TARGETS_V1",
                                          "MYC_TARGETS_V2","P53_PATHWAY",
                                          "MITOTIC_SPINDLE","ANDROGEN_RESPONSE",
                                          "ESTROGEN_RESPONSE_EARLY",
                                          "ESTROGEN_RESPONSE_LATE","IL2_STAT5_SIGNALING",
                                          "KRAS_SIGNALING_UP","KRAS_SIGNALING_DN",
                                          "MTORC1_SIGNALING","NOTCH_SIGNALING",
                                          "PI3K_AKT_MTOR_SIGNALING",
                                          "HEDGEHOG_SIGNALING","TGF_BETA_SIGNALING",
                                          "TNFA_SIGNALING_VIA_NFKB",
                                          "WNT_BETA_CATENIN_SIGNALING"),
                     Process.Category = c("cellular component",
                                          "cellular component","cellular component",
                                          "development","development",
                                          "development","development","development",
                                          "development","DNA damage","DNA damage",
                                          "DNA damage","immune","immune",
                                          "immune","immune","immune","immune",
                                          "immune","metabolic","metabolic",
                                          "metabolic","metabolic","metabolic",
                                          "metabolic","metabolic","pathway",
                                          "pathway","pathway","pathway",
                                          "pathway","proliferation","proliferation",
                                          "proliferation","proliferation",
                                          "proliferation","proliferation",
                                          "signaling","signaling","signaling",
                                          "signaling","signaling","signaling",
                                          "signaling","signaling","signaling",
                                          "signaling","signaling","signaling",
                                          "signaling"),
                          Description = c("apical junction complex consisting of adherens and tight junctions",
                                          "membrane proteins in the apical domain","peroxisomes",
                                          "adipocyte development","blood vessel formation",
                                          "epithelial mesenchymal transition",
                                          "muscle differentiation",
                                          "sperm development and male fertility",
                                          "genes specific to pancreatic beta cells",
                                          "DNA repair","UV response: downregulated genes",
                                          "UV response: upregulated genes",
                                          "allograft rejection",
                                          "blood coagulation cascade","complement cascade",
                                          "interferon alpha response",
                                          "interferon gamma response",
                                          "IL6 STAT3 signaling during acute phase response",
                                          "inflammation","biosynthesis of bile acids",
                                          "cholesterol homeostasis",
                                          "fatty acid metabolism",
                                          "glycolysis and gluconeogenesis",
                                          "heme metabolism and erythroid lineage",
                                          "oxidative phosphorylation and citric acid cycle",
                                          "metabolism of xenobiotics",
                                          "programmed cell death; caspase pathway",
                                          "response to hypoxia; HIF1A targets","protein secretion",
                                          "unfolded protein response; ER stress","reactive oxygen species pathway",
                                          "cell cycle progression: E2F targets",
                                          "cell cycle progression: G2/M checkpoint","MYC targets, variant 1",
                                          "MYC targets, variant 2","p53 pathway",
                                          "cell cycle progression: mitotic spindle assembly","androgen response",
                                          "early estrogen response",
                                          "late estrogen response","IL2 STAT5 signaling",
                                          "KRAS signaling, upregulated genes",
                                          "KRAS signaling, downregulated genes",
                                          "mTORC1 signaling","Notch signaling",
                                          "PI3K signaling via AKT to mTORC1",
                                          "Hedgehog signaling","TGF beta signaling",
                                          "TNFA signaling via NFκB",
                                          "canonical beta catenin pathway"),
               Number.of.Founder.Sets = c(37L,12L,28L,36L,14L,107L,
                                          64L,24L,24L,44L,17L,16L,190L,
                                          71L,71L,82L,82L,24L,120L,28L,28L,
                                          53L,87L,36L,93L,124L,80L,87L,
                                          74L,22L,13L,420L,420L,404L,6L,
                                          85L,108L,8L,61L,61L,13L,14L,16L,
                                          487L,49L,591L,79L,29L,132L,49L),
                      Number.of.Genes = c(200L,44L,107L,200L,36L,
                                          200L,200L,135L,40L,150L,144L,158L,
                                          200L,138L,200L,97L,200L,87L,200L,
                                          112L,74L,158L,200L,200L,200L,
                                          200L,161L,200L,96L,113L,49L,200L,
                                          200L,200L,58L,200L,200L,117L,
                                          200L,200L,200L,200L,200L,200L,32L,
                                          105L,36L,54L,200L,42L)
             )

hallmark <- as.data.frame(readxl::read_xlsx("Hallmarks.xlsx"))
colnames(hallmark) <- c("Number", "Hallmark Name", "Process Category", "Description",
                         "Number of Founder Sets", "Number of Genes")
rownames(hallmark) <- hallmark$`Hallmark Name`

### Preparing GSEA results data.frame for plotting
for(x in 1:5){
  print(subtype[x])
  lnc <- list(NULL)
  if(subtype[x]=="BRCA.Basal"){
    z=2
    cols <- c("LY6E-DT","U62317.1","MIR3142HG","LINC01871","PELATON","ETV7-AS1",
              "U62317.3","AC007991.2","AC008760.2","APCDD1L-DT")
  for(i in 1:length(cols)){
    lnc[[i]] <- readxl::read_xlsx("BRCA.Basal_GSEA_all_H.xlsx", sheet=i)
    lnc[[i]]$Lnc <- rep(cols[i], nrow(lnc[[i]]))
    lnc[[i]]$pathway <- substr(lnc[[i]]$pathway,10, nchar(lnc[[i]]$pathway))
    lnc[[i]]$pathway.class <- hallmark[lnc[[i]]$pathway,"Process Category"]
    lnc[[i]] <- lnc[[i]][lnc[[i]]$padj<0.05 & complete.cases(lnc[[i]]),]
  }
  names(lnc) <- cols
  lnc.df <- as.data.frame(lnc[[1]])
  for(i in 2:length(cols)){lnc.df <- rbind(lnc.df, lnc[[i]])}
}
if(subtype[x]=="BRCA.LumA"){
  z=4
  cols <- c("MAP3K4-AS1", "AL049838.1", "DNM3OS","MEG3", "LINC02544",
            "LINC01711" , "LINC01638","AP001189.1","EWSAT1","LINC00271" , "AC105285.1")
  for(i in 1:length(cols)){
    lnc[[i]] <- readxl::read_xlsx("BRCA.LumA_GSEA_all_H.xlsx", sheet=i)
    lnc[[i]]$Lnc <- rep(cols[i], nrow(lnc[[i]]))
    lnc[[i]]$pathway <- substr(lnc[[i]]$pathway,10, nchar(lnc[[i]]$pathway))
    lnc[[i]]$pathway.class <- hallmark[lnc[[i]]$pathway,"Process Category"]
    lnc[[i]] <- lnc[[i]][lnc[[i]]$padj<0.05 & complete.cases(lnc[[i]]),]
  }
  names(lnc) <- cols
  lnc.df <- as.data.frame(lnc[[1]])
  for(i in 2:length(cols)){lnc.df <- rbind(lnc.df, lnc[[i]])}
  }
if(subtype[x]=="BRCA.LumB"){
  z=5
  cols <- c("EBLN3P", "USP27X-AS1","AC096733.2","FAM160A1-DT","AL035661.1",
            "LINC00957","AC092718.3","AC009119.1","LINC02620","AL445490.1",
            "AC104984.4")
  for(i in 1:length(cols)){
    lnc[[i]] <- readxl::read_xlsx("BRCA.LumB_GSEA_all_H.xlsx", sheet=i)
    lnc[[i]]$Lnc <- rep(cols[i], nrow(lnc[[i]]))
    lnc[[i]]$pathway <- substr(lnc[[i]]$pathway,10, nchar(lnc[[i]]$pathway))
    lnc[[i]]$pathway.class <- hallmark[lnc[[i]]$pathway,"Process Category"]
    lnc[[i]] <- lnc[[i]][lnc[[i]]$padj<0.05 & complete.cases(lnc[[i]]),]
  }
  names(lnc) <- cols
  lnc.df <- as.data.frame(lnc[[1]])
  for(i in 2:length(cols)){lnc.df <- rbind(lnc.df, lnc[[i]])}
  }
if(subtype[x]=="BRCA.Normal"){
  z=3
  cols <- c("HLX-AS1","AL133371.2","AC008957.1","RXYLT1-AS1","AL162171.1",
            "LINC02613","AC092164.1","LRRC8C-DT","AC107959.1","LINC02660","LINC02723")
  for(i in 1:length(cols)){
    lnc[[i]] <- readxl::read_xlsx("BRCA.Normal_GSEA_all_H.xlsx", sheet=i)
    lnc[[i]]$Lnc <- rep(cols[i], nrow(lnc[[i]]))
    lnc[[i]]$pathway <- substr(lnc[[i]]$pathway,10, nchar(lnc[[i]]$pathway))
    lnc[[i]]$pathway.class <- hallmark[lnc[[i]]$pathway,"Process Category"]
    lnc[[i]] <- lnc[[i]][lnc[[i]]$padj<0.05 & complete.cases(lnc[[i]]),]
    
  }
  names(lnc) <- cols
  lnc.df <- as.data.frame(lnc[[1]])
  for(i in 2:length(cols)){lnc.df <- rbind(lnc.df, lnc[[i]])}
  }
if(subtype[x]=="BRCA.Her2"){
  z=6
  cols <- c("LINC02694","XXYLT1-AS2","LINC02446","USP30-AS1", "AL365361.1",
            "LINC02384","LINC01857", "AC093583.1","LINC02528","LINC02073")
  for(i in 1:length(cols)){
    lnc[[i]] <- readxl::read_xlsx("BRCA.Her2_GSEA_all_H.xlsx", sheet=i)
    lnc[[i]]$Lnc <- rep(cols[i], nrow(lnc[[i]]))
    lnc[[i]]$pathway <- substr(lnc[[i]]$pathway,10, nchar(lnc[[i]]$pathway))
    lnc[[i]]$pathway.class <- hallmark[lnc[[i]]$pathway,"Process Category"]
    lnc[[i]] <- lnc[[i]][lnc[[i]]$padj<0.05 & complete.cases(lnc[[i]]),]
  }
  names(lnc) <- cols
  lnc.df <- as.data.frame(lnc[[1]])
  for(i in 2:length(cols)){lnc.df <- rbind(lnc.df, lnc[[i]])}
  
  }
## Organizing Hallmarks classes
lnc.df$pathway.class <- factor(lnc.df$pathway.class,
                               levels = c("immune", "signaling","proliferation", 
                                          "cellular component", "development" ,
                                          "metabolic" , "pathway","DNA damage"))
### Defining positive or negative enrichment
lnc.df$effect <- ifelse(lnc.df$NES>1,">1","<1")
### Preparing data.frame
lnc.df <- lnc.df[order(lnc.df$pathway.class, lnc.df$pathway),]
lnc.df$pathway <- factor(lnc.df$pathway, levels=c(unique(rev(lnc.df$pathway))))
lnc.df$Lnc <- factor(lnc.df$Lnc, levels=c(cols))
### Plotting GSEA tables
g<- ggplot(lnc.df) +
  geom_rug(aes(y=pathway, color=pathway.class, x=NULL, fill=NULL),size=3)+
  geom_point(aes(y=pathway, x=Lnc, size=-log10(padj), fill = effect ),
             show.legend = c(size=T, colour=T), alpha=0.7, shape = 21) +
  scale_fill_manual(values = c("blue","red"), aesthetics="fill", drop=F)+
  scale_color_brewer(type="qual",palette="Accent")+
  labs(size="-log10 (padj)",fill="Enrichment\n(NES)", color="Process Category",
       title= paste0(subtype[x]," lncRNAs GSEA" )) +
  scale_x_discrete(expand = expansion(mult = c(0.1,0.05))) +
  scale_y_discrete(expand = expansion(mult = c(0.02,0.02))) +
  #scale_size(breaks=c(1.3,  40, 80)) +
  ylab("Hallmarks gene sets") + 
  xlab("lncRNAs")+
  theme(
    plot.title=element_text(hjust=0),
    legend.position="right",
    legend.margin = margin(-5, -5, 5, -5),
    legend.spacing.y = unit(0.3, 'cm'),
    title = element_text(size=12), 
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
    axis.title.y = element_text(size=12), 
    axis.text.y=element_text(size=8), 
    axis.ticks = element_line(linetype=1, color="grey"),
    legend.text=element_text(size=7),
    legend.title=element_text(size=8, face = "bold"),
    # aspect.ratio=2/1,
    # axis.line = element_blank(),
    # panel.grid.major.y = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype=111, color="grey80", size = 0.4),
    panel.grid.major = element_line(linetype=3, color="grey", size = 0.2),
    panel.background = element_rect(fill = "grey98", colour = "grey50"),
    panel.border = element_rect(colour = "grey", fill=NA, size=1)) +
  guides(size = guide_legend(title.hjust=0.5)) +
  guides(fill = guide_legend( title.hjust=0.5))
g
### Saving plots
ggsave(plot= g, 
       filename= paste0("Supplementar Figure ",z, " ", subtype[x],"_GSEA_plot.pdf"),
       device="pdf",
        dpi=300, width=180,
        height =180, units="mm")  
}

```

Code snippet for lncRNAs panels. Figures 4-8 and Supplementary Figure 4
----
```r
#####################################################################################
# Boxplot and Kaplan-Meier for selected lncRNAs
#####################################################################################
for(x in 1:5){
  print(subtype[x])
  # Defining lncRNA for each molecular subtype
  if(subtype[x]=="BRCA.Basal"){varn <- "ENSG00000235576"}
  if(subtype[x]=="BRCA.LumA"){varn <- "ENSG00000214548"}
  if(subtype[x]=="BRCA.LumB"){varn <- "ENSG00000281649"}
  if(subtype[x]=="BRCA.Normal"){varn <-"ENSG00000231367"}
  if(subtype[x]=="BRCA.Her2"){varn <- "ENSG00000230266"}
  
  print(varn)
  ### Preparting data.frame with molecular subtype patient
  annot_sub <- thorsson[thorsson$TCGA.Subtype==subtype[x],]
  rownames(annot_sub) <- annot_sub$TCGA.Participant.Barcode
  longo <- as.data.frame(t(rna_lnc[varn,rownames(annot_sub), drop=F]))
  longo$Actv <- factor(ifelse(longo[,varn]>median(longo[,varn]), "High", "Low"),
                       levels=c("High", "Low"))
    
  identical(rownames(longo), colnames(gexp)) #TRUE
  
  identical(rownames(longo), rownames(annot_sub))
  annot <- cbind(annot_sub,longo)
  
#####################################################################################
# Boxplot for immune subtypes
#####################################################################################
  annot2 <- annot[annot$Immune.Subtype %in% c(immune_subtype),]
  annot2$LNG <- annot2[,varn]
  summary(annot2$Immune.Subtype)
  
  if(subtype[x]=="BRCA.Basal" | subtype[x]=="BRCA.Her2"){
    # Define comparison to plot
    comparison <- list(c("C1","C2") )
    annot2<- annot2[annot2$Immune.Subtype %in% c("C1","C2"),]
  
    ## Plotting
    library(ggpubr)
    ymax <- max(annot2[,varn])
    g2<- ggplot(annot2,aes(x=Immune.Subtype, y= LNG, fill=Immune.Subtype))+
      geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
      geom_violin(
        aes(fill=NULL, col=Immune.Subtype),
        alpha=0.4, 
        show.legend = F)+
      geom_jitter(size=0.2, width = 0.1)+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      scale_color_manual(values=c("C1"="red","C2"="yellow","C3"="green3",
                                  "C4"="cyan","C5"="blue","C6"="pink"))+
      scale_fill_manual(values=c("C1"="red","C2"="yellow","C3"="green3",
                                 "C4"="cyan","C5"="blue","C6"="pink"))+
                       # labels = c("C1 (n=42)", "C2 (n=46)", "C3 (n=31)",
                       #             "C4 (n=10)", "C6 (n=6)"))+ 
      labs(y=paste0("log2 ",varn,"\ngene expression"), 
           x= "Immune Subtype", 
           title = subtype[x]) +
      theme(legend.position = "bottom", 
            legend.title = element_text(face="bold", size=8),
            legend.key=element_rect(size=10, color=NA), 
            legend.key.size=unit(5,"mm"),
            legend.text=element_text(size=8), 
            legend.direction = "vertical",
            legend.box = "vertical" )+
      stat_compare_means(comparisons = comparison,
                         size=3)+
      #scale_y_continuous(expand = expansion(mult = c(0, 0.3)))+
      guides(fill=guide_legend(ncol= 3), 
             col=guide_legend(ncol=3, title = "Immune\nSubtype"))
    
    g2
  }
  
  if(subtype[x]=="BRCA.Normal" | subtype[x]=="BRCA.LumA"  ){
    # Define comparison to plot
    if(subtype[x]=="BRCA.LumA"  ){
      comparison <- list(c("C1","C2"), c("C2","C3"), c("C3","C4"), 
                         c("C2","C4"), c("C1","C4"), c("C4","C6"), c("C2","C6") )
    }
    if(subtype[x]=="BRCA.Normal"  ){
      comparison <- list(c("C2","C3"),c("C2","C4"), c("C1","C3"), c("C3","C4"), 
                          c("C1","C4"), c("C4","C6"), c("C3","C6") )
    }
    annot2<- annot2[annot2$Immune.Subtype %in% c("C1","C2", "C3","C4","C6"),]
    ## Plotting
    library(ggpubr)
    ymax <- max(annot2[,varn])
    g2<- ggplot(annot2,aes(x=Immune.Subtype, y= LNG, fill=Immune.Subtype))+
      geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
      geom_violin(
        aes(fill=NULL, col=Immune.Subtype),
        alpha=0.4, 
        show.legend = F)+
      geom_jitter(size=0.2, width = 0.1)+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      scale_color_manual(values=c("C1"="red","C2"="yellow","C3"="green3",
                                  "C4"="cyan","C5"="blue","C6"="pink"))+
      scale_fill_manual(values=c("C1"="red","C2"="yellow","C3"="green3",
                                 "C4"="cyan","C5"="blue","C6"="pink"))+
                        #labels = c("C1 (n=42)", "C2 (n=46)", "C3 (n=31)",
                        #           "C4 (n=10)", "C6 (n=6)"))+
      labs(y=paste0("log2 ",varn,"\ngene expression"), 
           x= "Immune Subtype", 
           title = subtype[x]) +
      theme(legend.position = "bottom", 
            legend.title = element_text(face="bold", size=8),
            legend.key=element_rect(size=10, color=NA), 
            legend.key.size=unit(5,"mm"),
            legend.text=element_text(size=8), 
            legend.direction = "vertical",
            legend.box = "vertical" )+
      stat_compare_means(label.x = c(1), 
                        label.y = ymax*1.7,
                        label="p.format",
                         size=3)+
      stat_compare_means(comparisons = comparison,
                         size=3)+
      guides(fill=guide_legend(ncol= 3), 
             col=guide_legend(ncol=3, title = "Immune\nSubtype"))
    
    g2
  }
  
  if(subtype[x]=="BRCA.LumB"  ){
    annot2<- annot2[annot2$Immune.Subtype %in% c("C1","C2", "C3","C4"),]
    # Define comparison to plot
    comparison <- list(c("C2","C3"), c("C1","C4"), 
                       c("C1","C3") ,c("C2","C4")  )
    
      ## Plotting
    library(ggpubr)
    ymax <- max(annot2[,varn])
    g2<- ggplot(annot2,aes(x=Immune.Subtype, y= LNG, fill=Immune.Subtype))+
      geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
      geom_violin(
        aes(fill=NULL, col=Immune.Subtype),
        alpha=0.4, 
        show.legend = F)+
      geom_jitter(size=0.2, width = 0.1)+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      scale_color_manual(values=c("C1"="red","C2"="yellow","C3"="green3",
                                  "C4"="cyan","C5"="blue","C6"="pink"))+
      scale_fill_manual(values=c("C1"="red","C2"="yellow","C3"="green3",
                                 "C4"="cyan","C5"="blue","C6"="pink"))+
                       # labels = c("C1 (n=59)", "C2 (n=87)", "C3 (n=7)",
                       #             "C4 (n=28)", "C6 (n=0)"))+
      labs(y=paste0("log2 ",varn,"\ngene expression"), 
           x= "Immune Subtype", 
           title = subtype[x]) +
      theme(legend.position = "bottom", 
            legend.title = element_text(face="bold", size=8),
            legend.key=element_rect(size=10, color=NA), 
            legend.key.size=unit(5,"mm"),
            legend.text=element_text(size=8), 
            legend.direction = "vertical",
            legend.box = "vertical" )+
      stat_compare_means(label.x = c(1), 
                         label.y = ymax*1.3,
                         label="p.format",
                         size=3)+
      stat_compare_means(comparisons = comparison,
                         size=3)+
      guides(fill=guide_legend(ncol= 3), 
             col=guide_legend(ncol=3, title = "Immune\nSubtype"))
    
    g2
  }
  ggsave(plot= g2, file= paste0(subtype[x],"_boxplot_immune_subtype.svg"), device="svg", 
      dpi=300, height = 120, width=100, units = "mm")

#####################################################################################
#Kaplan Meier analysis and plot
#####################################################################################

  #For PFI
fit <- survfit(Surv(time=PFI.Time, event = PFI)~Actv, data=annot)
g1 <- ggsurvplot(fit,
                 data=annot,
                 title= paste0("Progression Free Interval - ", varn),
                 xlab="Time in days",
                 legend.title="",
                 legend.labs=c(paste0("lncRNA High (n=",summary(annot$Actv)[1],")"),
                               paste0("lncRNA Low (n=",summary(annot$Actv)[2],")")),
                 conf.int = T, 
                 conf.int.style= "ribbon",
                 conf.int.alpha=0.15,
                 pval = T,
                 pval.size=3, 
                 risk.table = TRUE,
                 risk.table.y.text=F, 
                 risk.table.fontsize=3, 
                 risk.table.title="",
                 ggtheme = theme_grey() )
g1$plot <- ggpar(g1$plot, font.x = 0, font.title = 10)
#g1$table <- ggpar(g1$table, font.y = 0)

# For OS
fit <- survfit(Surv(time=OS.Time, event = OS)~Actv, data=annot)  
g0 <- ggsurvplot(fit,
                 data=annot,
                 title= paste0("Overall Survival - ", varn),
                 xlab="Time in days",
                 legend.title="",
                 legend.labs=c(paste0("lncRNA High (n=",summary(annot$Actv)[1],")"),
                               paste0("lncRNA Low (n=",summary(annot$Actv)[2],")")),
                 conf.int = T, 
                 conf.int.style= "ribbon",
                 conf.int.alpha=0.15,
                 pval = T,
                 pval.size=3, 
                 risk.table = TRUE,
                 risk.table.y.text=F, 
                 risk.table.fontsize=3, 
                 risk.table.title="",
                 ggtheme = theme_grey() )
  g0$plot <- ggpar(g0$plot, font.x = 0, font.title = 10)
  #g0$table <- ggpar(g0$table, font.y = 0)

  #Plotting
row1<- plot_grid(g1$plot, 
                 g1$table, 
                 ncol = 1, 
                 rel_heights = c(1,0.3), 
                 align="v")
plot_grid(row1)
row2<- plot_grid(g0$plot, 
                 g0$table, 
                 ncol = 1, 
                 rel_heights = c(1,0.3), 
                 align="v")
plot_grid(row2)
g <- plot_grid(row2, row1, nrow=2)
# Saving results
ggsave(plot = g,file= paste0(subtype[x],"_KM.svg"), device="svg", 
       dpi=300, height = 240, width=100, units = "mm")

#####################################################################################
# Boxplot in molecular subtypes
#####################################################################################

type_nao <- c(subtype[-x],"Non-tumoral")
ymax <- max(gexp_lnc[,varn])
comparison<- list(NULL)
for(i in 1:5){
  comparison[[i]] <- c(subtype[x], type_nao[i])
}
g<- ggplot(gexp_lnc,aes(x=TCGA.Subtype, y= gexp_lnc[,varn], fill=TCGA.Subtype))+
  geom_boxplot(width=c(0.5), lwd=0.2, outlier.size = 0.5)+  
  geom_violin(
    aes(fill=NULL, col=TCGA.Subtype),
    alpha=0.4, 
    show.legend = F)+
  geom_jitter(size=0.2, width = 0.1)+
  scale_color_manual(values=c("BRCA.Basal"="sienna1", "BRCA.Normal"="purple3",
                              "BRCA.Her2" ="hotpink", "BRCA.LumA"="cornflowerblue",
                              "BRCA.LumB"="springgreen2", "Non-tumoral"="red"))+
  scale_fill_manual(values=c("BRCA.Basal"="sienna1", "BRCA.Normal"="purple3",
                             "BRCA.Her2" ="hotpink", "BRCA.LumA"="cornflowerblue",
                             "BRCA.LumB"="springgreen2", "Non-tumoral"="red"),
                    labels = c("Basal (n=169)","LumA (n=499)" ,
                               "LumB (n=184)" ,"Normal (n=136)",
                               "Her2 (n=72)", "Non-tumoral (n=113)"))+
  scale_x_discrete(labels=c("Basal",  "LumA","LumB","Normal", "Her2","Non-tumoral" ))+
  labs(y=paste0("log2 ",varn,"\ngene expression"), 
       x= "TCGA Subtype"
       #title = subtype[x]
  ) +
  #coord_cartesian(ylim=c(0,1.3))+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=8),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(5,"mm"),
        legend.text=element_text(size=8), 
        legend.direction = "vertical",
        legend.box = "vertical" )+
  stat_compare_means(label.x = c(1), 
                     label.y = ymax*1.7,
                     size=3)+
  stat_compare_means(comparisons = comparison,
                     label="p.format",
                     size=3)+
  guides(fill=guide_legend(nrow = 2), 
         col=guide_legend(nrow=2, title = "TCGA\nSubtype"))

g
ggsave(plot=g, file= paste0(subtype[x],"_boxplot_molecular_subtype.svg"), device = "svg",
       dpi=300, width = 110, height = 120 , units="mm")
}
```

R Session
----
```r
sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.2 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#  [1] LC_CTYPE=pt_BR.UTF-8       LC_NUMERIC=C               LC_TIME=pt_BR.UTF-8       
#  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=pt_BR.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=pt_BR.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=pt_BR.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] msigdbr_7.2.1        fgsea_1.16.0         cowplot_1.1.1        survminer_0.4.9      ggpubr_0.4.0        
#  [6] survival_3.2-11      ComplexHeatmap_2.6.2 forcats_0.5.1        stringr_1.4.0        dplyr_1.0.5         
# [11] purrr_0.3.4          readr_1.4.0          tidyr_1.1.3          tibble_3.1.1         ggplot2_3.3.3       
# [16] tidyverse_1.3.1     
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-0    ggsignif_0.6.1      rjson_0.2.20        ellipsis_0.3.2      rio_0.5.26         
#   [6] circlize_0.4.12     markdown_1.1        GlobalOptions_0.1.2 fs_1.5.0            gridtext_0.1.4     
#  [11] ggtext_0.1.1        clue_0.3-59         rstudioapi_0.13     farver_2.1.0        fansi_0.4.2        
#  [16] lubridate_1.7.10    xml2_1.3.2          splines_4.0.5       knitr_1.33          ggunchained_0.0.1  
#  [21] jsonlite_1.7.2      Cairo_1.5-12.2      broom_0.7.6         km.ci_0.5-2         cluster_2.1.2      
#  [26] dbplyr_2.1.1        png_0.1-7           compiler_4.0.5      httr_1.4.2          backports_1.2.1    
#  [31] assertthat_0.2.1    Matrix_1.3-2        cli_2.5.0           htmltools_0.5.1.1   tools_4.0.5        
#  [36] igraph_1.2.6        gtable_0.3.0        glue_1.4.2          ggthemes_4.2.4      fastmatch_1.1-0    
#  [41] Rcpp_1.0.6          carData_3.0-4       cellranger_1.1.0    vctrs_0.3.8         writexl_1.4.0      
#  [46] svglite_2.0.0       xfun_0.22           networkD3_0.4       openxlsx_4.2.3      rvest_1.0.0        
#  [51] lifecycle_1.0.0     rstatix_0.7.0       zoo_1.8-9           scales_1.1.1        hms_1.0.0          
#  [56] parallel_4.0.5      RColorBrewer_1.1-2  yaml_2.2.1          curl_4.3            gridExtra_2.3      
#  [61] KMsurv_0.1-5        stringi_1.5.3       S4Vectors_0.28.1    BiocGenerics_0.36.1 zip_2.1.1          
#  [66] BiocParallel_1.24.1 shape_1.4.5         rlang_0.4.10        pkgconfig_2.0.3     systemfonts_1.0.1  
#  [71] matrixStats_0.58.0  evaluate_0.14       lattice_0.20-41     htmlwidgets_1.5.3   labeling_0.4.2     
#  [76] tidyselect_1.1.1    DataExplorer_0.8.2  magrittr_2.0.1      R6_2.5.0            IRanges_2.24.1     
#  [81] generics_0.1.0      DBI_1.1.1           pillar_1.6.0        haven_2.4.1         foreign_0.8-81     
#  [86] withr_2.4.2         abind_1.4-5         modelr_0.1.8        crayon_1.4.1        car_3.0-10         
#  [91] survMisc_0.5.5      utf8_1.2.1          rmarkdown_2.7       GetoptLong_1.0.5    readxl_1.3.1       
#  [96] data.table_1.14.0   reprex_2.0.0        digest_0.6.27       xtable_1.8-4        stats4_4.0.5       
# [101] munsell_0.5.0      
```

