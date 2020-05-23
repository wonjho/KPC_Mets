#### RNASEQ ANALYSIS OF SITE-SPECIFIC TME FEATURES FOR PDAC
#### LAsT UPDATED 05/19/2020

######################################################
#### PREPARE WORKSPACE ####

rm(list = ls())
setwd("C:/Users/who10/Desktop/Research/KPC Immunophenotyping/RNAseq")

knitr::opts_chunk$set(echo = F, message = FALSE,
                      warning = FALSE, cache = T)

#all libraries needed
{
library("digest")
library("rmarkdown")
library(simpleaffy)
library("DT")
library("limma")
library(DESeq2)
library(edgeR)
library(ggplot2)
library("org.Mm.eg.db")
library('GeneOverlap')
library('GO.db')
library('sva')
library('KEGG.db')
library('ComplexHeatmap')
library('pca3d')
library('CoGAPS')
library('GSA')
library('biomaRt')
library('gtools')
library('ggrepel')
library('pheatmap')
library('RColorBrewer')
}
  
sessionInfo()

getCols=function(class,subclass=NULL,bahman=FALSE){ # make it also work if there is only one class or subclass
  # set the six hues
  mycols=c(0,120/360,200/360,60/360,300/360)
  # if there are subclasses
  if(!is.null(subclass)){
    # if the subclasses are a numeric range
    if(is.numeric(subclass)){
      
      # get the range to be 0-1
      testrange2=(subclass-min(subclass))/max(subclass-min(subclass))
      #set the color range to be 0.4 to 1, as darker colors can look black despite the hue
      temp=0.4+testrange2*0.6
      
      # get indexes for classes of samples
      uclas=unique(class)
      whichclass=class
      for(i in 1:length(class)){
        whichclass[which(class==uclas[i])]=i
      }
      
      # assign the colors to the classes vector
      classcols=mycols[as.numeric(whichclass)]
      ColResult=class
      
      # get the color codes into the vector
      ColResult=hsv(h=classcols,s=1,v=temp)
      
      # if subclass is not numeric
    }else{
      
      # create the list for storing the color codes for each class, relative to the number of subclasses
      temp=list()
      
      # loop through the unique classes, for each sample put in the appropriate range of colors
      for(i in 1:length(unique(class))){
        n=length(unique(subclass[which(class==unique(class)[i])])) # how many subclasses are in this class
        
        if(n==1){
          temp[[i]]=hsv(h=mycols[i],s=1,v=1)
        }else{
          temp[[i]]=hsv(h=mycols[i],s=1,v=seq(1,0.3,-0.7/(n-1))) # for that list element, create the color range
        }
        
      }
      
      # will need to get the numeric rendition of which class and subclass each sample is
      whichsub=subclass
      whichclass=class
      
      uclas=unique(class)
      
      for(i in 1:length(class)){
        # which samples are each of the unique classes
        whichclass[which(class==uclas[i])]=i
        # the unique subclasses for each class
        usub=unique(subclass[which(class==unique(class)[i])])
        
        for(j in 1:length(usub))
          # which samples are each of the unique subclasses
          whichsub[which(subclass==usub[j])]=j
      }
      
      whichclass=as.numeric(whichclass)
      whichsub=as.numeric(whichsub)
      
      ColResult=class
      for(i in 1:length(class)){
        ColResult[i]=temp[[whichclass[i]]][whichsub[i]]
      }
      
    }
  }else{ #if there is no subclass, a rainbow is sufficient
    mycols=rainbow(length(unique(class)))
    uclas=unique(class)
    whichclass=class
    for(i in 1:length(class)){
      whichclass[which(class==uclas[i])]=i
    }
    ColResult=mycols[as.numeric(whichclass)]
  }
  if(bahman==TRUE){
    bahmanlist=as.list(unique(class))
    names(bahmanlist)=unique(class)
    names(ColResult)=subclass
    for(i in class){
      dup=duplicated(ColResult[class==i])
      bahmanlist[[i]]=ColResult[class==i][which(!dup)]
    }
    ColResult=bahmanlist
  }
  return(ColResult)
} #end getCols

######################################################
#### RETRIEVE ANNOTATIONS AND COUNTS DATA ####

anno=read.table("Won_RNAseq_sampleAnnot.txt",sep="\t",header=T,as.is=T)
rownames(anno)=anno$sampNames
anno$Design <- apply(anno[,c('type','tissue')],1,paste,collapse=' ')

countsDat <- readRDS("txi_KPC_RNAseq_data_prot.rds")
storage.mode(countsDat[["counts"]]) <- "integer"
countsDat <- as.data.frame(countsDat[["counts"]])

STCDataSet <- DESeqDataSetFromMatrix(countData = countsDat, colData = anno[colnames(countsDat),], 
                                     design = ~0+Design) 
STCDataSet <- STCDataSet[rowSums(counts(STCDataSet))>1,]
dat=countsDat 

logSTCDataSet <- rlog(STCDataSet) #log transform counts data
mat <- assay(logSTCDataSet) #create a SummarizedExperiment object

#boxplot(log2(dat+1),las=2, ylab='log transformed read counts')
#title('before rlog transformation')
#boxplot(mat,las=2, ylab='rlog transformed read counts')
#title('after rlog transformation')

#Renaming columns/sample names and ordering them by what I want
newcolnames<-c("KPC1","NLLung2","NLLung3","NLLung4","NLLung5",
               "NLLiver1","NLLiver2","NLLiver3","NLLiver4","NLLiver5",
               "KPCLung1","KPC2","KPCLung2","KPCLung3","KPCLung4",
               "KPCLung5","KPCLiver1","KPCLiver2","KPCLiver3","KPCLiver4",
               "KPCLiver5","KPCLiver6","KPC3","KPCLiver7","KPC4",
               "KPC5","KPC6","KPC7","KPC8","NLLung1")
colnames(mat)<-newcolnames
mat <- mat[,c("KPC1","KPC2","KPC3","KPC4","KPC5","KPC6","KPC7","KPC8",
                    "NLLung1","NLLung2","NLLung3","NLLung4","NLLung5",
                    "NLLiver1","NLLiver2","NLLiver3","NLLiver4","NLLiver5",
                    "KPCLung1","KPCLung2","KPCLung3","KPCLung4","KPCLung5",
                    "KPCLiver1","KPCLiver2","KPCLiver3","KPCLiver4","KPCLiver5","KPCLiver6","KPCLiver7")]

#Exploratory analysis for technical artifacts

pcaplot <- plotPCA(logSTCDataSet, intgroup = c("type", "tissue"))
#saveRDS(file="pcaplot.rds",pcaplot)
#pcaplot2<- readRDS("pcaplot.rds")
#rm(pcaplot2)

pdf("PCAplot.pdf",width=4,height=4)
pcaplot + 
  geom_text_repel(aes(label=newcolnames), show.legend = FALSE, size=2.5) +
  scale_color_manual(values=c("black","red","blue","magenta","green"))+
  theme(legend.position = "none")
dev.off()

plot(standard.pearson(assay(logSTCDataSet)))

anno=anno[colnames(countsDat),]

######################################################
#### Run DESeq ####

dds <- DESeq(STCDataSet, betaPrior = F)
resultsNames(dds)

# Get the results for KPC lung vs KPC cells
res2DvLG <- results(dds,contrast=c(1,0,-1,0,0))
sig2DvLG <- row.names(res2DvLG)[which(abs(res2DvLG$log2FoldChange) > 1 & res2DvLG$padj < 0.05)]

# Get the results for KPC liver vs KPC cells
res2DvLV <- results(dds,contrast=c(1,-1,0,0,0))
sig2DvLV <- row.names(res2DvLV)[which(abs(res2DvLV$log2FoldChange) > 1 & res2DvLV$padj < 0.05)]

# Get the results for KPC lung vs. liver
resLGvLV <- results(dds,contrast=c(0,1,-1,0,0))
sigLGvLV <- row.names(resLGvLV)[which(abs(resLGvLV$log2FoldChange) > 1 & resLGvLV$padj < 0.05)]
padjLGvLV <- resLGvLV$padj;names(padjLGvLV)<-rownames(resLGvLV)#need this for labeling heatmaps with significance later


# Visualize and print out the genes that are different across groups by Venn diagram
allGene <- unique(rownames(resLGvLV))
allStats <- matrix(0, nrow=length(allGene), ncol=3,
                   dimnames = list(allGene,c('2D - Lung', 
                                             '2D - Liver','Lung - Liver')))
allStats[sig2DvLG,'2D - Lung'] <- 1
allStats[sig2DvLV,'2D - Liver'] <- 1
allStats[sigLGvLV,'Lung - Liver'] <- 1
vennDiagram(allStats) #plot
allStats2=data.frame(allStats)
allStats2$sum = rowSums(allStats)
write.table(cbind(GeneSymbol=row.names(allStats), allStats), #print out all the genes
            file=sprintf('DEGenes%s.csv',Sys.Date()), sep=",", row.names = F) 
anysig=rownames(allStats2)[which(allStats2$sum>0)]   

######################################################
#### PATHWAY ANALYSIS ####

#KEGG gene sets

SYMBOLToKEGG <- mapIds(org.Mm.eg.db,keys=rownames(mat),column='PATH',keytype = 'SYMBOL',multiVals = list)
KEGGToSYMBOL <- sapply(reverseSplit(SYMBOLToKEGG),unique)
KEGGToSYMBOL <- KEGGToSYMBOL[sapply(KEGGToSYMBOL,length)>5]

stat2DvLG <- res2DvLG$stat
names(stat2DvLG) <- row.names(res2DvLG)
KEGGStat2DvLG <- p.adjust(sapply(KEGGToSYMBOL,function(x){wilcoxGST(index = x, statistics = stat2DvLG,
                                                                  alternative='either')}),
                        method='BH')


stat2DvLV <- res2DvLV$stat
names(stat2DvLV) <- row.names(res2DvLV)
KEGGStat2DvLV <- p.adjust(sapply(KEGGToSYMBOL,function(x){wilcoxGST(index = x, statistics = stat2DvLV,
                                                                    alternative='either')}),
                          method='BH')


statLGvLV <- resLGvLV$stat
names(statLGvLV) <- row.names(resLGvLV)
KEGGStatLGvLV <- p.adjust(sapply(KEGGToSYMBOL,function(x){wilcoxGST(index = x, statistics = statLGvLV,
                                                                    alternative='either')}),
                          method='BH')


KEGGTable <- data.frame(Lungvs2D=KEGGStat2DvLG, 
                        Livervs2D=KEGGStat2DvLV,
                        LungvsLiver=KEGGStatLGvLV)

KEGGTable <- cbind(KEGG.Name=unlist(as.list(KEGGPATHID2NAME)[row.names(KEGGTable)]),
                   KEGGTable)

saveWidget(datatable(KEGGTable),"KEGGresult.html")

#KEGG pathways

KEGGSignaling <- KEGGToSYMBOL[row.names(KEGGTable)[KEGGTable$KEGG.Name %in% 
                                                     grep('signaling pathway',KEGGTable$KEGG.Name,value=T)]]


stat2DvLG <- res2DvLG$stat
names(stat2DvLG) <- row.names(res2DvLG)
KEGGStat2DvLG <- p.adjust(sapply(KEGGSignaling,function(x){wilcoxGST(index = x, statistics = stat2DvLG,
                                                                    alternative='either')}),
                          method='BH')


stat2DvLV <- res2DvLV$stat
names(stat2DvLV) <- row.names(res2DvLV)
KEGGStat2DvLV <- p.adjust(sapply(KEGGSignaling,function(x){wilcoxGST(index = x, statistics = stat2DvLV,
                                                                    alternative='either')}),
                          method='BH')


statLGvLV <- resLGvLV$stat
names(statLGvLV) <- row.names(resLGvLV)
KEGGStatLGvLV <- p.adjust(sapply(KEGGSignaling,function(x){wilcoxGST(index = x, statistics = statLGvLV,
                                                                    alternative='either')}),
                          method='BH')


KEGGTable <- data.frame(Lungvs2D=KEGGStat2DvLG, 
                        Livervs2D=KEGGStat2DvLV,
                        LungvsLiver=KEGGStatLGvLV)

KEGGTable <- cbind(KEGG.Name=unlist(as.list(KEGGPATHID2NAME)[row.names(KEGGTable)]),
                   KEGGTable)

saveWidget(datatable(KEGGTable),"KEGGresult2.html")



######################################################
#### HEATMAP PLOTS ####

#PLOTTER FUNCTIONS

#KEGG gene set heatmaps: Ver1 will allow me to plot only the genes that are specific to one KEGG set

selectgenestoplot1 <- function (matrixcounts, keggcode){
  
  matkeggs <- matrixcounts[which(rownames(matrixcounts) %in% KEGGToSYMBOL[[keggcode]]),]

return(matkeggs)

}

#KEGG gene set exploratory heatmaps: Ver2 will then select genes that are of highest or lowest mean expression in either KPC liver or KPC lung 
#this allows me to focus on genes that are not inadvertently demonstrating a "mixed cell" effect
#this however does not help me to avoid potential non-tumor "cell enrichment" effect

selectgenestoplot2 <- function (matrixcounts, keggcode){
  
  kegggenelist <- matrixcounts[which(rownames(matrixcounts) %in% KEGGToSYMBOL[[keggcode]]),]
  
  KPC<-rowMeans(kegggenelist[,1:8])
  NLLung<-rowMeans(kegggenelist[,9:13])
  NLLiver<-rowMeans(kegggenelist[,14:18])
  KPCLung<-rowMeans(kegggenelist[,19:23])
  KPCLiver<-rowMeans(kegggenelist[,24:30])
  
  KPCLungmax <-
  (KPCLung > KPCLiver) & 
  (KPCLung > KPC) & 
  (KPCLung > NLLung) & 
  (KPCLung > NLLiver)
  KPCLungmin <-
  (KPCLung < KPCLiver) & 
  (KPCLung < KPC) & 
  (KPCLung < NLLung) & 
  (KPCLung < NLLiver)
  KPCLivermax <- 
  (KPCLiver > KPCLung) & 
  (KPCLiver > KPC) & 
  (KPCLiver > NLLung) & 
  (KPCLiver > NLLiver)
  KPCLivermin <- 
  (KPCLiver < KPCLung) & 
  (KPCLiver < KPC) & 
  (KPCLiver < NLLung) & 
  (KPCLiver < NLLiver)

matkeggs <- kegggenelist[(KPCLivermax | KPCLivermin | KPCLungmax | KPCLungmin),]

return(matkeggs)

}

#Function for generating the pdfs for particular KEGG gene set and code

pheatmapspdf <- function (mat=mat, name, code, rowfontsize=8, colfontsize=8){
  pdf(paste0("kegg",name,"all.pdf"), width=6.5, height=11)
  pheatmap(selectgenestoplot1(mat,code),
           main=paste0("All Genes in ",name," KEGG"),
           fontsize_row = rowfontsize,
           fontsize_col = colfontsize)
  dev.off()
  pdf(paste0("kegg",name,"select.pdf"), width=6, height=6)
  pheatmap(selectgenestoplot2(mat,code),
           main=paste0("Select Genes in ",name," KEGG"),
           fontsize_row = rowfontsize,
           fontsize_col = colfontsize,
           border_color=NA,
           cluster_cols =TRUE)
  dev.off()
}

pheatmapspdf(mat,"TP53", "04115")
pheatmapspdf(mat,"PPAR", "03320")
pheatmapspdf(mat,"WNT", "04310")
pheatmapspdf(mat,"ErbB", "04012", rowfontsize = 4)
pheatmapspdf(mat,"Neurotropin", "04722")
pheatmapspdf(mat,"Pancreatic Cancer", "05212")
pheatmapspdf(mat,"Cytokine", "04060",rowfontsize = 4)
pheatmapspdf(mat,"Adhesion", "04514")
pheatmapspdf(mat,"HepC", "05160", rowfontsize = 4)

#Heatmaps for manually constructed list of genes
#Ver3 allows for plotting all 5 groups as means (one column per group)
selectgenestoplot3 <- function(mat=mat,manuallist){
  matmanual <- mat[which(rownames(mat) %in% manuallist),]
  KPC<-rowMeans(matmanual[,1:8])
  NLLung<-rowMeans(matmanual[,9:13])
  NLLiver<-rowMeans(matmanual[,14:18])
  KPCLung<-rowMeans(matmanual[,19:23])
  KPCLiver<-rowMeans(matmanual[,24:30])
  matmanualmean <- cbind(NLLung,NLLiver,KPC,KPCLung,KPCLiver)
  return(matmanualmean)
} 
#Ver4 allows for plotting of only KPC, KPCLG, KPCLV
selectgenestoplot4 <- function(mat=mat,manuallist){
  matmanual <- mat[which(rownames(mat) %in% manuallist),]
  KPC<-rowMeans(matmanual[,1:8])
  KPCLung<-rowMeans(matmanual[,19:23])
  KPCLiver<-rowMeans(matmanual[,24:30])
  matmanualmean <- cbind(KPC,KPCLung,KPCLiver)
  return(matmanualmean)
} 
#Ver5 allows for plotting of all samples
selectgenestoplot5 <- function(mat=mat,manuallist){
  matmanual <- mat[which(rownames(mat) %in% manuallist),]
  return(matmanual)
}
#Function for annotating statistical significance between KPCLG and KPCLV on the gene names
phmlabelingstat <- function(res_hmap,genelist){
  phmlabels<- genelist[order(genelist)][res_hmap$tree_row$order] #pheatmap function provides an output of tree_row$order
  padjLGvLV[phmlabels][is.na(padjLGvLV[phmlabels])]<-1
  sig1<-(padjLGvLV[phmlabels]<=0.05)&(padjLGvLV[phmlabels]>0.01)
  sig2<-(padjLGvLV[phmlabels]<=0.01)&(padjLGvLV[phmlabels]>0.005)
  sig3<-(padjLGvLV[phmlabels]<=0.005)
  phmlabels[sig1]<-paste0(phmlabels[sig1],"*")
  phmlabels[sig2]<-paste0(phmlabels[sig2],"**")
  phmlabels[sig3]<-paste0(phmlabels[sig3],"***")
  annphmlabels<-sort(phmlabels[order(phmlabels)]) #pheatmap stores a sort ordered vector
  return(annphmlabels)
} 


#PLOTS

#For genes that allow us to see the end result of processing of samples
Controlgenes<-c("Epb42","Ptprc","Pdpn","Vim","Epcam","B2m","Pecam1","Pdx1","Krt18", "Krt19", "Sox9", "Muc1", "Muc4", "Muc5ac", "Cdh1", "Cd44", "F2", "Fgb", "Sftpa1", "Sftpb","Sema5a",
                "Ptprc","Cd3e","Cd19","Klrb1","Itgam", "Itgax")
pdf("EnrichedProfileCtrlSet.pdf",height=3.5,width=4)
pheatmap(selectgenestoplot5(mat,Controlgenes), 
         color = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100),
         cluster_cols = FALSE,
         scale = "row",
         border_color = NA,
         main="Expression Profiles of Enriched Samples",
         fontsize_row = 7, fontsize_col = 5,
         legend=TRUE)
dev.off()


Chemokines <- c("Ccl1","Ccl2","Ccl3","Ccl4", "Ccl5","Ccl6", "Ccl7","Ccl8","Ccl9","Ccl11", "Ccl17", "Ccl20", "Ccl21a", "Ccl22", "Ccl24", "Ccl25", "Ccl27a", "Ccl28",
                "Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcl10", "Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl17",
                "Cx3cl1",
                "Xcl1")
chemok<-pheatmap(selectgenestoplot3(mat,Chemokines), 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         main="Chemokine Profiles",
         fontsize_row = 7, fontsize_col = 9,
         labels_row = Chemokines[order(Chemokines)],
         legend=TRUE)
dev.off()
Chemokineslab<-phmlabelingstat(chemok,Chemokines)
pdf("EnrichedChemokines.pdf",height=4,width=4)
pheatmap(selectgenestoplot5(mat,Chemokines), 
         color = rev(colorRampPalette(brewer.pal(n = 11, name = "PuOr"))(100)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         scale="row",
         main="Chemokine Profiles",
         fontsize_row = 7, fontsize_col = 5,
         labels_row = sort(Chemokineslab[order(Chemokines)]),
         legend=TRUE)
dev.off()


Immsupp <- c("Tgfb1","Tgfb2","Tgfb3","Il10"
             , "Cd274","Pdcd1lg2", "Fgl1"
             , "Ido1","Csf2", "Fasl" 
             #, "H2-K1","H2-D1","H2-M1","H2-Q1"
             )
immsup <- pheatmap(selectgenestoplot3(mat,Immsupp), 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         main="Tumor Derived Immune Suppressors",
         fontsize_row = 7, fontsize_col = 9,
         #labels_row = sort(Chemokineslab[order(Chemokines)]),
         legend=TRUE)
dev.off()
Immsupplab<-phmlabelingstat(immsup,Immsupp)
pdf("EnrichedImmsupp.pdf",height=2,width=4)
pheatmap(selectgenestoplot5(mat,Immsupp), 
         color = rev(colorRampPalette(brewer.pal(n = 11, name = "PuOr"))(100)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         scale = "row",
         main="Immune Suppressors in Enriched Samples",
         fontsize_row = 7, fontsize_col = 5,
         labels_row = sort(Immsupplab[order(Chemokines)]),
         cellwidth = 5, cellheight = 6,
         legend=TRUE)
dev.off()

#Need to load PD1 pathway here
{
PDL1pathway<-c("Hif1a",
              "Egf",
              "Egfr",
              "Hras",
              "Kras",
              "Nras",
              "Raf1",
              "Map2k1",
              "Map2k2",
              "Mapk1",
              "Mapk3",
              "Fos",
              "Jun",
              "Eml4",
              "Alk",
              "Pik3r1",
              "Pik3r2",
              "Pik3r3",
              "Pik3ca",
              "Pik3cd",
              "Pik3cb",
              "Pten",
              "Akt1",
              "Akt2",
              "Akt3",
              "Mtor",
              "Rps6kb1",
              "Rps6kb2",
              "Chuk",
              "Ikbkb",
              "Ikbkg",
              "Nfkbia",
              "Nfkbib",
              "Nfkbie",
              "Nfkb1",
              "Rela",
              "Ifng",
              "Ifngr1",
              "Ifngr2",
              "Jak1",
              "Jak2",
              "Stat1",
              "Stat3",
              "Tlr2",
              "Tlr4",
              "Tlr9",
              "Tirap",
              "Myd88",
              "Traf6",
              "Nfatc1",
              "Nfatc2",
              "Nfatc3",
              "Ticam1",
              "Ticam2",
              "Cd274",
              "Pdcd1lg2")
}


pdl1tumor <- pheatmap(selectgenestoplot3(mat,PDL1pathway), 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         main="PDL1 Pathway",
         fontsize_row = 7, fontsize_col = 9,
         #labels_row = sort(Chemokineslab[order(Chemokines)]),
         legend=TRUE)
dev.off()
PDL1tumorpathlab<-phmlabelingstat(pdl1tumor,PDL1pathway)
pdf("EnrichedPDL1path.pdf",height=5.5,width=4)
pheatmap(selectgenestoplot5(mat,PDL1pathway), 
         color = rev(colorRampPalette(brewer.pal(n = 11, name = "PuOr"))(100)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         scale = "row",
         main="PDL1 Pathway in Enriched Samples",
         fontsize_row = 7, fontsize_col = 5,
         labels_row = sort(PDL1tumorpathlab[order(PDL1pathway)]),
         legend=TRUE)
dev.off()


#Need to load Fgl1pathway here
{
  Fgl1pathway<-c(
"Fgg",
"Fgb",
"Fgl1",
"Cd83",
"Cxcr5",
"Itih4",
"Cd79b",
"Vav3",
"Bach2",
"Vpreb1",
"Vpreb2",
"Pirb",
"Crlf2",
"Cxcl15",
"Clec10a",
"Capg",
"Mbl1",
"Cdh17",
"Ptch2",
"Slmap",
"Epb41",
"Cd40lg",
"Galr2",
"Lag3",
"Ccr1",
"Fcgr3",
"Il7",
"Tnfrsf9",
"Il7r",
"Stat5a",
"Cenpj",
"Stat5b",
"Dapp1",
"Rac2",
"Inpp5d",
"Cd69",
"Cd274",
"Klhl6",
"Rassf5",
"Stk17b",
"Ifnar1",
"Tap2",
"Grap2",
"Sit1",
"Tcf7",
"Cd8b1",
"Elk4"
)
}
fgl1path <- pheatmap(selectgenestoplot3(mat,Fgl1pathway), 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         main="Fgl1 Pathway",
         fontsize_row = 7, fontsize_col = 9,
         legend=TRUE)
dev.off()
Fgl1pathlab<-phmlabelingstat(fgl1path,Fgl1pathway)
pdf("EnrichedFgl1path.pdf",height=5.5,width=4)
pheatmap(selectgenestoplot5(mat,Fgl1pathway), 
         color = rev(colorRampPalette(brewer.pal(n = 11, name = "PuOr"))(100)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         scale = "row",
         main="FGL1 Network in Enriched Samples",
         fontsize_row = 7, fontsize_col = 5,
         labels_row = sort(Fgl1pathlab[order(Fgl1pathway)]),
         legend=TRUE)
dev.off()

#to plot the statistical results of Wilcoxon GST of KEGG signaling pathways
pathwayres<- data.frame(
  KEGGids <- unlist(as.list(KEGGPATHID2NAME)[names(KEGGStatLGvLV)]),
  Statistics <- round(-log2(KEGGStatLGvLV),2),
  Sig <- ifelse(KEGGStatLGvLV <=0.05,"Significant","NS")
)
colnames(pathwayres) <- c("KEGGids","Statistics","Sig")

#to add in the PDL1 and FGL1 pathways, first create dfs
#run wilcoxGST first
pPDL1<- wilcoxGST(PDL1pathway, statistics = statLGvLV, alternative='either');pPDL1
pFlg1 <- wilcoxGST(Fgl1pathway, statistics = statLGvLV, alternative='either');pFlg1
#then add in the wilcoxGST results
pathwayresadd<- data.frame(
  c("PDL1 pathway*", "FGL1-Lag3 network (GeneMania)"),
  round(-log2(c(pPDL1, pFlg1)),2),
  c("Significant","Significant")
)
colnames(pathwayresadd) <- c("KEGGids","Statistics","Sig")
pathwayresall <- rbind(pathwayres,pathwayresadd)
pathwayresall <- pathwayresall[order(pathwayresall$Statistics),]
pathwayresall$KEGGids <- factor(pathwayresall$KEGGids, levels = pathwayresall$KEGGids)
#plot
ggp <- ggplot(
  data=pathwayresall, 
  aes(x=KEGGids, y=Statistics, fill=Sig, width=0.75)) + 
  geom_bar(stat="identity") +
  xlab("") +
  ylab("-log2(q)") +
  scale_y_continuous(expand = c(0,0))+
  #geom_text(aes(label=Statistics), hjust=-0.2, vjust=0.2, size=2) +
  coord_flip() + 
  theme_classic() +
  theme(
    axis.text.y=element_text(size=6, face="bold", color="black"),
    axis.text.x=element_text(size=6, face="bold", color="black"),
    axis.title.x=element_text(size=6, face="bold"),
    legend.title = element_blank(),
    legend.spacing.x = unit(.15,'lines'),
    legend.text =element_text(size=5),
    legend.justification=c(1,0), 
    legend.position=c(0.99,0.01),
    legend.key.size = unit(0.5,'lines'))+
  scale_fill_manual(
    values=c("black","red"),
    labels=c("N.S.",expression(" q "<=" 0.05")))+
  guides(fill = guide_legend(reverse=TRUE))
pdf("Pathwaystatsresults.pdf",width=3.5,height=3); ggp; dev.off()

######################################################
#### Exporting tables for pathway visualization in Cytoscape ####

write.csv(resLGvLV[union(union(PDL1pathway,Chemokines),Immsupp),c("log2FoldChange","padj")],"valuesforpathway.csv")
write.csv(resLGvLV[Fgl1pathway,c("log2FoldChange","padj")],"valuesforFlg1pathway.csv")
write.csv(resLGvLV[unique(c(KEGGToSYMBOL[["04060"]],KEGGToSYMBOL[["04514"]],KEGGToSYMBOL[["05212"]],KEGGToSYMBOL[["04012"]],KEGGToSYMBOL[["04115"]],KEGGToSYMBOL[["04310"]]))
                   ,c("log2FoldChange","padj")],"valuesforagnosticpathway.csv")
  #write.csv(resLGvLV[KEGGToSYMBOL[["04722"]],c("log2FoldChange","padj")],"valuesforNeurpathway.csv")
write.csv(resLGvLV[KEGGToSYMBOL[["04012"]],c("log2FoldChange","padj")],"valuesforErbBpathway.csv")



######################################################
#### Boxplot generation ####

#a select set of chemokines
selectcmk<-c("Ccl2","Ccl5","Ccl20","Ccl22","Ccl28","Cxcl5","Cxcl12","Cxcl17",
"Cxcl9","Cxcl10","Cxcl11","Cxcl14")
selectcmkdf<-melt(mat[selectcmk,c(1:8,19:30)])
categ<-rep(c(rep("Pro-tumor",8),rep("Pro-immune",4)),20)
ttype<-c(rep(rep("KPC2D",length(selectcmk)),8),rep(rep("KPCLG",length(selectcmk)),5),rep(rep("KPCLV",length(selectcmk)),7))
selectcmkdf<-cbind(selectcmkdf,categ,ttype)
colnames(selectcmkdf) <- c("cmk","sample","logcount","cmktype","sgroup")

cmkplot <- ggplot(selectcmkdf, aes(x=sgroup, y=logcount, fill=sgroup))+
  geom_boxplot(outlier.size=0, lwd=0.1, color="white")+
  geom_jitter(width=0.2, size=0.2)+
  scale_fill_manual(values=c("gray20","blue","red"))+
  ylab("Log(count)")+
  ylim(-0.2,12)+
  facet_grid(~cmk)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(size=10),
        axis.ticks.y = element_line(size=0.5, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7),
        panel.border = element_rect(fill=NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(size=5)
        )
pdf("Chemokineboxplot.pdf", width=8.8, height=2.2) ; cmkplot ; dev.off()


#a select set of immune suppressors
selectimm<-c("Cd274","Pdcd1lg2","Tgfb1","Tgfb2","Tgfb3","Fasl","Il10","Fgl1","Csf2","Ido1")
selectimmdf<-melt(mat[selectimm,c(1:8,19:30)])
ttype<-c(rep(rep("KPC2D",length(selectimm)),8),rep(rep("KPCLG",length(selectimm)),5),rep(rep("KPCLV",length(selectimm)),7))
selectimmdf<-cbind(selectimmdf,ttype)
colnames(selectimmdf) <- c("imm","sample","logcount","sgroup")

immplot <- ggplot(selectimmdf, aes(x=sgroup, y=logcount, fill=sgroup))+
  geom_boxplot(outlier.size=0, lwd=0.1, color="white")+
  geom_jitter(width=0.2, size=0.2)+
  scale_fill_manual(values=c("gray20","blue","red"))+
  ylab("Log(count)")+
  ylim(-0.2,12)+
  facet_grid(~imm)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(size=10),
        axis.ticks.y = element_line(size=0.5, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7),
        panel.border = element_rect(fill=NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(size=5)
  )
pdf("Immsuppboxplot.pdf", width=8.8, height=2.2) ; immplot ; dev.off()


#a select set of ErbB signaling genes
selectErbb<-c("Egf","Tgfa","Areg","Hbegf","Btc","Hbegf","Ereg","Nrg1","Nrg2","Nrg3","Nrg4","Egfr","Erbb2","Erbb3","Erbb4")
selectErbbdf<-melt(mat[selectErbb,c(1:8,19:30)])
ttype<-c(rep(rep("KPC2D",length(selectErbb)),8),rep(rep("KPCLG",length(selectErbb)),5),rep(rep("KPCLV",length(selectErbb)),7))
selectErbbdf<-cbind(selectErbbdf,ttype)
colnames(selectErbbdf) <- c("erbb","sample","logcount","sgroup")

erbbplot <- ggplot(selectErbbdf, aes(x=sgroup, y=logcount, fill=sgroup))+
  geom_boxplot(outlier.size=0, lwd=0.1, color="white")+
  geom_jitter(width=0.2, size=0.2)+
  scale_fill_manual(values=c("gray20","blue","red"))+
  ylab("Log(count)")+
  ylim(-0.2,12)+
  facet_grid(~erbb)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(size=10),
        axis.ticks.y = element_line(size=0.5, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7),
        panel.border = element_rect(fill=NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(size=5)
  )
pdf("Erbbboxplot.pdf", width=8.8, height=2.2) ; erbbplot ; dev.off()
