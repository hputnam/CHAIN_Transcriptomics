#Title: P. acuta nutrient exposure
#Project: CHAIN
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20210125
#See Readme file for details

#call the DESeq2 library
#BiocManager::install("BiocUpgrade")
#BiocManager::install("DESeq2")
library("DESeq2")
#install.packages("fdrtool")
library(fdrtool)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("reshape2")
library(reshape2)
library(ggplot2)
#install.packages("tm")
#library(tm)
library(genefilter)
#install.packages("pheatmap")
library(pheatmap)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("limma")
library(limma)
library(spdep)
library(adegenet)
#BiocManager::install("GSEABase")
library(GSEABase)
#BiocManager::install("pathview")
#library("pathview")
#BiocManager::install("rtracklayer")
library("rtracklayer")
#BiocManager::install("goseq")
library("goseq")
library("GO.db")
#devtools::install_github("hms-dbmi/UpSetR")
library("UpSetR")
library(lattice)
library(latticeExtra)
library(IRanges)
library(GenomicRanges)
library(sqldf)
library(rentrez)


#Load annotation file
Annot <- read.csv("Data/final.trinotate_annotation_report.xls", header=TRUE, sep="\t", na.string="NA", stringsAsFactors = F)
colnames(Annot)[1] <- "Trinity.gene" 
Annot <- Annot[!duplicated(Annot$Trinity.gene), ]

#####
#Read in data
Ahya_TranscriptCountData <- as.data.frame(read.csv(file= "Data/Trinity/RSEM.gene.counts.matrix", sep="\t", header = TRUE))
head(Ahya_TranscriptCountData)
rownames(Ahya_TranscriptCountData) <- Ahya_TranscriptCountData[,1]
Ahya_TranscriptCountData <- Ahya_TranscriptCountData[,-1]
Ahya_TranscriptCountData <- round(Ahya_TranscriptCountData,0)

###filtering values for PoverA
#set filter values for PoverA, P percent of the samples have counts over A
filt <- filterfun(pOverA(0.90,10))

#create filter for the counts data
tfil <- genefilter(Ahya_TranscriptCountData, filt)

#identify transcripts to keep by count filter
keep <- Ahya_TranscriptCountData[tfil,]

#identify transcript list
gn.keep <- rownames(keep)

#data filtered in PoverA, P percent of the samples have counts over A
counts.5x <- as.matrix(Ahya_TranscriptCountData[which(rownames(Ahya_TranscriptCountData) %in% gn.keep),])
#write.csv(counts.5x, file="Output/filtered_counts.csv")

Ahya_sample_ColData <- read.csv("Data/Ahya_Seq_Sample_Info.csv", header=TRUE, sep=",")
print(Ahya_sample_ColData)

#change rownames to match
rownames(Ahya_sample_ColData) <- Ahya_sample_ColData$Tube_ID
colnames(counts.5x) <- colnames(counts.5x)
head(Ahya_sample_ColData)
head(counts.5x)

# Check all sample IDs in Ahya_sample_ColData are also in Ahya_TranscriptCountData and match their orders
all(rownames(Ahya_sample_ColData) %in% colnames(counts.5x))  #Should return TRUE
# returns TRUE
all(rownames(Ahya_sample_ColData) == colnames(counts.5x))    # should return TRUE
#returns TRUE

#give the treatment column levels
Ahya_sample_ColData$Site <- factor(Ahya_sample_ColData$Site)
levels(Ahya_sample_ColData$Site)


#### Construct DESeq dataset from matrix
data <- DESeqDataSetFromMatrix(countData = counts.5x, colData = Ahya_sample_ColData, design = ~ Site) #create a DESeqDataSet object


#apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
rld <- rlog(data, blind=FALSE)
head(assay(rld), 3) #view data

#calculate distance matix
sampleDists <- dist(t(assay(rld)))

#distance matrix
sampleDistMatrix <- as.matrix(sampleDists)

#assign row names
rownames(sampleDistMatrix) <- colnames(rld)

#assign col names
colnames(sampleDistMatrix) <- NULL

#assign colors
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


#### Plotting heatmap for Expression Visualization

pdf(file="Output/heatmap_all.pdf")
pheatmap(sampleDistMatrix, #plot matrix of expression similarity
         clustering_distance_rows=sampleDists, #cluster rows
         clustering_distance_cols=sampleDists, #cluster columns
         col=colors) #set colors
dev.off()

#### Plotting PCA from Expression visualization

pdf(file="Output/PCA_all.pdf")
plotPCA(rld, intgroup = c("Site")) #plot PCA of samples with all data
dev.off()

#### Differential Gene Expression Analysis
# Interaction Test: test of the factor of "Site" 

#run differential expression test by group using the wald test
DEG.int <- DESeq(data)

#save DE results
DEG.int.res <- results(DEG.int)

# Volcano Plot (significance as a function of fold change)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- as.data.frame(DEG.int.res)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value (pajd))))

#Color the significant points with log2 fold change >2 red ()
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#Add lines for absolute FC>2
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
#abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

#view DE results
resultsNames(DEG.int)

log2_changes <- DEG.int.res[,c(2,6)]


#identify the number of significant p values with 5%FDR (padj<0.05)
sig.num <- sum(DEG.int.res$padj <0.05, na.rm=T)

#identify signficant pvalues with 5%FDR
sig <- subset(DEG.int.res, padj<0.05,)

#subset list of sig transcripts from original count data
sig.list <- data[which(rownames(data) %in% rownames(sig)),]

#apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
rsig <- rlog(sig.list, blind=FALSE)
write.csv(counts(sig.list), file="Output/DEG_5FDR.all.csv")

#### Plotting DEG PCA

PCA.plot <- plotPCA(rsig, intgroup = c("Site")) #Plot PCA of all samples for DEG only
PCA.plot #view plot
PC.info <-PCA.plot$data #extract plotting data
pdf(file="Output/PCA.DEG.pdf")
plot(PC.info$PC1, PC.info$PC2, xlim=c(-16,20), ylim=c(-9,9), xlab="PC1 52%", ylab="PC2 12%",  pch =16, col = c("gray", "green")[as.numeric(PC.info$Site)], cex=2.5)
legend(x="topright",
       bty="n",
       legend = c("Forereef", "BackReef"),
       text.col = c("black","black"),
       pch = c(16, 16),
       col = c("gray","green"),
       cex=2)
dev.off()

#### Plotting DEG heatmap

#sort by decreasing sig
topVarGenes <- head(order(rowVars(assay(rsig)),decreasing=TRUE),sig.num)

#make an expression object
mat <- assay(rsig)[ topVarGenes, ]

#difference in expression compared to average across all samples
mat <- mat - rowMeans(mat)

#ordering columns by treatment
col.order <- c("A179","A19","A191","A227","A760", 
               "A167","A203","A215","A784","A796")
mat <- mat[,col.order]

#make dataframe
df <- data.frame(colData(rsig)[c("Site")])

dev.off()

pdf(file="Output/DEG_Heatmap.pdf")
pheatmap(mat, annotation_col = df, clustering_method = "average",
                clustering_distance_rows="euclidean", show_rownames =TRUE, cluster_cols=TRUE,
                show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()

DEGs <- as.data.frame(counts(sig.list))
DEGs$Trinity.gene <- rownames(DEGs)
Output <- merge(DEGs,Annot,  by="Trinity.gene", all = FALSE)

Web <- read.csv("Data/DEG_Alignment-HitTable_blastn.csv", header=FALSE, sep=",", na.string="NA", stringsAsFactors = T)
Web <- Web[!duplicated(Web$V1),]
colnames(Web)[1] <- "Trinity.gene"
Web$Trinity.gene <- gsub("_[^_]+$","",Web$Trinity.gene) #remove extra characters
#Output <- merge(Output, Web,  by="Trinity.gene", all = TRUE)

#for loop to extract IDs and titles from accession numbers and write out accession and title
results <- data.frame()

for(i in 1:length(Web$V2)){
r_search <- entrez_search(db="nucleotide", term=Web$V2[i])
ids <- r_search$id
esummary <- entrez_summary(db="nucleotide", id=ids)

# capture summary stats to data frame
df <- data.frame(Trinity.gene =Web$Trinity.gene[i],
                 description =Web$V2[i],
                 accession = esummary$title)


# bind rows of temporary data frame to the results data frame
results <- rbind(results, df)
}

DEG.Output <- merge(Output, results,  by="Trinity.gene", all = TRUE)
write.csv(Output, file="Output/DEG_Output_Annot.csv")



sig.log2.changes <- subset(log2_changes, padj<0.05,)
sig.log2.changes <- as.data.frame(sig.log2.changes)
sig.log2.changes$Trinity.gene <- rownames(sig.log2.changes)
sig.log2.changes$log2FoldChange <- as.numeric(as.character(sig.log2.changes$log2FoldChange))
sig.log2.changes <- sig.log2.changes[order(sig.log2.changes$log2FoldChange),]

DEG.FC <- merge(DEG.Output, sig.log2.changes,  by="Trinity.gene", all = FALSE)
DEG.FC$accession <- gsub("PREDICTED: Acropora millepora ", "", DEG.FC$accession)
DEG.FC$accession <- gsub("PREDICTED: Acropora digitifera ", "", DEG.FC$accession)
DEG.FC$accession<- gsub("\\(.*","",DEG.FC$accession)
DEG.FC <- DEG.FC[!grepl("uncharacterized", DEG.FC$accession), ]
DEG.FC<- droplevels(DEG.FC[!grepl("uncharacterized", DEG.FC$accession), ])
DEG.FC <- DEG.FC[order(DEG.FC$log2FoldChange),]

cols <- c("black", "orange", "blue","grey")

annotation <- data.frame(
  x = c(10,35),
  y = c(8,8),
  label = c("Higher in the Forereef", "Higher in the Backreef"))

pdf(file="Output/Egg_Expression_FC.pdf")
ggplot(DEG.FC, aes(reorder(x=Trinity.gene,-(-log2FoldChange)), y=log2FoldChange), 
       fill=factor(ifelse(DEG.FC$accession=="histone H2A.V",T,F))) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 25.5, linetype="dotted") +
  scale_fill_manual(name = "accession",values = cols) +
  geom_label(data=annotation, aes(x=x, y=y, label=c("Higher in the Forereef", "Higher in the Backreef")),                 , 
             color="black", size=3 , angle=0, fontface="bold") +
  scale_x_discrete(labels = DEG.FC$accession) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.6))
dev.off()

geom_label(data=annotation, aes( x=x, y=y, label=label),                 , 
           color="orange", 
           size=7 , angle=45, fontface="bold" )

# ##### GO enrichment of maternally provisioned genes using MWU #####
# #GO MWU approach https://github.com/z0on/GO_MWU which shows GO categories significantly enriched by either up or down-regulated genes#####
# 
# GO.terms <-read.table("Data/go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
# colnames(GO.terms) <- c("Trinity.gene", "GO.terms") #label columns for merging
# GO.terms$Trinity.gene <-gsub("_i.*","",GO.terms$Trinity.gene) #remove transcript information to retain gene only
# Ahya_egg_GO.terms <- GO.terms[!duplicated(GO.terms$Trinity.gene), ]
# 
# #identify MWU data
# ALL.MWU <- Ahya_egg_GO.terms
# ALL.MWU <- ALL.MWU[,c(1,2)]
# ALL.MWU$GO.terms <-gsub(",","; ",ALL.MWU$GO.terms)
# #write MWU data to file
# write.table(ALL.MWU, 'GO_MWU/Ahya_egg.GO.MWU', sep='\t', row.names=FALSE)
# #merge all genes with DMGs
# 
# log2_changes <- cbind(rownames(DEG.int.res), DEG.int.res$log2FoldChange)
# colnames(log2_changes) <- c("Trinity.gene", "log2FoldChange") 
# write.table(log2_changes, 'GO_MWU/Ahya_egg.foldchange.MWU', sep=',', row.names=FALSE)


##### GO enrichment of DE genes #####
ALL<-row.names(data) #set the all transcripts list
DEG <- as.character(row.names(sig.list)) #set the enrichment test list
IDs <- row.names(data) #set ID names
LENGTH <-read.table("Data/Trinity.gene_lengths.txt", sep = "", header=F) #reads in a list of gene lengths
gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
colnames(gn.keep) <- c("V1") #name columns
lens <- as.data.frame(LENGTH) #set lengths
LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
colnames(LENGTH) <- c("Trinity.gene", "length")
colnames(gn.keep) <- "Trinity.gene" #label column for merging

#load GO terms from Trinotate
GO.terms <-read.table("Data/go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
colnames(GO.terms) <- c("Trinity.gene", "GO.terms") #label columns for merging
GO.terms$Trinity.gene <-gsub("_i.*","",GO.terms$Trinity.gene) #remove transcript information to retain gene only
GO.terms <- merge(gn.keep, GO.terms, by="Trinity.gene", all=TRUE) #subset list to genes that are used for full analysis
GO.terms$GO.terms[is.na(GO.terms$GO.terms)] <- "unknown" #list NA as unknown
splitted <- strsplit(as.character(GO.terms$GO.terms), ",") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.terms$Trinity.gene, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript

ALL <- GO.terms

#change contig lists to vectors
ALL.vector <-c(t(ALL))
DEG.vector <-c(t(DEG))
ID.vector <- LENGTH$Trinity.gene
LENGTH.vector <-LENGTH$length

ID.vector <- LENGTH$Trinity.gene
gene.vector <-  as.integer(ID.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
names(gene.vector)=ID.vector #set names
DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
tail(GO.wall)
nrow(GO.wall)

#Find only enriched GO terms that are statistically significant at cutoff 
enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
enriched.GO.05<-data.frame(enriched.GO.05.a)
colnames(enriched.GO.05) <- c("category")
enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")


MF.INT <- subset(enriched.GO.05, ontology=="MF")
MF.INT <- MF.INT[order(MF.INT$numDEInCat),]
CC.INT <- subset(enriched.GO.05, ontology=="CC")
CC.INT <- CC.INT[order(CC.INT$numDEInCat),]
BP.INT <- subset(enriched.GO.05, ontology=="BP")
BP.INT <- BP.INT[order(BP.INT$numDEInCat),]
write.csv(ALL , file = "Output/Transcriptome_GOTerms.csv")
write.csv(MF.INT , file = "Output/MF_Sig_Enriched_GO.05_INT.csv")
write.csv(CC.INT , file = "Output/CC_Sig_Enriched_GO.05_INT.csv")
write.csv(BP.INT , file = "Output/BP_Sig_Enriched_GO.05_INT.csv")

pdf(file="Output/GO_enrichment.pdf")
par(mfrow=c(3,1))
par(mar=c(3,1,1,0), oma = c(0.1, 15, 0.1, 0.5))
barplot(BP.INT$numDEInCat, horiz = T, col="black", names.arg=BP.INT$term, las=1, cex.names = 0.4, xlim=c(0,10), main = "Biological Process")
text(labels = 'Biological Process', xpd = NA, srt = 90, x=-26, y=15, cex=1)
barplot(MF.INT$numDEInCat, horiz = T, col="black", names.arg=MF.INT$term, las=1, cex.names = 0.4, xlim=c(0,10), main = "Molecular Function")
text(labels = 'Molecular Function', xpd = NA, srt = 90, x=-26, y=15, cex=1)
barplot(CC.INT$numDEInCat, horiz = T, col="black", names.arg=CC.INT$term, las=1, cex.names = 0.4, xlim=c(0,10), main = "Cellular Component")
text(labels = 'Cellular Component', xpd = NA, srt = 90, x=-26, y=15, cex=1)
dev.off()

pdf(file="Output/GO_enrichment_BPMF.pdf")
par(mfrow=c(2,1))
par(mar=c(3,1,1,0), oma = c(0.1, 15, 0.1, 0.5))
barplot(BP.INT$numDEInCat, horiz = T, col="black", names.arg=BP.INT$term, las=1, cex.names = 0.4, xlim=c(0,6), main = "Biological Process")
text(labels = 'Biological Process', xpd = NA, srt = 90, x=-26, y=20, cex=0.4)
barplot(MF.INT$numDEInCat, horiz = T, col="black", names.arg=MF.INT$term, las=1, cex.names = 0.4, xlim=c(0,6), main = "Molecular Function")
text(labels = 'Molecular Function', xpd = NA, srt = 90, x=-26, y=20, cex=1)
dev.off()

BP.INT$numDEInCat
BP.INT$prop <- (BP.INT$numDEInCat/sum(BP.INT$numDEInCat))*100
pdf(file="Output/Transcriptome_Bar_GO_BP.pdf")
ggplot() + 
  geom_bar(aes(y = prop, x=1,  fill = term), data = BP.INT, stat="identity") +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
dev.off()

##### Slim for ALL #####
All.Ids <- as.character(GO.wall$category)
myCollection <- GOCollection(All.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
All.BP.slim <- goSlim(myCollection, slim, "BP")
All.BP.slim  <-All.BP.slim[,c(1,3)]
All.BP.slim  <- All.BP.slim[order(-All.BP.slim$Count),]
All.BP.slim <- All.BP.slim[-1,]
All.BP.slim  <- All.BP.slim[apply(All.BP.slim!=0, 1, all),]
All.BP.slim$prop <- (All.BP.slim$Count/sum(All.BP.slim$Count))*100

All.MF.slim <- goSlim(myCollection, slim, "MF")
All.MF.slim  <-All.MF.slim[,c(1,3)]
All.MF.slim  <- All.MF.slim[order(-All.MF.slim$Count),]
All.MF.slim <- All.MF.slim[-1,]
All.MF.slim  <- All.MF.slim[apply(All.MF.slim!=0, 1, all),]
All.MF.slim$prop <- (All.MF.slim$Count/sum(All.MF.slim$Count))*100

All.CC.slim <- goSlim(myCollection, slim, "CC")
All.CC.slim  <-All.CC.slim[,c(1,3)]
All.CC.slim  <- All.CC.slim[order(-All.CC.slim$Count),]
All.CC.slim <- All.CC.slim[-1,]
All.CC.slim  <- All.CC.slim[apply(All.CC.slim!=0, 1, all),]
All.CC.slim$prop <- (All.CC.slim$Count/sum(All.CC.slim$Count))*100

pdf(file="Output/Transcriptome_GO.pdf")
par(mfrow=c(2,1))
par(mar=c(4,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
barplot(All.BP.slim$prop, horiz = T, col="black", names.arg=All.BP.slim$Term, las=1, cex.names = 0.4, xlim=c(0,20))
text(labels = 'Biological Process', xpd = NA, srt = 90, x=-5, y=40, cex=2)
barplot(All.MF.slim$prop, horiz = T, col="black", names.arg=All.MF.slim$Term, las=1, cex.names = 0.4, xlim=c(0,20))
text(labels = 'Molecular Function', xpd = NA, srt = 90, x=-5, y=25, cex=2)
#barplot(All.CC.slim$prop, horiz = T, col="black", names.arg=All.CC.slim$Term, las=1, cex.names = 0.4, xlim=c(0,20), xlab="Proportion of Genes in Category")
#text(labels = 'Cellular Component', xpd = NA, srt = 90, x=-5, y=25, cex=2)
dev.off()

pdf(file="Output/Transcriptome_Bar_GO.pdf")
ggplot() + 
  geom_bar(aes(y = prop, x=1,  fill = Term), data = All.BP.slim, stat="identity") +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
dev.off()

pdf(file="Output/Transcriptome_Pie_GO.pdf")
ggplot(All.BP.slim, aes(x = 1, y = prop, fill = Term)) +
  geom_col() +
  coord_polar(theta = "y") +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
dev.off()

##### Slim for Biological Process #####
BP.Ids <- as.character(BP.INT$category)
myCollection <- GOCollection(BP.Ids)
slim <- getOBOCollection(fl, evidenceCode="TAS")
BP.slim.INT <- goSlim(myCollection, slim, "BP")
BP.slim.INT <-BP.slim.INT[,c(1,3)]
BP.slim.INT <- BP.slim.INT[order(-BP.slim.INT$Count),]
BP.slim.INT <- BP.slim.INT[apply(BP.slim.INT!=0, 1, all),]
BP.slim.INT <- BP.slim.INT[-1,]
  
##### Slim for Molecular Function #####
MF.Ids <- as.character(MF.INT$category)
myCollection <- GOCollection(MF.Ids)
slim <- getOBOCollection(fl, evidenceCode="TAS")
MF.slim.INT <- goSlim(myCollection, slim, "MF")
MF.slim.INT <-MF.slim.INT[,c(1,3)]
MF.slim.INT <- MF.slim.INT[order(-MF.slim.INT$Count),]
MF.slim.INT <- MF.slim.INT[apply(MF.slim.INT!=0, 1, all),]
MF.slim.INT <- MF.slim.INT[-1,]

##### Slim for Cellular Component #####
CC.Ids <- as.character(CC.INT$category)
myCollection <- GOCollection(CC.Ids)
slim <- getOBOCollection(fl, evidenceCode="TAS")
CC.slim.INT <- goSlim(myCollection, slim, "CC")
CC.slim.INT <-CC.slim.INT[,c(1,3)]
CC.slim.INT <- CC.slim.INT[order(-CC.slim.INT$Count),]
CC.slim.INT <- CC.slim.INT[apply(CC.slim.INT!=0, 1, all),]
CC.slim.INT <- CC.slim.INT[-1,]

pdf(file="Output/GOSlim_enrichment.pdf")
#par(mfrow=c(3,1))
par(mar=c(3,1,1,0), oma = c(0.1, 15, 0.1, 0.5))
barplot(BP.slim.INT$Count, horiz = T, col="black", names.arg=BP.slim.INT$Term, las=1, cex.names = 0.4, xlim=c(0,15), main = "Biological Process")
text(labels = 'Biological Process', xpd = NA, srt = 90, x=-26, y=15, cex=1)
# barplot(MF.slim.INT$Count, horiz = T, col="black", names.arg=MF.slim.INT$Term, las=1, cex.names = 0.4, xlim=c(0,4), main = "Molecular Function")
# text(labels = 'Molecular Function', xpd = NA, srt = 90, x=-26, y=15, cex=1)
# barplot(CC.slim.INT$Count, horiz = T, col="black", names.arg=CC.slim.INT$Term, las=1, cex.names = 0.4, xlim=c(0,21), main = "Cellular Component")
# text(labels = 'Cellular Component', xpd = NA, srt = 90, x=-26, y=15, cex=1)
dev.off()
