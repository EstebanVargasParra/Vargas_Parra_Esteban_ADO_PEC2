if(!require(BiocManager)) install.packages("BiocManager")
if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")
if(!require(Rsamtools)) BiocManager::install("Rsamtools")
if(!require(GenomicFeatures)) BiocManager::install("GenomicFeatures")
if(!require(DESeq2)) BiocManager::install("DESeq2")
if(!require(apeglm)) BiocManager::install("apeglm")
if(!require(BiocParallel)) BiocManager::install("BiocParallel")
if(!require(genefilter)) BiocManager::install("genefilter")
if(!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
if(!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
if(!require(ReportingTools)) BiocManager::install("ReportingTools")
if(!require(RUVSeq)) BiocManager::install("RUVSeq")
if(!require(sva)) BiocManager::install("sva")
if(!require(Gviz)) BiocManager::install("Gviz")

if(!require(magrittr)) install.packages("magrittr", dep=TRUE)
if(!require(dplyr)) install.packages("dplyr", dep=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dep=TRUE)
if(!require(pheatmap)) install.packages("pheatmap", dep=TRUE)
if(!require(RColorBrewer)) install.packages("RColorBrewer", dep=TRUE)
if(!require(ggbeeswarm)) install.packages("ggbeeswarm", dep=TRUE)

setwd(".")
targets<- read.csv("./targets.csv", header = TRUE, sep = ",")
head(targets)

counts <- read.csv("./counts.csv", header = FALSE, sep = ";") 
counts[1,1]<-NA
counts$V1<-gsub("\\..*", "", counts$V1, fixed = FALSE)
counts$V1<-as.character(counts$V1)
counts$V1[is.na(counts$V1)]<-"Sample_Name"
head(counts)
female<-subset(targets, targets$sex=="female")
male<-subset(targets, targets$sex=="male")
#-----------------------------------------------------------------------------
library(dplyr)
set.seed(19950605)
NITf <- sample_n(female[female$Group=="NIT",], size = 5)
NITm <- sample_n(male[male$Group=="NIT",], size = 5)
SFIf <- sample_n(female[female$Group=="SFI",], size = 5)
SFIm <- sample_n(male[male$Group=="SFI",], size = 5)
ELIf <- sample_n(female[female$Group=="ELI",], size = 5)
ELIm <- sample_n(male[male$Group=="ELI",], size = 5)
target<- data.frame(rbind(NITf, NITm, SFIf, SFIm, ELIf, ELIm))
#-----------------------------------------------------------------------------
dim(counts)
rownames(counts)<-counts[,1]
transpuesta<-t(counts)
transpuesta<-data.frame(transpuesta)
dim(transpuesta)
transpuesta<-transpuesta[-1,]
rownames(transpuesta)<- 1:nrow(transpuesta)

library(plyr)
Sample_Name<-target$Sample_Name
ShortName <-target$ShortName 
Sample_Name<-data.frame(Sample_Name, ShortName)
objetivos<- join(Sample_Name, transpuesta)
objetivos<-objetivos[,c(-1)]
rownames(objetivos)<-objetivos[,1]
objetivos<-objetivos[,-1]
objetivos[sapply(objetivos, is.factor)]<- lapply(objetivos[sapply(objetivos, is.factor)],function(x) as.numeric(as.character(x)))
tobjetivos<- t(objetivos)
rownames(target)<-target[,length(target)]
target<-target[,-length(target)]

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = tobjetivos,
                                 colData = target,
                                 design = ~Group)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
colData(rld)

library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

sampleDists <- dist(t(assay(vsd)))
sampleDists


library("pheatmap")
library("RColorBrewer")
library("ggbeeswarm")

sampleDistMatrix <- as.matrix( sampleDists)
rownames(sampleDistMatrix) <- paste(rownames(target))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

plotPCA(vsd, intgroup ="Group")


mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Group, shape = Group)) +
  geom_point(size = 3) + coord_fixed()

dds <- DESeq(dds, parallel =TRUE)

NIT_SFI <- results(dds, contrast = c("Group", "NIT", "SFI"))
mcols(NIT_SFI, use.names = TRUE)
summary(NIT_SFI)
sum(NIT_SFI$padj < 0.1, na.rm=TRUE)
resSig_NIT_SFI <- subset(NIT_SFI, padj < 0.1)
head(resSig_NIT_SFI[ order(resSig_NIT_SFI$log2FoldChange, decreasing = TRUE), ])
topGene_NIT_SFI <- rownames(NIT_SFI)[which.min(NIT_SFI$padj)]
geneCounts_NIT_SFI <- plotCounts(dds, gene = topGene_NIT_SFI, intgroup =c("Group"),
                         returnData = TRUE)
geneCounts_NIT_SFI <- geneCounts_NIT_SFI[1:20,]
ggplot(geneCounts_NIT_SFI, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + ggtitle("NIT vs SFI")+
  theme(plot.title = element_text(hjust = 0.5))

NIT_ELI <- results(dds, contrast = c("Group", "NIT", "ELI"))
mcols(NIT_ELI, use.names = TRUE)
summary(NIT_ELI)
sum(NIT_ELI$padj < 0.1, na.rm=TRUE)
resSig_NIT_ELI <- subset(NIT_ELI, padj < 0.1)
head(resSig_NIT_ELI[ order(resSig_NIT_ELI$log2FoldChange), ])
head(resSig_NIT_ELI[ order(resSig_NIT_ELI$log2FoldChange, decreasing = TRUE), ])
topGene_NIT_ELI <- rownames(NIT_ELI)[which.min(NIT_ELI$padj)]
geneCounts_NIT_ELI <- plotCounts(dds, gene = topGene_NIT_ELI, intgroup = c("Group"),
                                 returnData = TRUE)
geneCounts_NIT_ELI <- geneCounts_NIT_ELI[-c(11:20),]
ggplot(geneCounts_NIT_ELI, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + ggtitle("NIT vs ELI")+
  theme(plot.title = element_text(hjust = 0.5))

SFI_ELI <- results(dds, contrast = c("Group", "SFI", "ELI"))
mcols(SFI_ELI, use.names = TRUE)
summary(SFI_ELI)
sum(SFI_ELI$padj < 0.1, na.rm=TRUE)
resSig_SFI_ELI <- subset(SFI_ELI, padj < 0.1)
head(resSig_SFI_ELI[ order(resSig_SFI_ELI$log2FoldChange, decreasing = TRUE), ])
topGene_SFI_ELI <- rownames(SFI_ELI)[which.min(SFI_ELI$padj)]
geneCounts_SFI_ELI <- plotCounts(dds, gene = topGene_SFI_ELI, intgroup = c("Group"),
                                 returnData = TRUE)
geneCounts_SFI_ELI <- geneCounts_SFI_ELI[-c(1:10),]
ggplot(geneCounts_SFI_ELI, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + ggtitle("SFI vs ELI")+
  theme(plot.title = element_text(hjust = 0.5))

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)


resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

sum(res$padj < 0.1, na.rm=TRUE)

library("apeglm")

res_NIT_SFI <- lfcShrink(dds, contrast =c("Group", "NIT", "SFI"), res =NIT_SFI)
plotMA(res_NIT_SFI, ylim = c(-5,5))
topGene_NIT_SFI <- rownames(res_NIT_SFI)[which.min(res_NIT_SFI$padj)]
with(res[topGene_NIT_SFI, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene_NIT_SFI, pos=2, col="dodgerblue")
})

res_NIT_ELI <- lfcShrink(dds, contrast =c("Group", "NIT", "ELI"), res =NIT_ELI)
plotMA(res_NIT_ELI, ylim = c(-5, 5))
topGene_NIT_ELI <- rownames(res_NIT_ELI)[which.min(res_NIT_ELI$padj)]
with(res[topGene_NIT_ELI, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene_NIT_ELI, pos=2, col="dodgerblue")
})

res_SFI_ELI <- lfcShrink(dds, contrast =c("Group", "SFI", "ELI"), res =SFI_ELI)
plotMA(res_SFI_ELI, ylim = c(-5, 5))
topGene_SFI_ELI <- rownames(res_SFI_ELI)[which.min(res_SFI_ELI$padj)]
with(res[topGene_SFI_ELI, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene_SFI_ELI, pos=2, col="dodgerblue")
})

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[topVarGenes,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group","body_site")])
pheatmap(mat, annotation_col =anno)


library("org.Hs.eg.db")
library("AnnotationDbi")
columns(org.Hs.eg.db)

res_NIT_SFI$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_NIT_SFI),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res_NIT_SFI$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_NIT_SFI),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered_NIT_SFI <- res_NIT_SFI[order(res_NIT_SFI$pvalue),]
head(resOrdered_NIT_SFI)
resOrderedDF_NIT_SFI <- as.data.frame(resOrdered_NIT_SFI)

res_NIT_ELI$symbol <- mapIds(org.Hs.eg.db,
                             keys=row.names(res_NIT_ELI),
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
res_NIT_ELI$entrez <- mapIds(org.Hs.eg.db,
                             keys=row.names(res_NIT_ELI),
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
resOrdered_NIT_ELI <- res_NIT_ELI[order(res_NIT_ELI$pvalue),]
head(resOrdered_NIT_ELI)
resOrderedDF_NIT_ELI <- as.data.frame(resOrdered_NIT_ELI)


res_SFI_ELI$symbol <- mapIds(org.Hs.eg.db,
                             keys=row.names(res_SFI_ELI),
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
res_SFI_ELI$entrez <- mapIds(org.Hs.eg.db,
                             keys=row.names(res_SFI_ELI),
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
resOrdered_SFI_ELI <- res_SFI_ELI[order(res_SFI_ELI$pvalue),]
head(resOrdered_SFI_ELI)
resOrderedDF_SFI_ELI <- as.data.frame(resOrdered_SFI_ELI)



library("ReportingTools")
htmlRep_NIT_SFI <- HTMLReport(shortName="report_NIT_SFI", title="My report NIT vs SFI",
                      reportDirectory=".")
publish(resOrderedDF_NIT_SFI, htmlRep_NIT_SFI)
url_NIT_SFI <- finish(htmlRep_NIT_SFI)
browseURL(url_NIT_SFI)

htmlRep_NIT_ELI <- HTMLReport(shortName="report_NIT_ELI", title="My report NIT vs ELI",
                              reportDirectory=".")
publish(resOrderedDF_NIT_ELI, htmlRep_NIT_ELI)
url_NIT_ELI <- finish(htmlRep_NIT_ELI)
browseURL(url_NIT_ELI)

htmlRep_SFI_ELI <- HTMLReport(shortName="report_SFI_ELI", title="My report SFI vs ELI",
                              reportDirectory=".")
publish(resOrderedDF_SFI_ELI, htmlRep_SFI_ELI)
url_SFI_ELI <- finish(htmlRep_SFI_ELI)
browseURL(url_SFI_ELI)

library(clusterProfiler)
library(DOSE)
OrgDB<- org.Hs.eg.db
NIT_SFI_plot<- subset(resOrderedDF_NIT_SFI, padj < 0.1)
enrich_NIT_SFI <- groupGO(gene =NIT_SFI_plot$entrez,
                          OrgDb = OrgDB,
                          keyType = "ENTREZID",
                          ont = "BP",
                          readable = TRUE)
head(enrich_NIT_SFI)
cnetplot(enrich_NIT_SFI, categorySize = "geneNum")
barplot(enrich_NIT_SFI, showCategory = 15, font.size = 7, title = "Pathway Analysis for ",)


NIT_ELI_plot<-subset(resOrderedDF_NIT_ELI, padj < 0.1)
enrich_NIT_ELI <- groupGO(gene =NIT_ELI_plot$entrez,
                          OrgDb = OrgDB,
                          keyType = "ENTREZID",
                          ont = "BP",
                          readable = TRUE)
head(enrich_NIT_ELI)
#cnetplot(enrich_NIT_ELI, categorySize = "geneNum")
barplot(enrich_NIT_ELI, showCategory = 15, font.size = 7, title = "Pathway Analysis for ",)


SFI_ELI_plot<-subset(resOrderedDF_SFI_ELI, padj < 0.1)
enrich_SFI_ELI <- groupGO(gene =SFI_ELI_plot$entrez,
                          OrgDb = OrgDB,
                          keyType = "ENTREZID",
                          ont = "BP",
                          readable = TRUE)
head(enrich_SFI_ELI)
#cnetplot(enrich_SFI_ELI, categorySize = "geneNum")
barplot(enrich_SFI_ELI, showCategory = 15, font.size = 7, title = "Pathway Analysis for ",)
