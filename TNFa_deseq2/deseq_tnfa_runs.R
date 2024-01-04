library(DESeq2)
library(tidyverse)
library(data.table)

tiss<-'lung'
wd <-paste0('/scratch/Users/tajo5912/perturbation_master/deseq2/')
setwd(wd)
outdir<-paste0(wd,tiss,'/')
dir.create(outdir)

metadatafile=paste0("/scratch/Users/tajo5912/perturbation_master/deseq2/",tiss,"_metadata.txt")

indir="/scratch/Users/tajo5912/perturbation_master/counts/"

coveragetablefile=paste0(indir,tiss,'_counts.txt')

cutoff=0.05

#read the metadata
metadata <- read.table(metadatafile, header=TRUE, sep="\t", fill=TRUE)
metadata

ids <- dplyr::pull(metadata, srr)

metadata$treatment <- factor(metadata$treatment)
metadata$treatment <- relevel(metadata$treatment, 'control')
metadata$treatment

########################################################################################################################
########################################################################################################################
#read counts
coveragetable <- read.table(coveragetablefile, sep='\t', header=TRUE)
rownames(coveragetable) <- coveragetable$ROI
countdata <- coveragetable %>% dplyr::select(as.vector(metadata$srr))
head(coveragetable)
colnames(coveragetable)

head(countdata)
dim(countdata)

########################################################################################################################
#########################################################################################################################
# run deseq
dds <- DESeqDataSetFromMatrix(countData =  countdata, 
                              colData = metadata, 
                              design = ~treatment) #use this for all samples
########################################################################################################################
########################################################################################################################
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds) #adds mu to dds <-estimateDispersionsFit(ddsr)
dds <- estimateDispersionsFit(dds)
dds <- estimateDispersionsMAP(dds)
dds <- nbinomWaldTest(dds)
########################################################################################################################
########################################################################################################################
normcounts <- as.data.frame(counts(dds, normalize=TRUE))
head(normcounts)

write.csv(as.data.frame(normcounts), file=paste(outdir,"normalizedcounts.csv",sep=""))

sink(paste0(outdir, "size_factors.txt"))
sizeFactors(dds)
sink()

png(paste0(outdir, "qc-dispersions.png", step=""), 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

rld <- rlogTransformation(dds)
head(assay(rld))

sampleDists <- dist(t(assay(rld) ) )
png(paste0(outdir,"qc-heatmap-samples.png"), w=1000, h=1000, pointsize=20)
heatmap(as.matrix(sampleDists), key=F, trace="none", main="Sample Distance Matrix")
dev.off()

png(paste0(outdir, "qc-pca_tx.png"), 1000, 1000, pointsize=20)
DESeq2::plotPCA(rld, intgroup="treatment")
dev.off()


################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
dds$treatment
resultsNames(dds)
#"Intercept"                 "treatment_TNFa_vs_control"

#(control then treatment, treatment/control)
condx <- "TNFa" #this one is the top of the ratio (treatment)
condy <- "control" #this one is the bottom of the ratio (control)

# Get differential expression results, specify which conditions to contrast (max=2)
res <- results(dds, contrast = c("treatment", condx, condy))

# merging res & the counts used (dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
# makinf row names a column called "ROI"
## Order by adjusted p-value
res <- res[order(res$padj), ]

names(resdata)[1] <- "ROI"

head(res, n=5)
#plotCounts(dds,'bidir_23167;DMSO_20_R2', intgroup=c('genotype','treatment'))

resdata <- mutate(resdata,ROI=as.character(ROI))
resdata <- resdata %>%
  select(ROI, everything())
head(resdata)
dim(resdata)

### Write results
write.table(resdata, paste0(outdir, file="diffexpr-results_DESeq2.txt"), append = FALSE, sep = "\t", row.names=FALSE)

resSig <- resdata[ resdata$padj < cutoff, ]
dim(resSig)
#head(resSig)
up <- subset(resSig, log2FoldChange > 0)
up <- na.omit(up)
up <- mutate(up,ROI=as.character(ROI))
up <- up %>%
  select(ROI, everything())

dim(up)
write.table(up[order(up$log2FoldChange, decreasing=TRUE),],
           paste0(outdir, "upregulated.txt"),row.names=FALSE)

down <- subset(resSig, log2FoldChange < 0)
down <- na.omit(down)
down <- mutate(down,ROI=as.character(ROI))
down <- down %>%
  select(ROI, everything())

dim(down)
write.table(down[order(down$log2FoldChange, decreasing=TRUE),],
           paste0(outdir, "downregulated.txt"),row.names=FALSE)
# ##############################################################################################################################################
# ##############################################################################################################################################
# ## Plotting
## Examine plot of p-values
png(paste0(outdir, 'hist_pval.png'))
hist(res$padj, breaks=50, col="grey")
dev.off()
# 
# ## Check independent filtering
# metadata(res)$filterThreshold
# 
# ## Get our data ready to plot... set thresholds
de_data <- resdata
rownames(de_data) <- de_data[,1]
de_data <- de_data[-c(1)]
# 
signifROI <- subset(de_data,  padj <= cutoff)
foldchange <- subset(de_data, log2FoldChange >= 2 | log2FoldChange <= -2)

ggplot(de_data, aes(log(baseMean,2), log2FoldChange)) +
  geom_point(color="gray20", alpha = 0.5) + ##change opacity of points
  theme_classic() + ##use the classic theme template
  xlab("log2 Average Normalized Counts") + ##label the x-axis
  ylab("log2 Fold Change") + ##label the y-axis
  ggtitle(paste0('MA plot')) + ##to add title
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 8)) + ##to center title
  geom_hline(aes(yintercept=0), colour="blue", linetype="dashed") + ##plot a horizontal line
  geom_point(data=signifROI,aes(log(baseMean,2), log2FoldChange), color="red",size=1, alpha=0.75)
ggsave("MA_plot.png", plot=last_plot(), path=outdir)

ggplot(de_data, aes(log2FoldChange, -log(padj))) +
  geom_point() +
  theme_classic() + ##use the classic theme template
  xlab("log2 Fold Change") + ##label the x-axis
  ylab("-log10 Adjusted P-value") + ##label the y-axis
  ggtitle(paste0('Volcano plot')) + ##to add title
  theme(plot.title = element_text(hjust = 0.5)) + ##to center title
  geom_point(data=signifROI,aes(log2FoldChange, -log(padj)), color="red",size=1.5, alpha=0.5) +
  geom_point(data=foldchange,aes(log2FoldChange, -log(padj)), color="orange",size=1.5, alpha=0.25)
ggsave("Volcano_plot.png", plot=last_plot(), path=outdir)
