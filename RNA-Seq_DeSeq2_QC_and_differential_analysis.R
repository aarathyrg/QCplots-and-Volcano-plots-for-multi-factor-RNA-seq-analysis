## RNA-seq analysis with DESeq2 (https://doi.org/doi:10.18129/B9.bioc.DESeq2)
## Aarathy R.G.

# RNA-seq data from GSE115435
# RNA-seq on Interferon (type I and type II )in WT and IRF9KO BMDMs

#count matrix obtained after running nf-core/rnaseq pipeline.


# set working directory
getwd()
setwd("/Users/aarathyrg/Dropbox/Aarathy_katrin_proseq/")

library("DESeq2")
library("dplyr")
library("RColorBrewer")
library("gplots")
#genes are rows, samples columns
#set work directory
setwd("/Users/aarathyrg/Dropbox/")
countdata <- read.table("salmon.merged.gene_counts.tsv",header = T,row.names = 1)
countdata <- as.data.frame(countdata)

# Remove first columns gene id
#rownames are ensemble ids.

countdata <- countdata[ ,2:ncol(countdata)]
names <- names(countdata)

# Convert to matrix
countdata <- as.matrix(countdata)

#Creating coldata based on 
##rep function repeats the string to the number provided
genotype <- factor(c(rep("ko", 9), rep("wt", 9)))
treatment <- factor(c(rep("beta", 3),rep("gamma", 3),rep("ut", 3),rep("beta", 3),rep("gamma", 3),rep("ut", 3)))
#combining genotype and treatment
condition <- factor(c(rep("ko_beta",3),rep("ko_gamma", 3),rep("ko_ut", 3),rep("wt_beta",3),rep("wt_gamma", 3),rep("wt_ut", 3)))

# Create a coldata dataframe and DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), genotype, treatment,condition)

#######################################
#setting design as condition###########
#######################################
#combined genotype and treatment tomake another factor -condition.
#alternative use intersect(check ?DESeqDataSetFromMatrix)
dds <- DESeqDataSetFromMatrix(countData=round(countdata), colData=coldata, design=~condition)
dds$condition

#####################################
#####*****QC plots*****####
###################################

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
rld_mat<-assay(rld)
###################################
## Sample correlation plot########
###################################
sampleDists <- dist( t( rld_mat ) ,method = "euclidean")
sampleDistMatrix <- as.matrix(sampleDists)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#hierarchical clustering
hc <- hclust(sampleDists)
png("qc-sample_correlation.png", w=1000, h=1000, pointsize=20)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )
dev.off()

###Principal components analysis########

plotPCA(rld, intgroup = c("genotype","treatment"))
png("qc-pca.png", w=1000, h=1000, pointsize=20)
plotPCA(rld,intgroup = c("genotype","treatment"))
dev.off()

##set the reference level*****IMPORTANT***###set the control to which you compare

dds$conditiion <- relevel(dds$conditiion,"wt_ut")


###*************************#####
# Run DESeq2 
###************************######
###
?results# check adjustable parameters
cooksCutoff <- TRUE                   
independentFiltering <- TRUE          
alpha <- 0.05                     
pAdjustMethod <- "BH"               


dds <-DESeq(dds)

# Get differential expression results
res <- results(dds,alpha = 0.05)
table(results$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
datares<-as.data.frame(res)
summary(res)

###########################################################
####save all comparisons you want in seperate tables###
############################################################
#to make all contrasts vs reference level
compare<-levels(condition)[1:length(levels(condition))-1]#excluding reference level to calculate length
table<-list()
for (i in 1:length(compare)){
  table[i]<-results(dds, alpha = 0.05, contrast = c("condition",compare[i],"wt_ut"))
  write.table(table[i],paste("/Users/aarathyrg/Desktop/CODES/",
                             compare[i],sep = ""),sep = "\t",col.names = T,row.names = T)
}

