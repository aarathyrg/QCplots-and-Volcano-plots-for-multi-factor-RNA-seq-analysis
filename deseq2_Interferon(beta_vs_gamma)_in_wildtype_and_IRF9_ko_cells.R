## RNA-seq analysis with DESeq2
## Aarathy, Thomas Decker group, Max Perutz Labs

# RNA-seq data from GSE115435
# RNA-seq on Interferon (type I and type II )in WT and IRF9KO BMDMs

#count matrix obtained after running nf-core/rnaseq pipeline.


# set working directory
getwd()
setwd("Dropbox/Aarathy_katrin_proseq/")

library("DESeq2")
library("dplyr")

countdata <- read.table("salmon.merged.gene_counts.tsv",header = T,row.names = 1)
head(countdata)

countdata <- as.data.frame(countdata)

# Remove first columns ensemble id
#rownames are ensemble ids.

countdata <- countdata[ ,2:ncol(countdata)]
names <- names(countdata)
nrow(countdata)

# Convert to matrix
countdata <- as.matrix(countdata)

head(countdata)

#Creating coldata
##rep function repeats the string to the number provided
genotype <- factor(c(rep("ko", 9), rep("wt", 9)))
treatment <- factor(c(rep("beta", 3),rep("gamma", 3),rep("ut", 3),rep("beta", 3),rep("gamma", 3),rep("ut", 3)))
condition <- factor(c(rep("ko_beta",3),rep("ko_gamma", 3),rep("ko_ut", 3),rep("wt_beta",3),rep("wt_gamma", 3),rep("wt_ut", 3)))
# Analysis with DESeq2 ----------------------------------------------------library(DESeq2)
##DESeq2 need a file with raw read counts, 
####*a ColData files describing the factors:rownames must match  coulnames of counts file 

# Create a coldata dataframe and DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), genotype, treatment,condition)
coldata 
all(rownames(coldata) == colnames(countdata))

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
rld<-assay(rld)


###################################
## Sample correlation plot########
###################################
# Colors for plots below
library(RColorBrewer)
(mycols <- brewer.pal(8, "Greys")[1:length(unique(condition))])
sampleDists <- dist( t( assay(rld) ) )
sampleDists

library(gplots)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$genotype, rld$treatment, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#hierarchical clustering
#alternative kmeans
hc <- hclust(sampleDists)
png("qc-heatmap-samples_blue_hc.png", w=1000, h=1000, pointsize=20)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )
dev.off()

###Principal components analysis########
## DESeq2::plotPCA(rld, intgroup="condition")
plotPCA(rld, intgroup = c("genotype","treatment"))
png("qc-pca.png", w=1000, h=1000, pointsize=20)
plotPCA(rld,intgroup = c("genotype","treatment"))
dev.off()

##set the reference level*****IMPORTANT***###set the control to which you compare

dds$conditiion <- relevel(dds$conditiion,"wt_ut")
colData(dds)
#################################
###*************************#####
# Run the DESeq pipeline
###************************######
#################################
dds <-DESeq(dds)

# Get differential expression results
res <- results(dds,alpha = 0.05)
table(results$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
datares<-as.data.frame(res)
summary(res)

##############################
## Merge with rlog transformed counts
rld_df<-as.data.frame(rld)
results_deseq2_all_rlog_counts <- merge(as.data.frame(res), rld,  by="row.names", sort=FALSE)
#######################################
#Obtain_gene_names
countdata_with_genename <- read.table("salmon.merged.gene_counts.tsv",header = T,row.names = 1)
gene_id_name<-countdata_with_genename["gene_name"]
head(gene_id_name)
####save all comparisons you want in seperate tables###
####***wt_beta***###

wt_beta <-results(dds, alpha = 0.05, contrast = c("condition","wt_beta","wt_ut"))
wt_beta<-merge(as.data.frame(wt_beta), as.data.frame(gene_id_name), by="row.names", sort = FALSE)
write.csv(wt_beta,"res_wt_beta.csv")
####***wt_gamma***###
wt_gamma <- results(dds,alpha = 0.05, contrast = c("condition","wt_gamma","wt_ut"))
wt_gamma<-merge(as.data.frame(wt_gamma), as.data.frame(gene_id_name), by="row.names", sort = FALSE)
write.csv(wt_gamma,"res_wt_gamma.csv")
####***ko_beta***###
ko_beta <-results(dds, alpha = 0.05, contrast = c("condition","ko_beta","wt_beta"))
ko_beta<-merge(as.data.frame(ko_beta), as.data.frame(gene_id_name), by="row.names", sort = FALSE)
write.csv(ko_beta,"res_ko_beta.csv")

####***ko_gamma***###
ko_wt_gamma <-results(dds,alpha = 0.05, contrast =c("condition","ko_gamma","wt_gamma"))
ko_gamma<-merge(as.data.frame(ko_gamma), as.data.frame(gene_id_name), by="row.names", sort = FALSE)
write.csv(ko_gamma,"res_ko_to_wt_gamma.csv")
########******ChIP annotated peaks******######
#annotated peak files from ChIP-seq data
###############################################
chip_peaks_ifnb_promoters_irf9 <- read.delim("WT_IRF9_IFNB_R1_peaks.annotatePeaks_PROMOTER_FILTERED.txt",header = T)
###############################################
#########################
###**Volcano plot***#####
#########################
library("ggplot2")
library("dplyr")
#rfiltering na values

wt_beta_volcano <- wt_beta %>% filter(! is.na(log2FoldChange))%>% filter(! is.na(padj))
wt_beta_volcano_down <- wt_beta_volcano%>% filter((log2FoldChange< -1))%>% filter((padj<0.05))
ko_beta_up_volcano <- ko_beta %>% filter(! is.na(log2FoldChange))%>% filter(! is.na(padj))%>% filter((log2FoldChange>1))%>% filter((padj<0.05))
ko_beta_down_volcano <- ko_beta %>% filter(! is.na(log2FoldChange))%>% filter(! is.na(padj))%>% filter((log2FoldChange< -1))%>% filter((padj<0.05))

wt_beta_volcano["group"] <- "a:n.s."
wt_beta_volcano[which(wt_beta_volcano['padj'] <= 0.05 & wt_beta_volcano['log2FoldChange'] >= 1.0 ),"group"] <- "b:Wt_IFNb_up"
wt_beta_volcano[which(wt_beta_volcano['padj'] <= 0.05 & wt_beta_volcano['log2FoldChange'] <= -1.0 ),"group"] <- "c:Wt_IFNb_down"
wt_beta_volcano[which(wt_beta_volcano$Row.names %in% ko_beta_down_volcano$Row.names & 
                        wt_beta_volcano['padj'] <= 0.05 & wt_beta_volcano['log2FoldChange'] >= 1.0 ),"group"] <- "d:IFNbUp_irf9_upregulated"
wt_beta_volcano[which(wt_beta_volcano$Row.names %in% ko_beta_up_volcano$Row.names &  
                        wt_beta_volcano['padj'] <= 0.05 & wt_beta_volcano['log2FoldChange'] >= 1.0 ),"group"] <- "e:IFNbUp_irf9_downregulated"
wt_beta_volcano[which(wt_beta_volcano$Row.names %in% ko_beta_down_volcano$Row.names  & 
                        wt_beta_volcano$Row.names %in% chip_peaks_ifnb_promoters_irf9$Nearest.PromoterID &
                        wt_beta_volcano['padj'] <= 0.05 & wt_beta_volcano['log2FoldChange'] >= 1.0 ),"group"] <- "f:IFNbUp_irf9_upregulated_IRF9bound"

wt_beta_volcano[which(wt_beta_volcano$Row.names %in% ko_beta_down_volcano$Row.names &  wt_beta_volcano['padj'] <= 0.05 &
                        wt_beta_volcano['log2FoldChange'] <= -1.0 ),"group"] <- "g:IFNbdown_irf9_upregulated"
wt_beta_volcano[which(wt_beta_volcano$Row.names %in% ko_beta_up_volcano$Row.names &  wt_beta_volcano['padj'] <= 0.05 &
                        wt_beta_volcano['log2FoldChange'] <= -1.0 ),"group"] <- "h:IFNbdown_irf9_downregulated"
wt_beta_volcano[which(wt_beta_volcano$Row.names %in% ko_beta_up_volcano$Row.names &
                        wt_beta_volcano$Row.names %in% chip_peaks_ifnb_promoters_irf9$Nearest.PromoterID & wt_beta_volcano['padj'] <= 0.05 &
                        wt_beta_volcano['log2FoldChange'] <= -1.0 ),"group"] <- "i:IFN_down_irf9_downregulated_irf9_bound"

################################################################
###***Volcano
##########################################################
#interactive plot from plotly
##+++++++++plotly++++++++++++###
###########################

library(plotly)


p <- plot_ly(data = wt_beta_volcano, x = wt_beta_volcano$log2FoldChange, y = -log10(as.numeric(wt_beta_volcano$padj)), text = wt_beta_volcano$gene_name, mode = "markers", color = wt_beta_volcano$group, colors = c("#5F5C5D","#2F5F6D",
                                                                                                                                                                                                                     "#77004F","#AF2804","#7F74AD","#77B7A5","#260077","#3181CE","#BFA5E0")) %>%
  layout(title = "Interferon_up_downregulated.pdf",
         xaxis = list(title = "log2(Fold change)"),
         yaxis = list(title = "-log10(padj)"))
p

htmlwidgets::saveWidget(as_widget(p), "plotly_ifnb.html")
##############################################
#Saving tables

Up_genes<-as.data.frame(wt_beta_volcano %>% filter((log2FoldChange>1))%>% filter((padj<0.05)))

write.table(Up_genes,"up_genes_IFNb_geneid.tsv", sep = "\t",row.names = F,col.names = T)

down_genes<-as.data.frame(wt_beta_volcano %>% filter((log2FoldChange< -1))%>% filter((padj<0.05)))
#irf9 as upregulator
Up_genes_irf9_upregulated <- as.data.frame(wt_beta_volcano %>% filter((log2FoldChange> 1))%>% filter((padj<0.05))%>%filter(Row.names %in% ko_beta_down_volcano$Row.names ))
#irf9 as downregulator
Down_genes_irf9_downregulated <-as.data.frame(wt_beta_volcano %>% filter((log2FoldChange<= -1))%>% filter((padj<0.05))%>%filter(Row.names %in% ko_beta_up_volcano$Row.names ))
#Considering IRF9 binding from ChIP -seq data

Down_genes_irf9_downregulated_IRF9_bound<-as.data.frame(subset(wt_beta_volcano, group == "i:IFN_down_irf9_downregulated_irf9_bound"))
Down_genes_irf9_downregulated_IRF9_nonbound<-as.data.frame(subset(wt_beta_volcano, group == "h:IFN_down_irf9_downregulated"))












