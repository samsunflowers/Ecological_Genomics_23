## Set your working directory; change to your working directory!
setwd("C:/Users/bjork/Documents/GitHub/Ecological_Genomics_23/PopGenomics/results/")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(eulerr)

# Import the counts matrix for Hudsonica
countsTableHudsonica <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTableHudsonica)

# Import the counts matrix for Tonsa
countsTableTonsaF1 <- read.table("DE_counts_F1.txt", header=TRUE, row.names=1)
head(countsTableTonsaF1)
countsTableTonsaF3 <- read.table("DE_counts_F3.txt", header=TRUE, row.names=1)
head(countsTableTonsaF3)

# Round the counts matrix for Hudsonica because DESeq2 does not like decimals
countsTableRoundHudsonica <- round(countsTableHudsonica)
head(countsTableRoundHudsonica)

# Round the counts matrix for Tonsa because DESeq2 does not like decimals
countsTableRoundTonsaF1 <- round(countsTableTonsaF1)
head(countsTableRoundTonsaF1)
countsTableRoundTonsaF3 <- round(countsTableTonsaF3)
head(countsTableRoundTonsaF3)

# Combining Tonsa F1 and F3 rounded counts data
common_row_names <- intersect(row.names(countsTableRoundTonsaF1), row.names(countsTableRoundTonsaF3))
countsTableRoundTonsaF1_common <- countsTableRoundTonsaF1[common_row_names, ]
countsTableRoundTonsaF3_common <- countsTableRoundTonsaF3[common_row_names, ]
countsTableRoundTonsa <- cbind(countsTableRoundTonsaF1_common, countsTableRoundTonsaF3_common)
head(countsTableRoundTonsa)

# Import the sample description table for Hudsonica and Tonsa
condsHudsonica <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(condsHudsonica)
condsTonsaF1 <- read.delim("RT_tonsa_F1_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
condsTonsaF3 <- read.delim("RT_tonsa_F3_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)

# Add a column that states the generation for each Tonsa dataset to make it easier to subset later
condsTonsaF1$generation <- c("F1")
condsTonsaF3$generation <- c("F3")

# Combine the condition data for F1 and F3 Tonsa data
condsTonsa <- rbind(condsTonsaF1, condsTonsaF3)
colnames(condsTonsa) <- c("treatment", "line", "environment", "generation")
head(condsTonsa)

# How many reads do we have in Hudsonica?
colSums(countsTableRoundHudsonica)
mean(colSums(countsTableRoundHudsonica))
barplot(colSums(countsTableRoundHudsonica), names.arg=colnames(countsTableRoundHudsonica),cex.names=0.5, las=3,ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRoundHudsonica)), col="blue", lwd=2)

# How many reads do we have in Tonsa?
colSums(countsTableRoundTonsa)
mean(colSums(countsTableRoundTonsa))
barplot(colSums(countsTableRoundTonsa), names.arg=colnames(countsTableRoundTonsa),cex.names=0.5, las=3,ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRoundTonsa)), col="blue", lwd=2)

# What is the average number of counts per gene in Hudsonica?
rowSums(countsTableRoundHudsonica)
mean(rowSums(countsTableRoundHudsonica)) # 8217.81
median(rowSums(countsTableRoundHudsonica)) # 377
apply(countsTableRoundHudsonica,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRoundHudsonica,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRoundHudsonica,1,mean),xlim=c(0,10000), ylim=c(0,60000),breaks=1000)

# What is the average number of counts per gene in Tonsa?
rowSums(countsTableRoundTonsa)
mean(rowSums(countsTableRoundTonsa)) # 24131.74
median(rowSums(countsTableRoundTonsa)) # 4776
apply(countsTableRoundTonsa,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRoundTonsa,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRoundTonsa,1,mean),xlim=c(0,10000), ylim=c(0,60000),breaks=1000)

# Create a DESeq object and define the experimental design here with the tilda
# HUDSONICA
ddsHudsonica <- DESeqDataSetFromMatrix(countData = countsTableRoundHudsonica, colData=condsHudsonica,
                                       design= ~ treatment + generation)
dim(ddsHudsonica)

# TONSA
ddsTonsa <- DESeqDataSetFromMatrix(countData = countsTableRoundTonsa, colData=condsTonsa, 
                                   design= ~ treatment + generation)
dim(ddsTonsa)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28) suggested by WGCNA on RNAseq FAQ
# HUDSONICA
ddsHudsonica <- ddsHudsonica[rowSums(counts(ddsHudsonica) >= 15) >= 28,]
nrow(ddsHudsonica) # 25260, that have at least 15 reads (aka counts) in 75% of the samples

# TONSA
ddsTonsa <- ddsTonsa[rowSums(counts(ddsTonsa) >= 15) >= 28,]
nrow(ddsTonsa) # 21353, that have at least 15 reads (aka counts) in 75% of the samples

# Subsetting data for F0 and F4
ddsHudsonica <- subset(ddsHudsonica, select = generation == "F0" | generation =="F4")
ddsHudsonica <- subset(ddsHudsonica, select = treatment == "OWA" | treatment =="AM")
dim(ddsHudsonica)

# Run the DESeq model to test for differential gene expression
ddsHudsonica$generation <- droplevels(ddsHudsonica$generation)
ddsHudsonica$treatment <- droplevels(ddsHudsonica$treatment)
ddsHudsonica <- DESeq(ddsHudsonica)
ddsTonsa <- DESeq(ddsTonsa)

# List the results you've generated
# HUDSONICA
resultsNames(ddsHudsonica)
# [1] "Intercept"           "treatment_OWA_vs_AM" "generation_F4_vs_F0"

# TONSA
resultsNames(ddsTonsa)
# [1] "Intercept"                  "treatment_AAtoHH_vs_AAtoAA" "treatment_HHtoAA_vs_AAtoAA"
# [4] "treatment_HHtoHH_vs_AAtoAA" "generation_F3_vs_F1"  

# Check the quality of the data by sample clustering and visualization
# The goal of transformation "is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low."
ntdHudsonica <- normTransform(ddsHudsonica)
meanSdPlot(assay(ntdHudsonica))

ntdTonsa <- normTransform(ddsTonsa)
meanSdPlot(assay(ntdTonsa))

# Variance stabilizing transformation
vsdHudsonica <- vst(ddsHudsonica, blind=FALSE)
meanSdPlot(assay(vsdHudsonica))
sampleDistsHudsonica <- dist(t(assay(vsdHudsonica)))
sampleDistMatrixHudsonica <- as.matrix(sampleDistsHudsonica)
rownames(sampleDistMatrixHudsonica) <- paste(vsdHudsonica$treatment, vsdHudsonica$generation, sep="-")
colnames(sampleDistMatrixHudsonica) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixHudsonica,
         clustering_distance_rows=sampleDistsHudsonica,
         clustering_distance_cols=sampleDistsHudsonica,
         col=colors)

vsdTonsa <- vst(ddsTonsa, blind=FALSE)
meanSdPlot(assay(vsdTonsa))
sampleDistsTonsa <- dist(t(assay(vsdTonsa)))
sampleDistMatrixTonsa <- as.matrix(sampleDistsTonsa)
rownames(sampleDistMatrixTonsa) <- paste(vsdTonsa$treatment, vsdTonsa$generation, sep="-")
colnames(sampleDistMatrixTonsa) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixTonsa,
         clustering_distance_rows=sampleDistsTonsa,
         clustering_distance_cols=sampleDistsTonsa,
         col=colors)

# PCA to visualize global gene expression patterns
# First transform the data for plotting using variance stabilization
vsdHudsonica <- vst(ddsHudsonica, blind=FALSE)
pcaDataHudsonica <- plotPCA(vsdHudsonica, intgroup=c("treatment","generation"), returnData=TRUE)
percentVarHudsonica <- round(100 * attr(pcaDataHudsonica,"percentVar"))
ggplot(pcaDataHudsonica, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarHudsonica[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarHudsonica[2],"% variance")) + 
  coord_fixed()

vsdTonsa <- vst(ddsTonsa, blind=FALSE)
pcaDataTonsa <- plotPCA(vsdTonsa, intgroup=c("treatment","generation"), returnData=TRUE)
percentVarTonsa <- round(100 * attr(pcaDataTonsa,"percentVar"))
ggplot(pcaDataTonsa, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarTonsa[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarTonsa[2],"% variance")) + 
  coord_fixed()

## Check on the DE results from the DESeq command way above
resultsNames(ddsHudsonica)
resHudsF4_F0<- results(ddsHudsonica, name="generation_F4_vs_F0", alpha=0.05)
resHudsF4_F0 <- resHudsF4_F0[order(resHudsF4_F0$padj),]
head(resHudsF4_F0)  
summary(resHudsF4_F0)

resultsNames(ddsTonsa)
resTonsaF3_F1<- results(ddsTonsa, name="generation_F3_vs_F1", alpha=0.05)
resTonsaF3_F1 <- resTonsaF3_F1[order(resTonsaF3_F1$padj),]
head(resTonsaF3_F1)  
summary(resTonsaF3_F1)

# We can make an MA plot for Hudsonica and Tonsa
plotMA(resHudsF4_F0, ylim=c(-4,4))
plotMA(resTonsaF3_F1, ylim=c(-4,4))

# Heatmap of top 20 genes sorted by pvalue
vsdHuds <- vst(ddsHudsonica, blind=FALSE)
topgenesHuds <- head(rownames(resHudsF4_F0),10)
matHuds <- assay(vsdHuds)[topgenesHuds,]
matHuds <- matHuds - rowMeans(matHuds)
dfHuds <- as.data.frame(colData(ddsHudsonica)[,c("generation","treatment")])
pheatmap(matHuds, annotation_col=dfHuds)
pheatmap(matHuds, annotation_col=dfHuds, cluster_cols = F)

vsdTons <- vst(ddsTonsa, blind=FALSE)
topgenesTons <- head(rownames(resTonsaF3_F1),20)
matTons <- assay(vsdTons)[topgenesTons,]
matTons <- matTons - rowMeans(matTons)
dfTons <- as.data.frame(colData(ddsTonsa)[,c("generation","treatment")])
pheatmap(matTons, annotation_col=dfTons)
pheatmap(matTons, annotation_col=dfTons, cluster_cols = F)

# Revist the Tonsa data
ddsTonsa <- DESeqDataSetFromMatrix(countData = countsTableRoundTonsa, colData=condsTonsa, 
                             design= ~ treatment)
dim(ddsTonsa)

# Filter 
ddsTonsa <- ddsTonsa[rowSums(counts(ddsTonsa) >= 15) >= 28,]
nrow(ddsTonsa) 

# Subset the DESeqDataSet to F1 generation
dds_Tonsa_F1 <- subset(ddsTonsa, select = generation == 'F1')
dim(dds_Tonsa_F1)

# Perform DESeq2 analysis on the subset
dds_Tonsa_F1 <- DESeq(dds_Tonsa_F1)
resultsNames(dds_Tonsa_F1)
res_Tonsa_F1_HHHHvsAAAA <- results(dds_Tonsa_F1, name="treatment_HHtoHH_vs_AAtoAA", alpha=0.05)
res_Tonsa_F1_HHHHvsAAAA <- res_Tonsa_F1_HHHHvsAAAA[order(res_Tonsa_F1_HHHHvsAAAA$padj),]
head(res_Tonsa_F1_HHHHvsAAAA)
summary(res_Tonsa_F1_HHHHvsAAAA)
res_Tonsa_F1_HHHHvsAAAA <- res_Tonsa_F1_HHHHvsAAAA[!is.na(res_Tonsa_F1_HHHHvsAAAA$padj),]
degs_Tonsa_F1_HHHHvsAAAA <- row.names(res_Tonsa_F1_HHHHvsAAAA[res_Tonsa_F1_HHHHvsAAAA$padj < 0.05,])
length(degs_Tonsa_F1_HHHHvsAAAA)  # 1139

# Subset the DESeqDataSet to F1 generation
dds_Tonsa_F3<- subset(ddsTonsa, select = generation == 'F3')
dim(dds_Tonsa_F3)

# Perform DESeq2 analysis on the subset
dds_Tonsa_F3 <- DESeq(dds_Tonsa_F3)
resultsNames(dds_Tonsa_F3)
res_Tonsa_F3_HHHHvsAAAA <- results(dds_Tonsa_F3, name="treatment_HHtoHH_vs_AAtoAA", alpha=0.05)
res_Tonsa_F3_HHHHvsAAAA <- res_Tonsa_F3_HHHHvsAAAA[order(res_Tonsa_F3_HHHHvsAAAA$padj),]
head(res_Tonsa_F3_HHHHvsAAAA)
summary(res_Tonsa_F3_HHHHvsAAAA)
res_Tonsa_F3_HHHHvsAAAA <- res_Tonsa_F3_HHHHvsAAAA[!is.na(res_Tonsa_F3_HHHHvsAAAA$padj),]
degs_Tonsa_F3_HHHHvsAAAA <- row.names(res_Tonsa_F3_HHHHvsAAAA[res_Tonsa_F3_HHHHvsAAAA$padj < 0.05,])
length(degs_Tonsa_F3_HHHHvsAAAA)  # 453

# Intersections
length(intersect(degs_Tonsa_F3_HHHHvsAAAA,degs_Tonsa_F1_HHHHvsAAAA))  # 73

# Tonsa venn diagram
# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("F1" = 1193, "F3" = 453, "F1&F3" = 73))
plot(fit1,  lty = 1:3, quantities = list(cex=2))
plot2 <- plot(fit1, quantities = list(cex=2), fill = wes_palette("GrandBudapest1"),
              lty = 1:3,
              labels = list(font = 25))

plot2

annotate_figure(plot2, top = text_grob("Shared Generational Genes in OWA vs AM in Tonsa",
                                       color = "black", face = "bold", size = 12))

# Revist the Hudsonica data
ddsHudsonica <- DESeqDataSetFromMatrix(countData = countsTableRoundTonsa, colData=condsTonsa, 
                                   design= ~ treatment)
dim(ddsHudsonica)

# Filter 
ddsHudsonica <- ddsHudsonica[rowSums(counts(ddsHudsonica) >= 15) >= 28,]
nrow(ddsHudsonica) 

# Subset the DESeqDataSet to F0 generation
dds_Hudsonica_F0 <- subset(ddsHudsonica, select = generation == 'F0')
dim(dds_Hudsonica_F0)

# Perform DESeq2 analysis on the subset
dds_Hudsonica_F0 <- DESeq(dds_Hudsonica_F0)
resultsNames(dds_Hudsonica_F0)
res_Hudsonica_F0_OWAvsAM <- results(dds_Hudsonica_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_Hudsonica_F0_OWAvsAM <- res_Hudsonica_F0_OWAvsAM[order(res_Hudsonica_F0_OWAvsAM$padj),]
head(res_Hudsonica_F0_OWAvsAM)
summary(res_Hudsonica_F0_OWAvsAM)
res_Hudsonica_F0_OWAvsAM <- res_Hudsonica_F0_OWAvsAM[!is.na(res_Hudsonica_F0_OWAvsAM$padj),]
degs_Hudsonica_F0_OWAvsAM <- row.names(res_Hudsonica_F0_OWAvsAM[res_Hudsonica_F0_OWAvsAM$padj < 0.05,])
length(degs_Hudsonica_F0_OWAvsAM)  # 1139

# Subset the DESeqDataSet to F4 generation
dds_Hudsonica_F4<- subset(ddsHudsonica, select = generation == 'F4')
dim(dds_Hudsonica_F4)

# Perform DESeq2 analysis on the subset
dds_Hudsonica_F4 <- DESeq(dds_Hudsonica_F4)
resultsNames(dds_Hudsonica_F4)
res_Hudsonica_F4_OWAvsAM <- results(dds_Hudsonica_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_Hudsonica_F4_OWAvsAM <- res_Hudsonica_F4_OWAvsAM[order(res_Hudsonica_F4_OWAvsAM$padj),]
head(res_Hudsonica_F4_OWAvsAM)
summary(res_Hudsonica_F4_OWAvsAM)
res_Hudsonica_F4_OWAvsAM <- res_Hudsonica_F4_OWAvsAM[!is.na(res_Hudsonica_F4_OWAvsAM$padj),]
degs_Hudsonica_F4_OWAvsAM <- row.names(res_Hudsonica_F4_OWAvsAM[res_Hudsonica_F4_OWAvsAM$padj < 0.05,])
length(degs_Hudsonica_F4_OWAvsAM)  # 453

# Intersections
length(intersect(degs_Hudsonica_F4_OWAvsAM,degs_Hudsonica_F0_OWAvsAM))  # 73

# Hudsonica venn diagram
# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("F0" = 3637, "F4" = 171, "F0&F4" = 55))
plot(fit1,  lty = 1:3, quantities = list(cex = 2))
plot2 <- plot(fit1, quantities = list(cex=2), fill = wes_palette("GrandBudapest1"),
              lty = 1:3,
              labels = list(font = 25))

annotate_figure(plot2, top = text_grob("Shared Generational Genes in OWA vs AM in Hudsonica",
                                       color = "black", face = "bold", size = 12))

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_Tonsa_F1_HHHHvsAAAA_df <- data.frame(transcriptID = rownames(res_Tonsa_F1_HHHHvsAAAA), res_Tonsa_F1_HHHHvsAAAA)
res_Tonsa_F3_HHHHvsAAAA_df <- data.frame(transcriptID = rownames(res_Tonsa_F3_HHHHvsAAAA), res_Tonsa_F3_HHHHvsAAAA)
res_Hudsonica_F0_OWAvsAM_df <- data.frame(transcriptID = rownames(res_Hudsonica_F0_OWAvsAM), res_Hudsonica_F0_OWAvsAM)
res_Hudsonica_F4_OWAvsAM_df <- data.frame(transcriptID = rownames(res_Hudsonica_F4_OWAvsAM), res_Hudsonica_F4_OWAvsAM)

# Split the "transcriptID" column by double colons and create new columns of the parts
res_Tonsa_F1_HHHHvsAAAA_df <- separate(res_Tonsa_F1_HHHHvsAAAA_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 
res_Tonsa_F3_HHHHvsAAAA_df <- separate(res_Tonsa_F3_HHHHvsAAAA_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 
res_Hudsonica_F0_OWAvsAM_df <- separate(res_Hudsonica_F0_OWAvsAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 
res_Hudsonica_F4_OWAvsAM_df <- separate(res_Hudsonica_F4_OWAvsAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 


# Create a new column by concatenating "part1" and "part2" with double colons in between
res_Tonsa_F1_HHHHvsAAAA_df$transcriptID_trim <- paste(res_Tonsa_F1_HHHHvsAAAA_df$part1, res_Tonsa_F1_HHHHvsAAAA_df$part2, sep = "::")
res_Tonsa_F3_HHHHvsAAAA_df$transcriptID_trim <- paste(res_Tonsa_F3_HHHHvsAAAA_df$part1, res_Tonsa_F3_HHHHvsAAAA_df$part2, sep = "::")
res_Hudsonica_F0_OWAvsAM_df$transcriptID_trim <- paste(res_Hudsonica_F0_OWAvsAM_df$part1, res_Hudsonica_F0_OWAvsAM_df$part2, sep = "::")
res_Hudsonica_F4_OWAvsAM_df$transcriptID_trim <- paste(res_Hudsonica_F4_OWAvsAM_df$part1, res_Hudsonica_F4_OWAvsAM_df$part2, sep = "::")

# Select the two columns we want to save for the GOMWU analysis
selected_res_Tonsa_F1_HHHHvsAAAA_df <- res_Tonsa_F1_HHHHvsAAAA_df[c("transcriptID_trim", "log2FoldChange")]
selected_res_Tonsa_F3_HHHHvsAAAA_df <- res_Tonsa_F3_HHHHvsAAAA_df[c("transcriptID_trim", "log2FoldChange")]
selected_res_Hudsonica_F0_OWAvsAM_df <- res_Hudsonica_F0_OWAvsAM_df[c("transcriptID_trim", "log2FoldChange")]
selected_res_Hudsonica_F4_OWAvsAM_df <- res_Hudsonica_F4_OWAvsAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file
write.csv(selected_res_Tonsa_F1_HHHHvsAAAA_df, file = "res_Tonsa_F1_HHHHvsAAAA_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(selected_res_Tonsa_F3_HHHHvsAAAA_df, file = "res_Tonsa_F3_HHHHvsAAAA_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(selected_res_Hudsonica_F0_OWAvsAM_df, file = "res_Hudsonica_F0_OWAvsAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(selected_res_Hudsonica_F4_OWAvsAM_df, file = "res_Hudsonica_F4_OWAvsAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU

# Save the selected columns WITH P VAL as a CSV file
write.csv(res_Tonsa_F1_HHHHvsAAAA_df, file = "res_pval_Tonsa_F1_HHHHvsAAAA_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(res_Tonsa_F3_HHHHvsAAAA_df, file = "res_pval_Tonsa_F3_HHHHvsAAAA_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(res_Hudsonica_F0_OWAvsAM_df, file = "res_pval_Hudsonica_F0_OWAvsAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(res_Hudsonica_F4_OWAvsAM_df, file = "res_pval_Hudsonica_F4_OWAvsAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU

