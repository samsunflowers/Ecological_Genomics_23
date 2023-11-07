## Set your working directory
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

# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] 8217.81
median(rowSums(countsTableRound)) # [1] 377
apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,10000), ylim=c(0,60000),breaks=1000)

## Create a DESeq object and define the experimental design here with the tilda
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment + generation)
dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28) suggested by WGCNA on RNAseq FAQ
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) # 25260, that have at least 15 reads (aka counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
# [1] "Intercept"            "treatment_OA_vs_AM"   "treatment_OW_vs_AM"  
# [4] "treatment_OWA_vs_AM"  "generation_F11_vs_F0" "generation_F2_vs_F0" 
# [7] "generation_F4_vs_F0"

###############################################################

# Check the quality of the data by sample clustering and visualization
# The goal of transformation "is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low."
# This gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
# AM_F11_Rep3 and OA_F2_Rep2 could be outliers

###############################################################
# PCA to visualize global gene expression patterns
# First transform the data for plotting using variance stabilization
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData,"percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# Strong gene response in the first generation and changed responses by the fourth through eleventh generation.
# Let's plot the PCA by generation in four panels
data <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

###############################################################
dataF0 <- subset(data, generation == 'F0')
F0 <- ggplot(dataF0, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-10, 25) + xlim(-40, 10)+ # zoom for F0 with new assembly
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "Ambient & Warming"))+
  scale_fill_manual(values=wes_palette("GrandBudapest1"), labels = c("Ambient", "Acidification","Warming", "Ambient & Warming"))+
  theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 4)))+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 2))+
  theme(text = element_text(size = 15)) +
  theme(legend.title = element_blank()) +
  ggtitle("F0 Gene Expression")
F0

###############################################################

dataF2 <- subset(data, generation == 'F2')
F2 <- ggplot(dataF2, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+
  scale_shape_manual(values=c(21,22,23), labels = c("Ambient", "Acidification","Warming"))+
  scale_fill_manual(values=wes_palette("GrandBudapest1"), labels = c("Ambient", "Acidification","Warming"))+
  theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(shape = guide_legend(override.aes = list(size = 4)))+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 2))+
  theme(text = element_text(size = 15)) +
  theme(legend.title = element_blank()) +
  ggtitle("F2 Gene Expression")
F2
# F2 is missing one ambient replicate

###############################################################

dataF4 <- subset(data, generation == 'F4')
F4 <- ggplot(dataF4, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+ # limits with filtered assembly
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "Ambient & Warming"))+
  scale_fill_manual(values=wes_palette("GrandBudapest1"), labels = c("Ambient", "Acidification","Warming", "Ambient & Warming"))+
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 4)))+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 2))+
  theme(text = element_text(size = 15)) +
  theme(legend.title = element_blank()) +
  ggtitle("F4 Gene Expression")
F4

###############################################################

dataF11 <- subset(data, generation == 'F11')
F11 <- ggplot(dataF11, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-45, 25) + xlim(-50, 55)+
  scale_shape_manual(values=c(21,24), labels = c("Ambient", "Ambient & Warming"))+
  scale_fill_manual(values=wes_palette("GrandBudapest1"), labels = c("Ambient", "Ambient & Warming"))+
  guides(shape = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 4)))+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 2))+
  theme(text = element_text(size = 15)) +
  theme(legend.title = element_blank()) +
  ggtitle("F11 Gene Expression")
F11

plot <- ggarrange(F0, F2, F4, F11, nrow = 2, ncol=2)
annotate_figure(plot, top = text_grob("Generational Gene Expression",
                                      color = "black", face = "bold", size = 20))

###############################################################
## Check on the DE results from the DESeq command way above
resultsNames(dds)
resAM_OWA <- results(dds, name="treatment_OWA_vs_AM", alpha=0.05)
resAM_OWA <- resAM_OWA[order(resAM_OWA$padj),]
head(resAM_OWA)  
summary(resAM_OWA)

resF11_F0 <- results(dds, name="generation_F11_vs_F0", alpha=0.05)
resF11_F0 <- resF11_F0[order(resF11_F0$padj),]
head(resF11_F0)  
summary(resF11_F0)

resF2_F0 <- results(dds, name="generation_F2_vs_F0", alpha=0.05)
resF2_F0 <- resF2_F0[order(resF2_F0$padj),]
head(resF2_F0)  
summary(resF2_F0)

resF4_F0 <- results(dds, name="generation_F4_vs_F0", alpha=0.05)
resF4_F0 <- resF4_F0[order(resF4_F0$padj),]
head(resF4_F0)  
summary(resF4_F0)

### Plot Individual genes ### 
# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN22748_c0_g1::TRINITY_DN22748_c0_g1_i4::g.47585::m.47585", intgroup = (c("treatment","generation")), returnData=TRUE)
d
p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

#################### 
# MODEL NUMBER 2 - subset to focus on effect of generation
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
dim(dds)

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to F0 generation
dds_sub_F0 <- subset(dds, select = generation == 'F0')
dim(dds_sub_F0)

# Perform DESeq2 analysis on the subset
dds_sub_F0 <- DESeq(dds_sub_F0)
resultsNames(dds_sub_F0)
res_F0_OWAvAM <- results(dds_sub_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
head(res_F0_OWAvAM)
summary(res_F0_OWAvAM)

# F2 has no OWA samples

# Subset the DESeqDataSet to F4 generation
dds_sub_F4 <- subset(dds, select = generation == 'F4')
dim(dds_sub_F4)

# Perform DESeq2 analysis on the subset
dds_sub_F4 <- DESeq(dds_sub_F4)
resultsNames(dds_sub_F4)
res_F4_OWAvAM <- results(dds_sub_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_F4_OWAvAM <- res_F4_OWAvAM[order(res_F4_OWAvAM$padj),]
head(res_F4_OWAvAM)
summary(res_F4_OWAvAM)

# Subset the DESeqDataSet to F11 generation
dds$treatment <- factor(dds$treatment, levels = c("OWA","AM"))
dds_sub_F11 <- subset(dds, select = generation == 'F11')
dim(dds_sub_F11)

# Perform DESeq2 analysis on the subset
dds_sub_F11 <- DESeq(dds_sub_F11)
resultsNames(dds_sub_F11)
res_F11_OWAvAM <- results(dds_sub_F11, name="treatment_AM_vs_OWA", alpha=0.05)
res_F11_OWAvAM <- res_F11_OWAvAM[order(res_F11_OWAvAM$padj),]
head(res_F11_OWAvAM)
summary(res_F11_OWAvAM)

########################################
# Plot Individual genes
# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds_sub, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment","generation")), returnData=TRUE)
d
p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

# We can make an MA plot
plotMA(res_F0_OWvAM, ylim=c(-4,4))

########################################

# Heatmap of top 20 genes sorted by pvalue
# By environment
vsd <- vst(dds_sub, blind=FALSE)
topgenes <- head(rownames(res_F0_OWAvAM),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_sub)[,c("generation","treatment")])
pheatmap(mat, annotation_col=df)
pheatmap(mat, annotation_col=df, cluster_cols = F)

########################################

# For OWA vs AM F0
res_F0_OWAvAM <- results(dds_sub_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
head(res_F0_OWAvAM)
summary(res_F0_OWAvAM)
res_F0_OWAvAM <- res_F0_OWAvAM[!is.na(res_F0_OWAvAM$padj),]
degs_F0_OWAvAM <- row.names(res_F0_OWAvAM[res_F0_OWAvAM$padj < 0.05,])
length(degs_F0_OWAvAM)  # 3918

# For OWA vs AM F4
res_F4_OWAvAM <- results(dds_sub_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_F4_OWAvAM <- res_F4_OWAvAM[order(res_F4_OWAvAM$padj),]
head(res_F4_OWAvAM)
summary(res_F4_OWAvAM)
res_F4_OWAvAM <- res_F4_OWAvAM[!is.na(res_F4_OWAvAM$padj),]
degs_F4_OWAvAM <- row.names(res_F4_OWAvAM[res_F4_OWAvAM$padj < 0.05,])
length(degs_F4_OWAvAM)  # 241

# For OWA vs AM F11
res_F11_OWAvAM <- results(dds_sub_F11, name="treatment_AM_vs_OWA", alpha=0.05)
res_F11_OWAvAM <- res_F11_OWAvAM[order(res_F11_OWAvAM$padj),]
head(res_F11_OWAvAM)
summary(res_F11_OWAvAM)
res_F11_OWAvAM <- res_F11_OWAvAM[!is.na(res_F11_OWAvAM$padj),]
degs_F11_OWAvAM <- row.names(res_F11_OWAvAM[res_F11_OWAvAM$padj < 0.05,])
length(degs_F11_OWAvAM)  # 1403

# Intersections
length(intersect(degs_F0_OWAvAM,degs_F4_OWAvAM))  # 57
length(intersect(degs_F0_OWAvAM,degs_F11_OWAvAM))  # 226
length(intersect(degs_F11_OWAvAM,degs_F4_OWAvAM))  # 15

intWA <- intersect(degs_F0_OWAvAM,degs_F4_OWAvAM)
length(intersect(degs_F11_OWAvAM,intWA)) # 2

# Number unique
3918-57-226+2 #3637 F0
241-57-15+2 #171 F4
1403-226-15+2 #1164 F11

57-2 #55 F0 and F4
226-2 #224 F0 and F11
15-2 #13 F11 and F4

# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("F0" = 3637, "F4" = 171, "F11" = 1164, "F0&F4" = 55, "F0&F11" = 224, "F11&F4" = 13, "F0&F4&F11" = 2))
plot(fit1,  lty = 1:3, quantities = TRUE)
plot2 <- plot(fit1, quantities = TRUE, fill = wes_palette("GrandBudapest1"),
     lty = 1:3,
     labels = list(font = 4))

annotate_figure(plot2, top = text_grob("Shared Genes per Generation in OWA vs AM",
                                      color = "black", face = "bold", size = 15))

#cross check
3637+224+55+2 #3918, total F0
171+55+2+13 #241, total F4
1164+224+2+13 #1403, total F11

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_F0_OWAvAM_df <- data.frame(transcriptID = rownames(res_F0_OWAvAM), res_F0_OWAvAM)
res_F4_OWAvAM_df <- data.frame(transcriptID = rownames(res_F4_OWAvAM), res_F4_OWAvAM)
res_F11_OWAvAM_df <- data.frame(transcriptID = rownames(res_F11_OWAvAM), res_F11_OWAvAM)

# Split the "transcriptID" column by double colons and create new columns of the parts
res_F0_OWAvAM_df <- separate(res_F0_OWAvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 
res_F4_OWAvAM_df <- separate(res_F4_OWAvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 
res_F11_OWAvAM_df <- separate(res_F11_OWAvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 

# Create a new column by concatenating "part1" and "part2" with double colons in between
res_F0_OWAvAM_df$transcriptID_trim <- paste(res_F0_OWAvAM_df$part1, res_F0_OWAvAM_df$part2, sep = "::")
res_F4_OWAvAM_df$transcriptID_trim <- paste(res_F4_OWAvAM_df$part1, res_F4_OWAvAM_df$part2, sep = "::")
res_F11_OWAvAM_df$transcriptID_trim <- paste(res_F11_OWAvAM_df$part1, res_F11_OWAvAM_df$part2, sep = "::")

# Select the two columns we want to save for the GOMWU analysis
selected_columns_OWA_F0 <- res_F0_OWAvAM_df[c("transcriptID_trim", "log2FoldChange")]
selected_columns_OWA_F4 <- res_F4_OWAvAM_df[c("transcriptID_trim", "log2FoldChange")]
selected_columns_OWA_F11 <- res_F11_OWAvAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file
write.csv(selected_columns_OWA_F0, file = "res_F0_OWAvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(selected_columns_OWA_F4, file = "res_F4_OWAvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
write.csv(selected_columns_OWA_F11, file = "res_F11_OWAvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
