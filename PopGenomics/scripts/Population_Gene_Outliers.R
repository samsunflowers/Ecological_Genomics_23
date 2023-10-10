###################################
#  Selection scans for red spruce #
###################################
#install.packages("RcppCNPy")

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("C:/Users/bjork/Documents/GitHub/Ecological_Genomics_23/PopGenomics/results/Homework1.2")

list.files()

### read in selection statistics (these are chi^2 distributed)

s<-npyLoad("allRS_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))
names(pval) = c("p_PC1","p_PC2")
pval

## read positions
p <- read.table("allRS_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
dim(p)

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

# get all the outliers with p-values below some cutoff
cutoff=1e-3   

outliers_PC1 <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]
outliers_PC2 <- p_filtered[which(pval$p_PC2<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1)[1]
dim(outliers_PC2)[1]

# write them out to a file
write.table(outliers_PC1,
            "allRS_poly_outliers_PC1.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)

COV <- as.matrix(read.table("allRS_poly.cov"))

PCA <- eigen(COV)

data=as.data.frame(PCA$vectors)
data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep
head(data)

write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)
