# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.
# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 
# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.
# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.
# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.
# Stretch the plot manually to match tree to text
# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("topGO")
# BiocManager::install("geneLenDataBase")
# install.packages("ape")
library(topGO)
library(geneLenDataBase)
library(quarts)
library(ape)

setwd("C:/Users/bjork/Documents/GitHub/Ecological_Genomics_23/PopGenomics/results/GOMWU")

# Can run for each of 3 inputs, and can do for each GO division BP, MF, CC
# Replaced values with respective files for F0, F4, and F11
input="res_F0_OWAvAM_LFC.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="trinotate_annotation_GOblastx_onlyanns_onlyGOs_justmergeGeneIDtab34.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
	largest=0.05,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
	smallest=10,   # a GO category should contain at least this many genes to be considered
	clusterCutHeight=0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# 585 GO terms at 10% FDR - for F0
# 524 GO terms at 10% FDR - for F4
# 882 GO terms at 10% FDR - for F11

# ----------- Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
 	# absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
  absValue=1, # un-remark this if you are using log2-fold changes
  level1=0.0001, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
 	level2=0.00005, # FDR cutoff to print in regular (not italic) font.
 	level3=0.00001, # FDR cutoff to print in large bold font.
 	txtsize=1.0,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
 	treeHeight=1.0, # height of the hierarchical clustering tree
  colors=c(wes_palette("GrandBudapest1")) # these are default colors, un-remar and change if needed
 )

png(filename = "dendroGOMWU_BP_lesssig_F11_OAvAM_LFC.png", width = 800, height = 1400, res = 200)
dev.off()
# Text representation of results, with actual adjusted p-values
results[[1]]

# ------- extracting representative GOs
# Chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
plot(results[[2]],cex=0.6)
abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
	rn=names(ct)[ct==ci]
	obs=grep("obsolete",rn)
	if(length(obs)>0) { rn=rn[-obs] }
	if (length(rn)==0) {next}
	rr=results[[1]][rn,]
	bestrr=rr[which(rr$pval==min(rr$pval)),]
	best=1
	if(nrow(bestrr)>1) {
		nns=sub(" .+","",row.names(bestrr))
		fr=c()
		for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
		best=which(fr==max(fr))
	}
	if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs

# F0
# delta.rank         pval       level nseqs                  term          name        p.adj
# 1004       2809 1.218108e-32    2   147 GO:0030239;GO:0010927 cellular component assembly involved in morphogenesis 3.443591e-29

# F4
# delta.rank         pval       level nseqs       term                     name        p.adj F4
# 346       1447 6.545853e-26     4   444 GO:0006811 monoatomic ion transport 1.855749e-22

# F11
# delta.rank         pval       level nseqs       term                     name        p.adj
# 347       2117 1.003461e-54     2   432 GO:0006811 monoatomic ion transport 2.791629e-51