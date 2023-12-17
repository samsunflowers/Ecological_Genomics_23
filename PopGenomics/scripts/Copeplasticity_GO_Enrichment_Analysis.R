# import package
library("ggplot2")
library("wesanderson")

setwd("C:/Users/bjork/Documents/GitHub/Ecological_Genomics_23/PopGenomics/results")

allgos <- read.csv("Shared_GO.csv")
allhudsgo <- read.csv("Shared_HudsonicaGO.csv")
top5hudsgo <- read.csv("Top_5_HudsGo.csv")
tonsgo <- read.csv("Tonsa_GO.csv")

colnames(allgos) <- c("GO_Category", "Species", "Generation", "Pval", "Gene_Ratio", "Gene Number", "Gene Total", "-Log10(P-value)")
colnames(allhudsgo) <- c("GO_Category", "Species", "Generation", "Pval", "Gene_Ratio", "Gene Number", "Gene Total", "-Log10(P-value)")
colnames(top5hudsgo) <- c("GO_Category", "Species", "Generation", "Pval", "Gene_Ratio", "Gene Number", "Gene Total", "-Log10(P-value)")
colnames(tonsgo) <- c("GO_Category", "Species", "Generation", "Pval", "Gene_Ratio", "Gene Number", "Gene Total", "-Log10(P-value)")

tonsago <- (allgos[allgos$Species %in% "tonsa", ])
hudsonicago <- (allgos[allgos$Species %in% "hudsonica", ])

# plot: dot plot
ggplot(data = top5hudsgo, aes(x = Generation, y = GO_Category, 
                        color = `-Log10(P-value)`, size = Gene_Ratio)) + 
  geom_point() +
  scale_color_gradient(low = "#FD6467", high = "#5B1A18") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO Enrichment Analysis A. hudsonica")

ggplot(data = allhudsgo, aes(x = Generation, y = GO_Category, 
                             color = `-Log10(P-value)`, size = Gene_Ratio)) + 
  geom_point() +
  scale_color_gradient(low = "#FD6467", high = "#5B1A18") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO Enrichment Analysis A. hudsonica")

ggplot(data = tonsgo, aes(x = Generation, y = GO_Category, 
                           color = `-Log10(P-value)`, size = Gene_Ratio)) + 
  geom_point() +
  scale_color_gradient(low = "#FD6467", high = "#5B1A18") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO Enrichment Analysis A. tonsa")

ggplot(data = hudsonicago, aes(x = Generation, y = GO_Category, 
                           color = `-Log10(P-value)`, size = Gene_Ratio)) + 
  geom_point() +
  scale_color_gradient(low = "#FD6467", high = "#5B1A18") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO Enrichment Analysis A. hudsonica")

