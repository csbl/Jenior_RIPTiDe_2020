#install.packages('gplots')
library(gplots)
#install.packages('RColorBrewer')
library(RColorBrewer)


# UVA colors
hot <- '#f47320'
cold <- '#0c275d'
uva_palette <- colorRampPalette(c(cold, hot))(n=100)

# Read in data
aerobic <- read.delim('/home/mjenior/Desktop/repos/Cdiff_modeling/data/aerobic_concordance.tsv', header=TRUE, sep='\t', row.names=1)
aerobic <- data.matrix(aerobic)
anaerobic <- read.delim('/home/mjenior/Desktop/repos/Cdiff_modeling/data/anaerobic_concordance.tsv', header=TRUE, sep='\t', row.names=1)
anaerobic <- data.matrix(anaerobic)

# Generate heatmap
pdf(file='~/Desktop/aerobic_heatmap.pdf', width=5, height=8)
heatmap.2(aerobic,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=uva_palette,       # use on color palette defined earlier
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")   
dev.off()

pdf(file='~/Desktop/anaerobic_heatmap.pdf', width=5, height=8)
heatmap.2(anaerobic,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=uva_palette,       # use on color palette defined earlier
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")   
dev.off()


