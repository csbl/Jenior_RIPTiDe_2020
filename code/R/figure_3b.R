
library(gplots)

# Load in data
essential <- as.data.frame(t(read.delim(file='~/Desktop/repos/Jenior_RIPTiDe_2019/data/essentiality_test2.tsv', sep='\t', header=TRUE, row.names=1)))

# Assess differential groups of essentiality
test <- as.data.frame(t(essential))
pfba_only <- subset(test, pfba == 2)
pfba_only <- subset(pfba_only, m9_aerobic == 0)
pfba_only <- subset(pfba_only, m9_anaerobic == 0)
pfba_only <- subset(pfba_only, lb_aerobic == 0)
pfba_only <- rownames(pfba_only)
lb_aerobic_only <- subset(test, lb_aerobic == 2)
lb_aerobic_only <- subset(lb_aerobic_only, m9_aerobic == 0)
lb_aerobic_only <- subset(lb_aerobic_only, m9_anaerobic == 0)
lb_aerobic_only <- subset(lb_aerobic_only, pfba == 0)
lb_aerobic_only <- rownames(lb_aerobic_only)
m9_aerobic_only <- subset(test, m9_aerobic == 2)
m9_aerobic_only <- subset(m9_aerobic_only, lb_aerobic == 0)
m9_aerobic_only <- subset(m9_aerobic_only, m9_anaerobic == 0)
m9_aerobic_only <- subset(m9_aerobic_only, pfba == 0)
m9_aerobic_only <- rownames(m9_aerobic_only)
m9_anaerobic_only <- subset(test, m9_anaerobic == 2)
m9_anaerobic_only <- subset(m9_anaerobic_only, lb_aerobic == 0)
m9_anaerobic_only <- subset(m9_anaerobic_only, m9_aerobic == 0)
m9_anaerobic_only <- subset(m9_anaerobic_only, pfba == 0)
m9_anaerobic_only <- rownames(m9_anaerobic_only)


test <- subset(test, m9_anaerobic == 0)
test <- subset(test, lb_aerobic == 2)
test <- subset(test, m9_aerobic == 2)
test <- subset(test, pfba == 2)
test <- rownames(test)


# Print unique gene analysis results
print(pfba_only)
print(lb_aerobic_only)
print(m9_aerobic_only)
print(m9_anaerobic_only)

# Format data for plotting
genres <- rownames(essential)
essential <- as.matrix(sapply(essential, as.numeric))
rownames(essential) <- genres
rm(genres)

# Generate figure
png(filename='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/essentiality2.png', 
    units='in', width=6, height=4, res=100)
heatmap.2(essential, col=c('black','white','forestgreen'), dendrogram='none', density.info='none', 
          trace='none', key=FALSE, margins=c(1,1), Rowv=FALSE, sepwidth=c(0.05,0.01), 
          sepcolor='black', colsep=1:ncol(essential), rowsep=1:nrow(essential), labCol=FALSE)
dev.off()

# Clean up
rm(list=ls())
gc()
