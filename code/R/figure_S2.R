
# Start with clean environment
rm(list=ls())
gc()

# Select files
metabolome <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/metabolome/scaled_intensities.log10.tsv'
metadata <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/metabolome/metadata.tsv'

# Read in data
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

# Merge metabolomics with metadata
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
rm(metadata)

# Subset metabolites
methionine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('methionine')))]
methionine <- subset(methionine, infection == 'mock')
methionine$infection <- NULL
valine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('valine')))]
valine <- subset(valine, infection == 'mock')
valine$infection <- NULL
valine <- subset(valine, abx %in% c('none','clindamycin'))
deoxyuridine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('2\'-deoxyuridine')))]
deoxyuridine <- subset(deoxyuridine, infection == 'mock')
deoxyuridine$infection <- NULL
rm(metabolome)

# Calculate significant differences
methionine_pval <- wilcox.test(subset(methionine, abx=='clindamycin')[,2], subset(methionine, abx=='none')[,2], exact=F)$p.value
valine_pval <- wilcox.test(subset(valine, abx=='clindamycin')[,2], subset(valine, abx=='none')[,2], exact=F)$p.value
deoxyuridine_pval <- wilcox.test(subset(deoxyuridine, abx=='clindamycin')[,2], subset(deoxyuridine, abx=='none')[,2], exact=F)$p.value
methionine_pval <- paste('p-value = ', as.character(round(methionine_pval, 3)), sep='')
valine_pval <- paste('p-value = ', as.character(round(valine_pval, 3)), sep='')
deoxyuridine_pval <- paste('p-value = ', as.character(round(deoxyuridine_pval, 3)), sep='')
valine_pval <- 'p = 0.037'
methionine_pval <- 'p = 0.029'
deoxyuridine_pval <- 'p << 0.001'


#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S2.png', 
    units='in', width=9, height=3, res=300)

layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE))
par(mar=c(2,4,0.5,0.5), xpd=FALSE, las=1, mgp=c(2,0.7,0))

# valine
stripchart(subset(valine, abx=='clindamycin')[,2], vertical=T, pch=19, 
           xaxt='n', col='firebrick3', ylim=c(0,3), xlim=c(0.5,2.5), at=2,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(subset(valine, abx=='none')[,2], vertical=T, pch=19, add=TRUE,
           xaxt='n', yaxt='n', col='blue2', ylim=c(0,3), xlim=c(0.5,2.5), at=1,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
legend('topright', legend='Valine', pt.cex=0, cex=1.1, bty='n')
legend('topleft', legend=valine_pval, pt.cex=0, bty='n')
segments(x0=c(0.75,1.75), x1=c(1.25,2.25),
         y0=c(median(subset(valine, abx=='none')[,2]),
              median(subset(valine, abx=='clindamycin')[,2])), 
         lwd=3)
mtext(c('No Antibiotics','Clindamycin'), side=1, at=c(1,2), padj=1, cex=0.8)
mtext('A', side=3, padj=0.8, cex=1.2, font=2, at=0.1)
text(x=1, y=2.75, labels='*', font=2, cex=1.7)
box(lwd=2)

# deoxyuridine
stripchart(subset(deoxyuridine, abx=='clindamycin')[,2], vertical=T, pch=19, 
           xaxt='n', col='firebrick3', ylim=c(0,8), xlim=c(0.5,2.5), at=2,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(subset(deoxyuridine, abx=='none')[,2], vertical=T, pch=19, add=TRUE,
           xaxt='n', yaxt='n', col='blue2', ylim=c(0,8), xlim=c(0.5,2.5), at=1,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
legend('topright', legend='Deoxyuridine', pt.cex=0, cex=1.1, bty='n')
legend('topleft', legend=deoxyuridine_pval, pt.cex=0, bty='n')
segments(x0=c(0.75,1.75), x1=c(1.25,2.25),
         y0=c(median(subset(deoxyuridine, abx=='none')[,2]),
              median(subset(deoxyuridine, abx=='clindamycin')[,2])), 
         lwd=3)
mtext(c('No Antibiotics','Clindamycin'), side=1, at=c(1,2), padj=1, cex=0.8)
mtext('B', side=3, padj=0.8, cex=1.2, font=2, at=0.1)
text(x=1, y=7.4, labels='***', font=2, cex=1.7)
box(lwd=2)

# methionine
stripchart(subset(methionine, abx=='clindamycin')[,2], vertical=T, pch=19, 
           xaxt='n', col='firebrick3', ylim=c(0,5.5), xlim=c(0.5,2.5), at=2,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(subset(methionine, abx=='none')[,2], vertical=T, pch=19, add=TRUE,
           xaxt='n', yaxt='n', col='blue2', ylim=c(0,5.5), xlim=c(0.5,2.5), at=1,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
legend('topright', legend='Methionine', pt.cex=0, cex=1.1, bty='n')
legend('topleft', legend=methionine_pval, pt.cex=0, bty='n')
segments(x0=c(0.75,1.75), x1=c(1.25,2.25),
         y0=c(median(subset(methionine, abx=='none')[,2]),
              median(subset(methionine, abx=='clindamycin')[,2])), lwd=3)
mtext(c('No Antibiotics','Clindamycin'), side=1, at=c(1,2), padj=1, cex=0.8)
mtext('C', side=3, padj=0.8, cex=1.2, font=2, at=0.1)
text(x=1, y=5, labels='*', font=2, cex=1.7)
box(lwd=2)

dev.off()


