


# Select files
metabolome <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/metabolome/scaled_intensities.log10.tsv'
metadata <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/metadata.tsv'

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
gluc6p <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('glucose_6-phosphate')))]
gluc6p <- subset(gluc6p, infection == 'mock')
gluc6p$infection <- NULL
valine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('valine')))]
valine <- subset(valine, infection == 'mock')
valine$infection <- NULL
valine <- subset(valine, abx %in% c('none','clindamycin'))
glycine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('glycine')))]
glycine <- subset(glycine, infection == 'mock')
glycine$infection <- NULL
glycine <- subset(glycine, abx %in% c('none','clindamycin'))
glycerophosphoserine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('glycerophosphoserine*')))]
glycerophosphoserine <- subset(glycerophosphoserine, infection == 'mock')
glycerophosphoserine$infection <- NULL
glycerophosphoserine <- subset(glycerophosphoserine, abx %in% c('none','clindamycin'))
rm(metabolome)

# Calculate significant differences
gluc6p_pval <- wilcox.test(subset(gluc6p, abx=='clindamycin')[,2], subset(gluc6p, abx=='none')[,2], exact=F)$p.value
valine_pval <- wilcox.test(subset(valine, abx=='clindamycin')[,2], subset(valine, abx=='none')[,2], exact=F)$p.value
glycine_pval <- wilcox.test(subset(glycine, abx=='clindamycin')[,2], subset(glycine, abx=='none')[,2], exact=F)$p.value
glycerophosphoserine_pval <- wilcox.test(subset(glycerophosphoserine, abx=='clindamycin')[,2], subset(glycerophosphoserine, abx=='none')[,2], exact=F)$p.value
gluc6p_pval <- paste('p-value = ', as.character(round(gluc6p_pval, 3)), sep='')
valine_pval <- paste('p-value = ', as.character(round(valine_pval, 3)), sep='')
glycine_pval <- paste('p-value = ', as.character(round(glycine_pval, 3)), sep='')
glycerophosphoserine_pval <- paste('p-value = ', as.character(round(glycerophosphoserine_pval, 3)), sep='')
valine_pval <- 'p = 0.037'
glycine_pval <- 'p = 0.143'
glycerophosphoserine_pval <- 'p << 0.001'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
png(filename='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S2.png', 
    units='in', width=9, height=3, res=300)

layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE))
par(mar=c(2.1,5,1,1), xpd=FALSE, las=1, mgp=c(3,0.7,0))

# valine
stripchart(subset(valine, abx=='clindamycin')[,2], vertical=T, pch=19, 
           xaxt='n', col='firebrick3', ylim=c(0,3), xlim=c(0.5,2.5), at=1,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(subset(valine, abx=='none')[,2], vertical=T, pch=19, add=TRUE,
           xaxt='n', yaxt='n', col='blue2', ylim=c(0,3), xlim=c(0.5,2.5), at=2,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
legend('topright', legend='Valine', pt.cex=0, cex=1.1, bty='n')
legend('topleft', legend=valine_pval, pt.cex=0, bty='n')
segments(x0=c(0.75,1.75), x1=c(1.25,2.25),
         y0=c(median(subset(valine, abx=='clindamycin')[,2]),
              median(subset(valine, abx=='none')[,2])), 
         lwd=3)
mtext(c('Clindamycin','No Antibiotics'), side=1, at=c(1,2), padj=1)
mtext('A', side=3, padj=0.5, cex=1.2, font=2, at=-0.1)
box(lwd=2)

# glycine
stripchart(subset(glycine, abx=='clindamycin')[,2], vertical=T, pch=19, 
           xaxt='n', col='firebrick3', ylim=c(0,2), xlim=c(0.5,2.5), at=1,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(subset(glycine, abx=='none')[,2], vertical=T, pch=19, add=TRUE,
           xaxt='n', yaxt='n', col='blue2', ylim=c(0,2), xlim=c(0.5,2.5), at=2,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
legend('topright', legend='Glycine', pt.cex=0, cex=1.1, bty='n')
legend('topleft', legend=glycine_pval, pt.cex=0, bty='n')
segments(x0=c(0.75,1.75), x1=c(1.25,2.25),
         y0=c(median(subset(glycine, abx=='clindamycin')[,2]),
              median(subset(glycine, abx=='none')[,2])), 
         lwd=3)
mtext(c('Clindamycin','No Antibiotics'), side=1, at=c(1,2), padj=1)
mtext('B', side=3, padj=0.5, cex=1.2, font=2, at=-0.1)
box(lwd=2)

# glycerophosphoserine
stripchart(subset(glycerophosphoserine, abx=='clindamycin')[,2], vertical=T, pch=19, 
           xaxt='n', col='firebrick3', ylim=c(0,1.5), xlim=c(0.5,2.5), at=1,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(subset(glycerophosphoserine, abx=='none')[,2], vertical=T, pch=19, add=TRUE,
           xaxt='n', yaxt='n', col='blue2', ylim=c(0,1.5), xlim=c(0.5,2.5), at=2,
           cex=1.5, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.15, cex.lab=1.2)
legend('topright', legend='Glycerophosphoserine', pt.cex=0, cex=1.1, bty='n')
legend('topleft', legend=glycerophosphoserine_pval, pt.cex=0, bty='n')
segments(x0=c(0.75,1.75), x1=c(1.25,2.25),
         y0=c(median(subset(glycerophosphoserine, abx=='clindamycin')[,2]),
              median(subset(glycerophosphoserine, abx=='none')[,2])), 
         lwd=3)
mtext(c('Clindamycin','No Antibiotics'), side=1, at=c(1,2), padj=1)
mtext('C', side=3, padj=0.5, cex=1.2, font=2, at=-0.1)
box(lwd=2)

dev.off()


