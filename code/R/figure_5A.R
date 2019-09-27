
# Start with clean environment
rm(list=ls())
gc()

# Load in data
transcription <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/transcript/clinda_k12.mapped.norm.tsv', header=TRUE, row.names=1)

# Format data
transcription <- transcription$normDepth

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_5A.png', units='in', width=4, height=6, res=300)
par(mar=c(4,4,1.5,1), las=1, mgp=c(2.6,1,0), xaxs='i', yaxs='i', lwd=2)
hist(transcription, main='', xlim=c(0,400), ylim=c(0,800), breaks=200, col='firebrick',
     xlab=as.expression(bquote(paste(italic('In vivo'),' Transcript Density',sep=''))), 
     ylab='Gene Frequency', cex.lab=1.2, lwd=2)
box(lwd=2)
legend('topright', legend=c(as.expression(bquote(italic('in vivo'))),'Reactions: 452','Metabolites: 453','Biomass Flux: 17.74%↓'), 
       cex=1.2, pt.cex=0, bty='n')
text(x=280, y=758, labels='-specific Model', cex=1.2)

dev.off()






