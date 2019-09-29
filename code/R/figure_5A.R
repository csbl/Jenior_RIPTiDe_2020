
# Start with clean environment
rm(list=ls())
gc()

# Load in data
transcription <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/transcript/clinda_k12.mapped.norm.tsv', header=TRUE, row.names=1)
weights <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/clinda_k12.weights.tsv')
colnames(weights) <- c('reaction', 'weight')

# Format data
transcription <- transcription$normDepth
weights <- unique(sort(weights$weight))
#weights <- 1.0 - weights

# Calculate density curves
transcript_density <- density(transcription)
weight_density <- density(weights)

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_5A.png', units='in', width=4, height=6, res=300)
par(mar=c(4,4,1,4), las=1, mgp=c(2.6,1,0), xaxs='i', yaxs='i', lwd=2)
hist(transcription, main='', xlim=c(0,500), ylim=c(0,800), breaks=200, col='firebrick',
     xlab=as.expression(bquote(paste(italic('in vivo'),' Transcript Density', ylab='', sep=''))), 
     ylab='', cex.lab=1.2, lwd=2, xaxt='n')
box(lwd=2)
legend('top', legend=c('Reactions: 452','Metabolites: 453','Biomass Flux: 17.74%↓'), cex=0.9, pt.cex=0, bty='n')
axis(side=1, at=c(0,250,500), labels=c(0,2500,5000), lwd=2)
axis(side=4, at=seq(0,800,160), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), lwd=2)
par(new=TRUE)
plot(weights, type='l', col='dodgerblue', lwd=5, xlim=c(1,790), xaxt='n', yaxt='n', xlab='', ylab='')
par(xpd=TRUE)
text(x=-200, y=0.5, 'Gene Frequency', srt=90, cex=1.2, col='firebrick', font=2)
text(x=1000, y=0.5, 'Reaction Weight', srt=270, cex=1.2, col='dodgerblue', font=2)
dev.off()
