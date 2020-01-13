
# Start with clean environment
rm(list=ls())
gc()

# Load in data
transcription <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/transcript/clinda_k12.mapped.norm.tsv', header=TRUE, row.names=1)
weights <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/clinda_k12.weights.tsv')
colnames(weights) <- c('reaction', 'weight')

# Format data
transcription <- transcription$normDepth
pruning_weights <- transcription / max(transcription)
sampling_weights <- unique(sort(weights$weight))

# Calculate density curves
transcript_density <- density(transcription)
pruning_weight_density <- density(pruning_weights)

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_5A.png', units='in', width=4, height=6, res=300)
par(mar=c(4,4,1,4), las=1, mgp=c(2.6,1,0), xaxs='i', yaxs='i', lwd=2, xpd=FALSE)
hist(transcription, main='', xlim=c(0,500), ylim=c(0,800), breaks=200, col='firebrick',
     xlab=as.expression(bquote(paste(italic('in vivo'),' Transcript Density', ylab='', sep=''))), 
     ylab='', cex.lab=1.2, lwd=2, xaxt='n')
legend('top', legend=c('Reactions: 452','Metabolites: 453','Biomass Flux: 17.74%↓'), cex=0.9, pt.cex=0, bty='n')
axis(side=1, at=c(0,250,500), labels=c(0,2500,5000), lwd=2)
axis(side=4, at=seq(0,800,160), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), lwd=2)
par(new=TRUE)
plot(rev(sampling_weights), type='l', col='dodgerblue', lwd=5, xlim=c(-120,790), xaxt='n', yaxt='n', xlab='', ylab='')
par(new=TRUE)
plot(sampling_weights, type='l', col='darkorchid4', lwd=5, xlim=c(1,790), xaxt='n', yaxt='n', xlab='', ylab='')
par(xpd=TRUE)
text(x=-200, y=0.5, 'Gene Frequency', srt=90, cex=1.2, col='firebrick', font=2)
text(x=1000, y=0.5, 'Reaction Weights (               &                 )', srt=270, cex=1.2, font=2)
text(x=1000, y=0.44, 'Pruning', srt=270, cex=1.2, col='dodgerblue', font=2)
text(x=1000, y=0.226, 'Sampling', srt=270, cex=1.2, col='darkorchid4', font=2)
par(xpd=FALSE)
box(lwd=2)
dev.off()




par(mar=c(5,6,1,3), las=1, mgp=c(3.5,1,0), xaxs='i', yaxs='i', xpd=FALSE)
plot(pruning_weight_density, xlim=c(-12,12), ylim=c(0,0.14), main='', xaxt='n', cex.lab=1.7, cex.axis=1.1,
     xlab='Simulated Metabolite Score', ylab='Score Density') 
axis(side=1, at=seq(-12,12,4), labels=seq(-12,12,4), cex.axis=1.2, tck=-0.03)
axis(side=1, at=c(-12:12), tck=-0.015, labels=FALSE)
polygon(pruning_weight_density, col='white', border='red', lwd=2) 
abline(v=score_median, lty=2, lwd=2, col='black') # Median
abline(v=c(lower_95,upper_95), lty=2, lwd=2, col='red') # 0.95 Confidence Interval
abline(v=7, lwd=2, col='blue') # Measured Score
legend('topleft', legend=c('Simulated Median','Simulated 95% Confidence Interval','Measured Metabolite Score'), 
       col=c('black','red','blue'), lty=c(2,2,1), lwd=2, bg='white', cex=1.1)



