
# Start with clean environment
rm(list=ls())
gc()

# Read in pFBA fluxes
fluxes <- read.delim('~/Desktop/repos/Cdiff_modeling/data/toy_pfba_fluxes.tsv', header=TRUE, sep='\t')

# Format flux tables
fluxes$reaction <- gsub('_', ' ', fluxes$reaction)
fluxes$flux <- as.numeric(fluxes$flux)
fluxes <- subset(fluxes, type != 'exchange')
a <- subset(fluxes, panel == 'a')
b <- subset(fluxes, panel == 'b')
c <- subset(fluxes, panel == 'c')
rm(fluxes)

# Generate figures
png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_1A.png', units='in', width=2.5, height=6, res=300)
par(mar=c(3,7,1,1.5), las=1, xpd=FALSE, mgp=c(1.5,0.5,0), xaxs='r', lwd=2)
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(-50,1050), ylim=c(0,4.5), pch=15, xlab='pFBA Solution Fluxes', 
     ylab='', cex.lab=0.8)
axis(1, at=c(0,500,1000), cex.axis=0.7, lwd=2)
box(lwd=2)
axis(2, at=seq(0.25,4.25,0.5), labels=rev(a$reaction), tck=0, cex.axis=1.1)
#barplot(rev(a$flux), ylim=c(0,4.5), width=0.41,
#        col='gray40', xaxt='n', yaxt='n', horiz=TRUE, add=TRUE)
dev.off()

png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_1B.png', units='in', width=3, height=3, res=300)
par(mar=c(3,7,1,1.5), las=1, xpd=FALSE, mgp=c(1.5,0.5,0), xaxs='r', lwd=2)
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(-50,1050), ylim=c(0,2.5), pch=15, xlab='pFBA Solution Fluxes', 
     ylab='', cex.lab=1.1)
axis(1, at=c(0,500,1000), cex.axis=0.8, lwd=2)
box(lwd=2)
axis(2, at=seq(0.25,2.25,0.5), labels=rev(b$reaction), tck=0, cex.axis=1.1)
barplot(rev(b$flux), ylim=c(0,2.5), width=0.4, 
        col='#8e7cc3', xaxt='n', yaxt='n', horiz=TRUE, add=TRUE)
dev.off()

png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_1C.png', units='in', width=3, height=4, res=300)
par(mar=c(3,7,1,1.5), las=1, xpd=FALSE, mgp=c(1.5,0.5,0), xaxs='r', lwd=2)
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(-50,1050), ylim=c(0,3.5), pch=15, xlab='pFBA Solution Fluxes', 
     ylab='', cex.lab=1.1)
axis(1, at=c(0,500,1000), cex.axis=0.8, lwd=2)
box(lwd=2)
axis(2, at=seq(0.25,3.25,0.5), labels=rev(c$reaction), tck=0, cex.axis=1.1)
barplot(rev(c$flux), ylim=c(0,3.5), width=0.4,
        col='#e69138', xaxt='n', yaxt='n', horiz=TRUE, add=TRUE)
dev.off()



# Clean up
rm(list=ls())
gc()
