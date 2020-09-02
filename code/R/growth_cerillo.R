

# Start with clean environment
rm(list=ls())
gc()

# Read in data
growth <- read.delim('~/Desktop/bhi_test.tsv', sep='\t', header=TRUE)

# Subset blank wells and find medians over time
library(matrixStats)
blank_wells <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12')
blank_series <- as.matrix(growth[, blank_wells])
blank_series <- rowMedians(blank_series)
growth <- growth[, !(colnames(growth) %in% blank_wells)]
rm(blank_wells)

# Normalize by blanks
for (x in colnames(growth)) {growth[,x] <- growth[,x] - blank_series}
growth[growth < 0] <- 0
rm(x, blank_series)

# Generate figure
par(mar=c(3,3,1,1), xpd=FALSE, mgp=c(2,1,0), las=1, lwd=2)
plot(0, xlim=c(1,nrow(growth)), ylim=c(0,1), type='n', xlab='Time point', ylab='OD589')
for (x in colnames(growth)) {lines(growth[,1], type='l', col='gray')}
lines(rowMedians(as.matrix(growth)), type='l', col='firebrick', lwd=3)
rm(x)




