args <- commandArgs(trailingOnly = TRUE)
cov_in    <- args[1]
out_table <- args[2]
out_plot  <- args[3]

cov <- read.table(cov_in, sep="\t", header=F, stringsAsFactors = FALSE)
colnames(cov) <- c('sample', 'chrM', 'genome')
cov$copies <- with(cov, chrM/genome)

# barplot of chrM coverage normalized for genome coverage
pdf(out_plot)
bp <- barplot(cov$copies, ylab='Copies of ChrM', xlab='Samples', las=2, ylim=c(1, max(cov$copies)*1.10))
offset <- max(cov$copies * 0.01)
text(x = bp[,1], y=cov$copies+offset, labels=round(cov$copies,0))
text(x = bp[,1], y=cov$copies-offset, labels=round(cov$genome,3))
dev.off()

# table corresponding to the chrM coverage
write.table(cov, file=out_table, sep="\t", col.names = TRUE, row.names = FALSE)
