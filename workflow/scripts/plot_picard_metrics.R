suppressPackageStartupMessages(library(optparse))
## TODO: Add in an argument to accept csv sample names (e.g. net-001a,net00-1b)
## and use only those samples to create plots.

####################
#### Parameters ####
option_list <- list(
  make_option(c("-q", "--qcdir"), type="character", default='results/qc',
              help="Directory containining picard metric files [default=%default]"),
  make_option(c("-o", "--output"), type="character", default='picard_qc.pdf',
              help="File name and path to output file [default=%default]"),
  make_option(c("-s", "--samples"), type="character", default='all',
              help="Comma-separated list of samples to plot [default=%default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
qcdir <- opt$qcdir
out_pdf <- opt$output
samples <- opt$samples

##########################
#### Helper Functions ####
readPicardFile <- function(f, skipped_rows=6, nrows=1){
  read.table(f, skip=skipped_rows, nrow=nrows,
             comment.char="", header=TRUE, fill=T,
             check.names=FALSE, stringsAsFactors=FALSE)
}

getMaxY <- function(x, range_y=c(0.1, 0.5, 1, 2, 5, 10, 50, 500, 1000)){
  max_x <- max(as.numeric(unlist(x)))*1.25
  idx   <- which.min(abs(max_x-range_y))
  range_y[idx]
}

#################################
#### List and read all files ####
wgs_files     <- list.files(qcdir, pattern="wgs_metrics")         # CollectWgsMetrics
insert_files  <- list.files(qcdir, pattern="insert_size_metrics") # CollectInsertSizeMetrics
if(samples=='all'){
  samples     <- gsub(".wgs_metrics$", "", wgs_files)
} else {
  samples     <- strsplit(samples, split=",")[[1]]
  idx         <- match(samples, gsub(".wgs_metrics", "", wgs_files))
  wgs_files   <- wgs_files[idx]
  insert_files <- insert_files[idx]
}


wgs_df    <- do.call(rbind, lapply(file.path(qcdir, wgs_files), readPicardFile))
insert_df <- do.call(rbind, lapply(file.path(qcdir, insert_files), readPicardFile))
rownames(wgs_df) <- rownames(insert_df) <- samples

#####################
#### Plotting QC ####
pdf(out_pdf, height=12, width=10)
split.screen(c(3,1))
# Mean Coverage
screen(1)
par(mar=c(1,4.1, 2, 2.1))
cov_col <- c(bait="gray")
barplot(t(as.matrix(wgs_df[,c('MEAN_COVERAGE'), drop=FALSE])),
        ylab="Mean coverage", xaxt='n', las=2, cex.names=0.75,
        ylim=c(0, getMaxY(wgs_df[,c('MEAN_COVERAGE')])), col=cov_col)

screen(2)
par(mar=c(1,4.1, 2, 2.1))
barplot(t(as.matrix(insert_df[,c('MEDIAN_INSERT_SIZE'), drop=FALSE])),
        ylab="Median insert size", xaxt='n', las=2, cex.names=0.75,
        ylim=c(0, getMaxY(insert_df[,c('MEDIAN_INSERT_SIZE')])), col=cov_col)

# Failed/filtered bases:
screen(3)
par(mar=c(6.1,4.1, 1, 2.1))
pct_col <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
pct <- c("PCT_EXC_MAPQ", "PCT_EXC_DUPE", "PCT_EXC_UNPAIRED",
         "PCT_EXC_BASEQ", "PCT_EXC_OVERLAP")
barplot(t(wgs_df[,pct]), col=pct_col, cex.names=0.75,
        ylim=c(0, getMaxY(wgs_df[,pct])), las=2,
        ylab="Fraction of filtered aligned bases")
legend("topright", fill=pct_col, cex=1,
       legend=c("MAPQ < 20", "Duplicates", "Singletons",
       "BASEQ < 20", "Overlapping reads"))
close.screen(all.screens=TRUE)
