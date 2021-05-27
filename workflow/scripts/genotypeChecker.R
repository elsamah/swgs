#### Author : A. Danesh
#### First created : Feb 2014
##### This code find the jaccard indices between samples from a vcf input file. It is made for Sample check #####
### The vcf file : The vcf file is created by merging all samples' vcfs from unifiedgenotyper using vcftools## the merged vcf has all the genotype data after the 10th column###
###########
#arguments##
###########

args <- commandArgs(trailingOnly = TRUE)
in_vcf           <- args[1]
sample_ids       <- args[2]
out_jacc_table   <- args[3]
out_n_table      <- args[4]
out_heatmap      <- args[5]


###########
#Functions#
###########
#### Genotype reader
getGT <- function(x, gt_idx=1){
  split_gt <- strsplit(as.vector(x), split = ":")
  sapply(split_gt, function(i) i[[gt_idx]])
}

cleanGenotype <- function(b) {
  # GT genotype, encoded as alleles values separated by either of ”/” or “|”,
  # e.g. The allele values are 0 for the reference allele (what is in the
  # reference sequence), 1 for the first allele listed in ALT, 2 for the
  # second allele list in ALT and so on. For diploid calls examples could be
  # 0/1 or 1|0 etc. For haploid calls, e.g. on Y, male X, mitochondrion,
  # only one allele value should be given. All samples must have GT call
  # information; if a call cannot be made for a sample at a given locus, ”.”
  # must be specified for each missing allele in the GT field (for example
  # ./. for a diploid). The meanings of the separators are:
  stopifnot(any(class(b)=='matrix'))
  
  b <- gsub("\\|", "/", b)                 # convert phased to unphased
  b[b == '.'] <- "0/0"                     # convert No-call to HOMREF
  
  # Deeper sanity checks
  splgt <- sapply(strsplit(as.vector(b), "/"), function(i) i[1:2])
  class(splgt) <- 'integer'
  
  # CHECK_1: Flips 1/0 to 0/1
  hetidx <- which(splgt[2,] != splgt[1,])
  rev_hetidx <- splgt[1,hetidx] > splgt[2,hetidx]
  if(any(rev_hetidx)){
    for(i in which(rev_hetidx)){
      splgt[,hetidx[i]] <- rev(splgt[,hetidx[i]])
    }
  }
  
  # CHECK_2: ...
  if(FALSE){
  
  }
  
  genotype <- matrix(apply(splgt, 2, paste, collapse="/"), ncol=ncol(b))
  return(genotype)
}

######
#Main##
#######
### Read the Sample genotype table and assign columns###
Mega_VCF = read.csv(in_vcf, sep = "\t", comment.char = "#")
rem_ref <- TRUE

sample_sub = Mega_VCF[,10:ncol(Mega_VCF)]
colnames(sample_sub) <- strsplit(sample_ids, ",")[[1]]

# Parse haplotype caller into the first column (i.e. 0/0, or 1/1, or 0/1).
sample_geno <- apply(sample_sub, 2, getGT)
sample_count <- cleanGenotype(as.matrix(sample_geno))

## Read all the samples and count the genotypes####
jacc_n <- apply(sample_count, 2, function(i){
  comp_jn <- apply(sample_count, 2, function(j){
    # Checks that each samples SNP matches the other samples SNP
    comp_idx <- c(1:length(i))
    if(rem_ref){
      #print("Removing SNPs which are ref in both samples")
      refs <- i == '0/0' & j=='0/0'
      comp_idx <- comp_idx[-which(refs)]
    }
    m <- sum(i[comp_idx] == j[comp_idx])
    n <- length(comp_idx)
    jacc <- sum(m)/n
    return(c(jacc,n))
  })
  return(list(round(comp_jn[1,],3), comp_jn[2,]))
})
jacc_mat <- sapply(jacc_n, function(i) i[[1]])
n_mat    <- sapply(jacc_n, function(i) i[[2]])

colnames(jacc_mat) <- rownames(jacc_mat) <- strsplit(sample_ids, ",")[[1]]
colnames(n_mat) <- rownames(n_mat) <- strsplit(sample_ids, ",")[[1]]


write.table(round(jacc_mat,3), file = out_jacc_table,
            row.names = TRUE, col.names = NA , sep = "\t")
write.table(round(n_mat,3), file = out_n_table,
            row.names = TRUE, col.names = NA , sep = "\t")


pdf(file=out_heatmap)
require(ggplot2)
require(reshape2)
cl <- hclust(dist(jacc_mat))

m_jacc   <- melt(jacc_mat[cl$order,cl$order])
m_n      <- melt(n_mat[cl$order,cl$order])
m_jacc_n <- merge(m_jacc, m_n, by=c('Var1', 'Var2'))
colnames(m_jacc_n) <- c('Var1', 'Var2', 'jacc', 'n')

ggplot(m_jacc_n, aes(y=factor(Var1), x=factor(Var2))) +
geom_point(aes(colour=jacc, size=n)) +
scale_color_gradient2(low='blue', mid='gray', high='red',
                      midpoint=0.5) +
theme_bw()
dev.off()
