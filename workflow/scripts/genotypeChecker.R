#### Author : A. Danesh
#### First created : Feb 2014
#### Last update : April 2017 (Rene: Sped up the code and removed fluff)
##### This code find the jaccard indices between samples from a vcf input file. It is made for Sample check #####
### The vcf file : The vcf file is created by merging all samples' vcfs from unifiedgenotyper using vcftools## the merged vcf has all the genotype data after the 10th column###
###########
#arguments##
###########

args <- commandArgs(trailingOnly = TRUE)
in_vcf      <- args[1]
sample_ids  <- args[2]
out_table   <- args[3]
out_heatmap <- args[4]

###########
#Functions#
###########
#### Genotype reader
genotype <- function(b) {
  c = b[[1]][1]
  if (c == "." | c == "0/0"){genotype = "0/0"}
  else {genotype = c}
  return(genotype)
}

######
#Main##
#######
### Read the Sample genotype table and assign columns###
Mega_VCF = read.csv(in_vcf, sep = "\t", comment.char = "#")

sample_sub = Mega_VCF[,10:ncol(Mega_VCF)]
colnames(sample_sub) <- strsplit(sample_ids, ",")[[1]]

# Parse haplotype caller into the first column (i.e. 0/0, or 1/1, or 0/1).  Stores as dataframe
sample_geno <- apply(sample_sub, 2, function(x) strsplit(as.vector(x), split = ":"))
sample_count <- lapply(sample_geno, function(i) sapply(i, genotype))
sample_count.data <- as.data.frame(do.call("cbind", sample_count))


## Read all the samples and count the genotypes####
mat <- apply(sample_count.data, 2, function(i){
  apply(sample_count.data, 2, function(j){
    # Checks that each samples SNP matches the other samples SNP
    sum_match <- sapply(c(1:length(i)), function(x){
      as.character(i[x]) == as.character(j[x])
    })
    sum_match <- table(sum_match)['TRUE']
    jaccard.val <- sum_match/length(i)
    return(jaccard.val)
  })
})
colnames(mat) <- rownames(mat) <- strsplit(sample_ids, ",")[[1]]


write.table(round(mat,3), file = out_table,
            row.names = TRUE, col.names = NA , sep = "\t")

pdf(file=out_heatmap)
heatmap(mat, Colv=NA, Rowv=NA)
dev.off()
