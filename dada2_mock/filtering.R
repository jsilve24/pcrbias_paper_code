# Applying dada2 pipeline to bioreactor time-series
## Following tutorial http://benjjneb.github.io/dada2/tutorial.html
## and here http://benjjneb.github.io/dada2/bigdata.html

setwd('/data/davidlab/users/jds/pcrbias/results/dada2_second_experiment/')

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")

path <- "/data/davidlab/users/jds/pcrbias/results/dada2_second_experiment/"

# Filtering and Trimming --------------------------------------------------

# Forward and Reverse Filenames
fnFs.s1 <-list.files(paste(path, "3_demultiplex/s1r1_split_samples/", sep=""))
fnRs.s1 <-list.files(paste(path, "3_demultiplex/s1r2_split_samples/", sep=""))

# Sort to ensure filenames are in the same order
fnFs.s1 <- sort(fnFs.s1)
fnRs.s1 <- sort(fnRs.s1)

sample.names.1 <- sapply(strsplit(fnFs.s1,".fastq", fixed=TRUE), `[`, 1)

# Remove extraneous samples ## EDIT THIS AS NEEDED
# to.keep <- grepl("^V", sample.names.1)
# sample.names.1 <- sample.names.1[to.keep]
# print(sample.names.1)
# fnFs.s1 <- fnFs.s1[to.keep] 
# fnRs.s1 <- fnRs.s1[to.keep]

# Fully Specify the path for the fnFs and fnRs
fnFs.s1 <- file.path(path, '3_demultiplex/s1r1_split_samples/', fnFs.s1)
fnRs.s1 <- file.path(path, '3_demultiplex/s1r2_split_samples/', fnRs.s1)

# Examine qulaity profiles of the forward and reverse reads
p <- plotQualityProfile(fnFs.s1[[1]]) 
ggsave('Forward_quality_profile_s1.png', plot=p)
p <- plotQualityProfile(fnRs.s1[[1]]) 
ggsave('Reverse_quality_profile_s1.png', plot=p)

# Perform filtering and trimming
filtpath <- file.path(path, "4_filter")
dir.create(filtpath)

# For the first sequencing run
dir.create(file.path(filtpath,'s1'))
filtFs.s1 <- file.path(filtpath, 's1', paste0(sample.names.1,"_F_filt.fastq.gz"))
filtRs.s1 <- file.path(filtpath, 's1', paste0(sample.names.1,"_R_filt.fastq.gz"))
for (i in seq_along(fnFs.s1)){
    fastqPairedFilter(c(fnFs.s1[i], fnRs.s1[i]), c(filtFs.s1[i], filtRs.s1[i]),
                      trimLeft=c(10,0), truncLen=c(150, 140),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE, verbose=TRUE,
                      rm.phix=TRUE)
}
