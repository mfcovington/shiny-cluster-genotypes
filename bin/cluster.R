library(cluster)
library(plyr)

input.file <- 'sample-data.snps.txt'

na.strings <- c('N', 'U')
genotypes <- read.table(file = input.file, sep = '\t', header = TRUE,
                        row.names = 1, na.strings = na.strings)
sample.ids <- colnames(genotypes)

genotypes <- as.data.frame(t(genotypes))
genotypes <- colwise(as.factor)(genotypes)

distances <- daisy(genotypes, metric = 'gower')    # gower for non-numeric data
