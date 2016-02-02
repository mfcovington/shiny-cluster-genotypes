library(cluster)
library(ggplot2)
library(ggrepel)
library(limma)
library(plyr)

experiment.id <- 'SNP Data'
input.file <- 'sample-data.snps.txt'
delimiter <- '\t'

na.strings <- c('N', 'U')
genotypes <- read.table(file = input.file, sep = delimiter, header = TRUE,
                        row.names = 1, na.strings = na.strings)
sample.ids <- colnames(genotypes)

mds <- plotMDS(data.matrix(genotypes))
mds.df <- as.data.frame(mds$cmdscale.out)
mds.df$idx <- 1:length(sample.ids)

ggplot(mds.df) +
  geom_point(aes(x = V1, y = V2)) +
  geom_label_repel(aes(x = V1, y = V2, label = idx)) +
  ggtitle(expression('Leading ' * Log[2] * ' fold change MDS')) +
  xlab('Dimension 1') +
  ylab('Dimension 2')

genotypes <- as.data.frame(t(genotypes))
genotypes <- colwise(as.factor)(genotypes)

distances <- daisy(genotypes, metric = 'gower')    # gower for non-numeric data
tree <- hclust(distances, method = 'complete')

plot.title <- paste('Samples Clustered by Markers', experiment.id, sep = '\n')
plot(tree, labels = sample.ids, main = plot.title, sub = '',
     xlab = 'Distances (complete method)')
