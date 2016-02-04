library(cluster)
library(ggplot2)
library(ggrepel)
library(limma)
library(plyr)
library(reshape)

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

genotypes.transposed <- as.data.frame(t(genotypes))
genotypes.transposed <- colwise(as.factor)(genotypes.transposed)

# gower for non-numeric data
distances <- daisy(genotypes.transposed, metric = 'gower')
tree <- hclust(distances, method = 'complete')

plot.title <- paste('Samples Clustered by Markers', experiment.id, sep = '\n')
plot(tree, labels = sample.ids, main = plot.title, sub = '',
     xlab = 'Distances (complete method)')

clustered.sample.order <- tree$order
genotypes.sorted <- genotypes[, clustered.sample.order]
genotypes.sorted$marker.id <- rownames(genotypes.sorted)

bin.genotypes.m <- melt(genotypes.sorted, id = c('marker.id'),
                        variable_name = 'sample')
bin.genotypes.m$value <- as.factor(bin.genotypes.m$value)

allele.colors <- c('skyblue', 'orange', 'black', 'green', 'yellow', 'plum')

ggplot(bin.genotypes.m) +
  geom_raster(aes(x = marker.id, y = sample, fill = value)) +
  ggtitle('Samples Clustered By Genotype') +
  xlab('Markers') +
  ylab('Samples') +
  labs(fill = 'Genotype') +
  scale_fill_manual(values = allele.colors) +
  theme(axis.text.x = element_text(angle = 90))
