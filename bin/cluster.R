library(cluster)
library(ggplot2)
library(ggrepel)
library(limma)
library(plyr)
library(reshape)


######################################################
# PROCESS SAMPLE DATA                                #
######################################################
# 1. Import all functions and libraries in this file #
#    or run: source('cluster.R')                     #
# 2. Run: processSampleData()                        #
######################################################

processSampleData <- function() {
  experiment.id <- 'SNP Data'
  input.file <- 'sample-data.snps.txt'
  delimiter <- '\t'
  na.strings <- c('N', 'U')
  allele.colors <- c('skyblue', 'orange', 'black', 'green', 'yellow', 'plum')

  genotypes <- getDataFromFile(input.file, delimiter, na.strings)
  sample.ids <- getSampleIds(genotypes)

  plotGenotypeMDS(genotypes, sample.ids)

  tree <- makeGenotypeTree(genotypes)
  plotGenotypeTree(tree, sample.ids, experiment.id)

  genotypes.sorted <- sortGenotypesByTree(genotypes, tree)
  plotGenotypeTile(genotypes.sorted, allele.colors)

  writeSortedGenotypes(genotypes.sorted, delimiter)
}


############
# GET DATA #
############

getDataFromFile <- function(input.file, delimiter, na.strings) {
  read.table(file = input.file, sep = delimiter, header = TRUE, row.names = 1,
             na.strings = na.strings)
}

getSampleIds <- function(genotypes) {
  colnames(genotypes)
}


############
# MDS PLOT #
############

plotGenotypeMDS <- function(genotypes, sample.ids) {
  mds <- plotMDS(data.matrix(genotypes))
  mds.df <- as.data.frame(mds$cmdscale.out)
  mds.df$idx <- 1:length(sample.ids)

  p <- ggplot(mds.df) +
    geom_point(aes(x = V1, y = V2)) +
    geom_label_repel(aes(x = V1, y = V2, label = idx)) +
    ggtitle(expression('Leading ' * Log[2] * ' fold change MDS')) +
    xlab('Dimension 1') +
    ylab('Dimension 2')

  print(p)
}


#############
# TREE PLOT #
#############

makeGenotypeTree <- function(genotypes) {
  genotypes.transposed <- as.data.frame(t(genotypes))
  genotypes.transposed <- colwise(as.factor)(genotypes.transposed)

  # gower for non-numeric data
  distances <- daisy(genotypes.transposed, metric = 'gower')
  tree <- hclust(distances, method = 'complete')
}

plotGenotypeTree <- function(tree, sample.ids, experiment.id) {
  plot.title <- paste('Samples Clustered by Markers', experiment.id, sep = '\n')
  plot(tree, labels = sample.ids, main = plot.title, sub = '',
       xlab = 'Distances (complete method)')
}


#############
# TILE PLOT #
#############

sortGenotypesByTree <- function(genotypes, tree) {
  clustered.sample.order <- tree$order
  genotypes.sorted <- genotypes[, clustered.sample.order]
  cbind(marker.id = rownames(genotypes.sorted), genotypes.sorted)
}

plotGenotypeTile <- function(genotypes.sorted, allele.colors) {
  bin.genotypes.m <- melt(genotypes.sorted, id = c('marker.id'),
                          variable_name = 'sample')
  bin.genotypes.m$value <- as.factor(bin.genotypes.m$value)
  bin.genotypes.m$marker.id <- factor(bin.genotypes.m$marker.id,
                                      levels = rownames(genotypes.sorted))

  p <- ggplot(bin.genotypes.m) +
    geom_raster(aes(x = marker.id, y = sample, fill = value)) +
    ggtitle('Samples Clustered By Genotype') +
    xlab('Markers') +
    ylab('Samples') +
    labs(fill = 'Genotype') +
    scale_fill_manual(values = allele.colors) +
    theme(axis.text.x = element_text(angle = 90))

  print(p)
}


################################
# SAVE DATA SORTED BY GENOTYPE #
################################

writeSortedGenotypes <- function(genotypes.sorted, delimiter) {
  write.table(genotypes.sorted[, 2:ncol(genotypes.sorted)],
              file = 'genotypes-sorted.txt', sep = delimiter, quote = FALSE,
              col.names = NA)
}
