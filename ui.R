library(shiny)

shinyUI(pageWithSidebar(

  headerPanel('Clustering Genotype Data'),

  sidebarPanel(

    textInput('experiment.id', 'Experiment ID:', ''),

    fileInput('input.file', 'Input Data File:'),

    selectInput('delimiter', 'File Type:',
                list('Tab-delimited' = '\t',
                     'Comma-separated values' = ',')),

    textInput('na.strings',
              'Strings that represent missing values in the data set (if >1, separate with commas):',
              ''),

    textInput('allele.colors', 'Color Palette (must have at least as many colors as the maximum number of alleles for a marker; separate with commas):',
              'skyblue,orange,black,green,yellow,plum')

  ),

  mainPanel()
))
