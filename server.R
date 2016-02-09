source('bin/cluster.R')
library(shiny)

shinyServer(function(input, output) {

  statusReport <- function(status.msg, app.out = TRUE, console.out = TRUE) {
    if (console.out)
      cat(status.msg)

    if (app.out)
      output$status <- renderText(status.msg)
  }

  getData <- eventReactive(input$run, {
    withProgress(message = 'Getting Data', {
      input.file <- input$input.file

      if (is.null(input.file)) {
        NULL
      } else {
        na.strings = unlist(strsplit(input$na.strings, ','))
        df <- getDataFromFile(input.file$datapath, input$delimiter, na.strings)

        statusReport(
          sprintf("Reading file '%s': %i markers (rows) x %i samples (columns).\n",
                  input.file$name, nrow(df), ncol(df)))

        df
      }
    })
  })

  sampleIds <- reactive({
    getSampleIds(getData())
  })

  output$genotypes <- renderDataTable({
    genotypes <- getData()
    cbind(marker.ids = rownames(genotypes), genotypes)
  }, options = list(pageLength = 10))

  output$mds.plot <- renderPlot({
    withProgress(message = 'Building MDS Plot', {
      plotGenotypeMDS(getData(), sampleIds(),
                      verbose.labels = input$verbose.labels)
    })
  })

  makeTree <- reactive({makeGenotypeTree(getData())})

  output$tree.plot <- renderPlot({
    withProgress(message = 'Building Tree Plot', {
      plotGenotypeTree(makeTree(), sampleIds(), input$experiment.id)
    })
  })

  output$tile.plot <- renderPlot({
    withProgress(message = 'Building Tile Plot', {
      genotypes.sorted <- sortGenotypesByTree(getData(), makeTree())
      allele.colors = unlist(strsplit(input$allele.colors, ','))
      plotGenotypeTile(genotypes.sorted, allele.colors)
    })
  })

  sortData <- reactive({
    withProgress(message = 'Sorting Genotype Data', {
      sortGenotypesByTree(getData(), makeTree())
    })
  })

  output$genotypes_sorted <- renderDataTable({sortData()},
                                              options = list(pageLength = 10))

  output$downloadSortedData <- downloadHandler(
    filename = function() {
      sub('(.+)\\.([^.]+)$', '\\1.sorted.\\2', input$input.file$name)
    },
    content = function(file) {
      writeSortedGenotypes(sortData(), input$delimiter, file)
    }
  )

})
