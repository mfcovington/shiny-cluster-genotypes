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

  sampleIds <- reactive({
    getSampleIds(getData())
  })

  output$genotypes <- renderDataTable({
    getData()
  }, options = list(pageLength = 10))

  output$mds.plot <- renderPlot({
    plotGenotypeMDS(getData(), sampleIds())
  })

  makeTree <- reactive({makeGenotypeTree(getData())})

  output$tree.plot <- renderPlot({
    plotGenotypeTree(makeTree(), sampleIds(), input$experiment.id)
  })

  output$tile.plot <- renderPlot({
    genotypes.sorted <- sortGenotypesByTree(getData(), makeTree())
    allele.colors = unlist(strsplit(input$allele.colors, ','))
    plotGenotypeTile(genotypes.sorted, allele.colors)
  })

  output$genotypes_sorted <- renderDataTable({
    sortGenotypesByTree(getData(), makeTree())
  }, options = list(pageLength = 10))
})
