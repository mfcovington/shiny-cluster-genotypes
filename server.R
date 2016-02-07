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
      df <- read.table(file = input.file$datapath, sep = input$delimiter,
                       header = TRUE, row.names = 1, na.strings = na.strings)

      statusReport(
        sprintf("Reading file '%s' (%i rows x %i columns).\n",
                input.file$name, nrow(df), ncol(df)))

      df
    }
  })

  output$genotypes <- renderTable({
    df <- getData()
    df[1:min(c(5, nrow(df))), 1:min(c(5, ncol(df)))]
  })

})
