library(ggvis)

shinyUI(fluidPage(
  fluidRow(column(width = 12, ggvisOutput('scatterplot'))),
  fluidRow(
    column(width = 3, selectInput(inputId = 'plot_type', label = h5('Plot Type'),
      choices = c('Bivariate'='bivariate', 'Partial Regression'='partial'))),

      column(width = 9,
             p(htmlOutput('description')))
  )
))
