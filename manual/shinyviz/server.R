library(ggvis)

# Read data
data.df <- read.delim('scatterplot.txt',
  stringsAsFactors=FALSE, colClasses=c('fips'='character'))

# Color options
stroke.colors <- c(bivariate = '#990000', partial = '#000099')
fill.colors <- c(bivariate = '#ffb2b2', partial = '#b2b2ff')
fig.panels <- c(bivariate = 'A', partial = 'B')

xlabels <- c(bivariate = 'Elevation (kilometers)', partial = 'Elevation Residuals')
ylabels <- c(bivariate = 'Lung Cancer Incidence (per 100,000)', partial = 'Lung Cancer Incidence Residuals')

axis.properties <- axis_props(title = list(fontSize = 15))

shinyServer(function(input, output, session) {

  scatter.df <- reactive({
    subset(data.df, plot_type == input$plot_type)
  })

  scatterplot <- reactive({
    xlabel <- as.character(xlabels[input$plot_type])
    ylabel <- as.character(ylabels[input$plot_type])
    #stroke.color <- stroke.colors[input$plot_type]
    #fill.color <- fill.colors[input$plot_type]

    scatter.df %>%
      ggvis(x = ~Elevation, y = ~Incidence) %>%
      layer_model_predictions(model = 'lm', se = TRUE,
        stroke := '#256325', strokeOpacity := 1,
        fill := '#4EA34E', fillOpacity := 1) %>%
      layer_points(opacity := ~alpha, key := ~tooltip) %>%
      add_axis(type = 'x', title = xlabel, properties = axis.properties) %>%
      add_axis(type = 'y', title = ylabel, properties = axis.properties) %>%
      add_tooltip(function(x) {x$tooltip}, on = 'click') %>%
      set_options(width = 510, height = 400)
  })

  scatterplot %>%
    bind_shiny('scatterplot')

  output$description <- renderText({
    paste0(
      'Points represent Western US counties shaded according to population. ',
      'Clicking a point shows county information. Deselect by clicking the 95% confidence band. ',
      'In the county information tables, <code>effect</code> indicates the estimated effect ',
      'that each variable had on incidence (in relation to the variable\'s mean). ',
      '',
      'This plot corresponds to <a href="http://dx.doi.org/10.7717/peerj.705/fig-4" target="_blank">manuscript Fig. 4',
      fig.panels[input$plot_type], '</a>.'

    )
    })

})
