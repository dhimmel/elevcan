library(RCurl)
library(plyr)
library(dplyr)
library(scales)
library(ggvis)


ReadFromGit <- function(...) {
  # Read a tab delimited text file from the
  # [project github repository](https://github.com/dhimmel/elevcan)
  # as a data.frame.
  raw.root <- 'https://raw.githubusercontent.com/dhimmel/elevcan/master'
  url <- file.path(raw.root, ...)
  response <- RCurl::getURL(url)
  df <- read.delim(text=response, check.names=FALSE, stringsAsFactors=FALSE, colClasses=c('fips'='character'))
  return(df)
}

# Read predictors
predictor.df <- ReadFromGit('output', 'figdata', 'predictors-best-subset.txt') %>%
  dplyr::filter(cancer == 'lung') %>%
  dplyr::arrange(dplyr::desc(abs(zcoef)))

covariates <- setdiff(predictor.df$predictor, '(Intercept)')

# Read county data
county.df <- ReadFromGit('data', 'county-data.txt') %>%
  dplyr::mutate(no_lung = all_cancer - lung)

county.select.df <- county.df %>%
  dplyr::select(fips, name, population, lung, one_of(covariates))

# Read scatterplot data
biv.df <- ReadFromGit('output', 'figdata', 'scatterplot-bivariate.txt') %>%
  dplyr::filter(Cancer == 'Lung') %>%
  dplyr::mutate(plot_type = 'bivariate') %>%
  plyr::join(county.select.df)

par.df <- ReadFromGit('output', 'figdata', 'scatterplot-partial.txt') %>%
  dplyr::filter(Cancer == 'Lung') %>%
  dplyr::mutate(plot_type = 'partial') %>%
  plyr::join(county.select.df)

mean.vec <- apply(biv.df[, covariates], 2, weighted.mean, w=biv.df$weight)


CreateHTML <- function(x) {
  # Function to create html for point tooltip with county information.

  effects <- NULL
  html.table.rows <- sapply(covariates, function(covariate) {
    coefficient <- subset(predictor.df, predictor == covariate, select='coef')
    effect <- (x[covariate] - mean.vec[covariate]) * coefficient
    effects <<- append(effects, effect)
    sprintf('<tr><td>%s</td><td>%.1f</td><td>%.1f</td></tr>',
            covariate, x[covariate], effect)
  })
  effects <- unlist(effects)
  html.table.rows <- html.table.rows[order(effects, decreasing=TRUE)]

  paste0(
    c(
      sprintf('<b>%s</b><br><i>pop. %s</i>', x['name'], x['population']),
      '<br>',
      '<br>',
      '<table>',
      '<tr><th>variable</th><th>value</th><th>effect</th></tr>',
       html.table.rows,
      '</table>'
    ), collapse='')
}


biv.df$tooltip <- plyr::daply(.data = biv.df, .variables = 'fips', .fun = CreateHTML)
par.df$tooltip <- plyr::daply(.data = par.df, .variables = 'fips', .fun = CreateHTML)

combined.df <- rbind(biv.df, par.df) %>%
  dplyr::mutate(alpha = scales::rescale(weight ^ 0.5, c(0.25, 1)), 4) %>%
  dplyr::select(fips:plot_type, alpha, tooltip) %>%
  format(digits = 3)

write.table(combined.df, 'scatterplot.txt', sep='\t', row.names=FALSE, quote=FALSE)

