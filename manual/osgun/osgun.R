library(ggplot2)
library(gridExtra)
library(RCurl)
library(plyr)
library(dplyr)

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

# Read dataset with counties included for lung cancer models, with dwyer smoking data added.
lung.df <- read.delim('lung-table-with-dwyer.txt', stringsAsFactors=FALSE)
lung.df <- na.omit(lung.df)

SetGGTheme <- function(gg) {
  # Function to modify ggplot theme to mansucript style
  gg <- gg + theme_bw()
  gg <- gg + theme(plot.margin=grid::unit(c(2, 2, 2, 2), 'points'))
  return(gg)
}

# Scatterplot options
biv.line.col   <- '#990000'
biv.fill.col   <- '#ffb2b2'

# Scatterplot -- County smoking in 2012 versus smoking in 1996
gg1 <- ggplot(lung.df, aes(smoking_1996, smoking_2012)) 
gg1 <- SetGGTheme(gg1) +
  geom_smooth(method='lm', aes(weight=weight), level=0.95,
    color=biv.line.col, fill=biv.fill.col, alpha=1) +
  geom_point(aes(alpha=weight ^ 2)) +
  scale_alpha(range=c(0.2, 1)) + guides(alpha=FALSE) +
  xlab('Smoking in 1996') + ylab('Smoking in 2012')

# Scatterplot -- Smoking in 2012 versus smoking in 1996
gg2 <- ggplot(lung.df, aes(elevation, dwyer_delta)) 
gg2 <- SetGGTheme(gg2) +
  geom_smooth(method='lm', aes(weight=weight), level=0.95,
    color=biv.line.col, fill=biv.fill.col, alpha=1) +
  geom_point(aes(alpha=weight ^ 2)) +
  scale_alpha(range=c(0.2, 1)) + guides(alpha=FALSE) +
  xlab('Elevation (km)') + ylab('Change in Smoking, 1996-2012')

# All variables included as predictors in lung cancer analyses plus dwyer's smoking deltas
global.covars <- c('metro', 'white', 'black', 'education', 'income', 'obesity')
lung.covars <- c('no_lung', 'male', 'smoking', 'radon', 'particulate')
envir.covars <- c('elevation', 'uvb', 'sunlight', 'precipitation', 'high_temp', 'diurnal_temp', 'radon', 'particulate')
variables <- unique(c('dwyer_delta', global.covars, lung.covars, envir.covars))

# Read entire county-level data for all variables.
county.df <- ReadFromGit('data', 'county-data.txt')

# Calculate variable correlations with dywer's smoking delta for counties in lung cancer model
joined.df <- plyr::join(lung.df, county.df)
cor.vec <- cor(joined.df[, variables], use='pairwise.complete.obs')[, 'dwyer_delta'][-1]
cor.vec <- cor.vec[order(abs(cor.vec), decreasing=TRUE)]
cor.df <- data.frame(variable=names(cor.vec), correlation=cor.vec)


# Tile Plot Options
div.cols <- c('#0000FF', '#f7eff7', '#FF0000')
ylimits <- rev(cor.df$variable)
ylabels <- gsub('_', ' ', ylimits)

# Plot variable correlation with dywer's smoking delta
gg3 <- SetGGTheme(ggplot(cor.df, aes(x=0, y=variable, fill=correlation))) +
  geom_tile(color='white') +
  scale_fill_gradientn(limits=c(-1, 1), colours=div.cols) + 
  xlab(NULL) + ylab(NULL) + guides(fill=FALSE) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    theme(plot.background=element_blank(), panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.line=element_line(color='grey50', size=0.25)) +
    scale_y_discrete(limits=ylimits, labels=ylabels)

# Save three-panel figure as a raster
png('smoking-lagtime.png', width=7, height=3 , bg='white', family='sans', units='in', res=300)
gridExtra::grid.arrange(gg1, gg2, gg3, nrow=1, widths=c(2.5, 2.5, 1))
dev.off()

# Save three-panel figure as a vector
cairo_pdf('smoking-lagtime.pdf', width=7, height=3 , bg='white', family='sans')
gridExtra::grid.arrange(gg1, gg2, gg3, nrow=1, widths=c(2.5, 2.5, 1))
dev.off()

