#daniel.himmelstein@gmail.com

AssertPackages <- function(packages) {
  # Halts execution if any packages are not installed
  all.packages <- rownames(installed.packages())
  missing <- setdiff(packages, all.packages)
  if (length(missing) != 0) {
    missing.char <- paste(missing, collapse=', ')
    stop(sprintf('Missing packages: %s', missing.char))
  }
}

SimpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
    sep="", collapse=" ")
}

CatDiv <- function(title.chr='') {
  chr <- format(title.chr, width=80, justify='centre')
  chr <- gsub(' ', '-', chr)
  chr <- sprintf('%s\n', chr)
  cat(chr)
}

AddNewLines <- function(chr) {
  newline.chr <- gsub('(.{1,80})(\\s|$)', '\\1\n', chr)
  return(newline.chr)
}

ChrRound <- function(x, digits=2) {
  # Round x to digits and return as a character with trailing zeros.
  x.rounded <- round(x, digits)
  sprintf.chr <- sprintf('%%.%sf', digits)
  x.chr <- sprintf(sprintf.chr, x.rounded)
  return(x.chr)
}

BarometricFormula <- function(
    z,             # elevation in meters
    P0=101325,     # pressure at sea level (units don't matter)
    T0=288,        # standard temperature in K, here 288K (15C)
    B=0.0065,      # lapse rate in K/m
    m=4.8E-26,     # molecular mass of dry air in kg
    k=1.38065E-23, # Boltzmann constant in J/K
    g=9.80665      # acceleration due to gravity in m/s^2
  ) {
  # Returns atmouspheric oxygen concentration as a percentage of the sea level concentration
  # Barometric formula adjusted for vertical linear temperature changes  
  # http://dx.doi.org/10.1119/1.18555
  
  Bz <- B * z
  mg <- m * g
  kB <- k * B
  exp <- mg/kB
  
  ratio <- (1 - Bz / T0)^exp
  Pz <- P0 * ratio  
  percent <- ratio * 100
  
  return(percent)
}


ReadCountyData <- function(path, 
  states=c('AZ', 'CA', 'CO', 'ID', 'MT', 'NV', 'NM', 'OR', 'UT', 'WA', 'WY'),
  min.population=1e4, max.native=25.0, max.migration=40.0, verbose=FALSE) {
  # Read the data file with United States county level measurements and
  # perform filtering. states is a vector of state abbreviations to keep.
  # Use states='all' to suppress state filtering.
  # TODO verbose

  data.df <- read.delim(path, na.strings='', 
    colClasses=c('fips'='character'), stringsAsFactors=FALSE)

  if (verbose) {
    elev.range <- range(data.df$elevation, na.rm=TRUE)
    min.county <- data.df$name[which.min(data.df$elevation)]
    max.county <- data.df$name[which.max(data.df$elevation)]
    sprintf.chr <- 'Range of U.S. county elevation (km) without filtering: %.4f (%s) to %.4f (%s)'
    percent.delta <- 100 - 100 * BarometricFormula(elev.range[2] * 1000) / BarometricFormula(elev.range[1] * 1000)
    cat(AddNewLines(sprintf(sprintf.chr, elev.range[1], min.county, elev.range[2], max.county)))
    cat(AddNewLines(sprintf('%.2f%% decrease in oxygen concentration across elevation range', percent.delta)))
    CatDiv()
  }

  # the square root of county population is used to weight the regression.
  weight <- sqrt(data.df$population)
  weight[weight > 500] <- 500
  data.df[, 'weight'] <- weight
  if (states[1] != 'all') {data.df <- subset(data.df, state %in% states)}
  data.df <- subset(data.df, population >= min.population)
  data.df <- subset(data.df, native <= max.native)
  data.df <- subset(data.df, migration <= max.migration)

  # Calculate incidence minus specific cancer
  data.df$no_lung <- data.df$all_cancer - data.df$lung
  data.df$no_colorectal <- data.df$all_cancer - data.df$colorectal
  data.df$no_breast <- data.df$all_cancer - data.df$breast / 2
  data.df$no_prostate <- data.df$all_cancer - data.df$prostate / 2

  return(data.df)
}

ReadAllCountyData <- function(path) {
  # No filtering
  data.df <- ReadCountyData(path, states='all', min.population=0,
    max.native=100, max.migration=100)
  return(data.df)
}


WeightedScale <- function(values, value.weights) {
  # Returns standardized values (z-scores) for weighted observations.
  mean.w <- Hmisc::wtd.mean(values, value.weights)
  sd.w <- sqrt(Hmisc::wtd.var(values, value.weights, normwt=TRUE))
  return((values - mean.w) / sd.w)
}


StandardizeBeta <- function(estimate, predictor.sd, response.sd) {
  return(estimate * predictor.sd / response.sd)
}

UnstandardizeBeta <- function(zestimate, predictor.sd, response.sd) {
  return(zestimate * response.sd / predictor.sd)
}

StandardizeBetas <- function(model) {
  # Return the standardized coefficients for an lm object.
  # This method works for weighted regressions and relies on the
  # Hmisc package for calculating the weighted variance.
  # Observations with missing values should be omitted before
  # the model fitting stage.
  model.df <- model.frame(model)
  if (! any('(weights)' == colnames(model.df))) {
    model.df[, '(weights)'] <- 1
  }
  obs.weights = model.df[, '(weights)']
  response.sd <- Hmisc::wtd.var(model.df[, 1], weights=obs.weights, normwt=TRUE) ^ 0.5
  model.coefs <- coef(model)
  predictors <- names(model.coefs)[-1]
  std.betas <- sapply(predictors, function(predictor) {
    predictor.sd <- Hmisc::wtd.var(model.df[, predictor], weights=obs.weights, normwt=TRUE) ^ 0.5
    return(as.numeric(model.coefs[predictor]) / (response.sd / predictor.sd))
  })
  
  return(std.betas)
}


PartialRegressionY <- function(model, term='elevation') {
  # Residuals from regressing the response against all predictors except term
  tab <- model.frame(model)
  all.terms <- attr(terms(model), 'term.labels')
  covariates <- setdiff(all.terms, term)
  reg.formula <- as.formula(paste(colnames(tab)[1], '~', paste(covariates, collapse=' + ')))
  y.model <- lm(reg.formula, data=tab, weights=tab[, '(weights)'])
  resids <- resid(y.model)# + weighted.mean(tab[, colnames(tab)[1]], tab[, '(weights)'])
  return(resids)
}

PartialRegressionX <- function(model, term='elevation') {
  # Residuals from regressing term against the remaining covariates
  tab <- model.frame(model)
  all.terms <- attr(terms(model), 'term.labels')
  covariates <- setdiff(all.terms, term)
  reg.formula <- as.formula(paste(term, '~', paste(covariates, collapse=' + ')))
  x.model <- lm(reg.formula, data=tab, weights=tab[, '(weights)'])
  resids <- resid(x.model)
  return(resids)
}


