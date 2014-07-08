# Daniel Himmelstein
# daniel.himmelstein@gmail.com

dropbox.dir <- '~/Dropbox/lung/'
mapping.dir <- file.path(dropbox.dir, 'data', 'mappings')

PadZeros <- function(x, length.out) {
  # Prepend zeros to x (or each element of vector x) to reach a character
  # number of length.out.
  # x must be capable of being converted by as.numeric(). 
  x <- as.numeric(x)
  x <- formatC(x, width=length.out, format='d', flag='0')
  return(x)
}

RemoveCommas <- function(s) {
  # Remove commas from string s
  gsub(",", "", s, fixed = TRUE) 
}

Trim <- function(s) {
  # strip leading and trailing whitespace from string s
  gsub("^[[:space:]]+|[[:space:]]+$", "", s)
}

TrimTable <- function(tab) {
  # Trim every element of every column of a table
  for (j in 1:ncol(tab)) {
    tab[, j] <- Trim(tab[, j])
  }
  return(tab)
}

# Read State information
state.mapping.path <- file.path(mapping.dir, 'states.txt')
state.tab <- read.delim(state.mapping.path, colClasses='character')
fips2state.name <- state.tab$name; names(fips2state.name) <- state.tab$fips
fips2state.abbrev <- state.tab$abbrev; names(fips2state.abbrev) <- state.tab$fips
state.abbrev2fips <- state.tab$fips; names(state.abbrev2fips) <- state.tab$abbrev
state.name2fips <- state.tab$fips; names(state.name2fips) <- state.tab$name

StateName2Abbrev <- function(state.names) {
  state.abbrevs <- state.abb[match(state.names, state.name)]
  return(state.abbrevs)
}

StateAbbrev2Name <- function(state.abbrevs) {
  state.names <- state.name[match(state.abbrevs, state.abbrevs)]
  return(state.names)
}


# Read county information
county.mapping.path <- file.path(mapping.dir, 'counties.txt')
county.tab <- read.delim(county.mapping.path, colClasses='character')
county.tab$fips <- paste(county.tab$state_fips, county.tab$county_fips, sep='')
county.tab$lowercase_county_name = tolower(county.tab$county_name)

ReadStateCancerProfileCSV <- function(path, skip=5) {
  # Read an NCI State Cancer Profile CSV file.
  # For 2005-2009 csv files use skip=6
  # For 2006-2010 csv files use skip=5
  # http://www.statecancerprofiles.cancer.gov/map/map.noimage.php
  scp.tab <- read.table(path, sep=',', comment.char='"',
    skip=skip, na.strings =c('\xa7', '\xb6', '*', 'Data Not Available'), header=TRUE,
    colClasses='character', strip.white=TRUE, quote='')
  return(scp.tab)
}

ReadStateCancerProfileIncidence <- function(path) {
  print(path)
  scp.tab <- ReadStateCancerProfileCSV(path)
  if(any(colnames(scp.tab) == 'Annual.Incidence.Rate')) {
  # Remove ' # ' string appended to incidence which indicate diognoses in other states omitted.
  scp.tab$Annual.Incidence.Rate <- gsub(pattern="[# ]", '', scp.tab$Annual.Incidence.Rate)
  }
  cancer.tab <- data.frame('fips'=scp.tab$FIPS,
    'incidence'=scp.tab$Annual.Incidence.Rate,
    'lower.95'=scp.tab$Lower.95..Confidence.Interval,
    'upper.95'=scp.tab$Upper.95..Confidence.Interval,
    'annual.cases'=scp.tab$Average.Annual.Cases, stringsAsFactors=FALSE)
  for (j in 2:5) suppressWarnings(cancer.tab[, j] <- as.numeric(cancer.tab[, j]))
  return(cancer.tab)
}

ReadSAE <- function(path, skip.first, skip.last) {
  # Read the csv file for NCI Small Area Estimates export.
  # http://sae.cancer.gov/estimates/tables/both_lifetime.html
  tab <- read.table(path, header=TRUE, skip=3, sep=',', quote='"',
    na.strings =c('!', '^'), colClasses='character',
    strip.white=TRUE, comment.char='', fill=TRUE,  blank.lines.skip=FALSE)
  # remove bottom skip.last rows
  tab <- tab[seq(nrow(tab) - skip.last), ]
  tab <- TrimTable(tab)
  return(tab)
}


ReadStateTobaccoSAE <- function(path) {
  tab <- ReadSAE(path, 3, 9)
  smoking <- mapply(function(x, y) mean(c(x, y), na.rm=TRUE),
    as.numeric(tab$Model.Based.Estimate....),
    as.numeric(tab$Model.Based.Estimate.....1))
  smoking.tab <- data.frame('state'=tab[, 1], 'smoking'=smoking, stringsAsFactors=FALSE)
  return(smoking.tab)
}

ReadCountyTobaccoSAE <- function(path, printing=FALSE) {
  # Read the csv file for tobacco information from the NCI Small Area Estimates
  # export.
  # http://sae.cancer.gov/estimates/tables/both_lifetime.html
  tobacco.tab <- ReadSAE(path, 3, 10)
  if (printing) cat('Counties in the NCI SAE tobacco data which do not map to our counties:\n')
  tobacco.fips <- sapply(tobacco.tab$County, function(county) {
    county.state <- strsplit(county, ', ', fixed=TRUE)[[1]]
    state.fips <- as.character(state.abbrev2fips[substr(county.state[2], 1, 2)])
    county.name <- county.state[1]
    county.name <- sub(" [(]County[)]$", '', county.name)
    county.name <- gsub("[.]", '', county.name)
    county.name <- tolower(county.name)
    keep.row <- county.tab$state_fips == state.fips & county.tab$lowercase_county_name == county.name
    if (any(keep.row)) {
      fips <- county.tab[keep.row, 'fips']
    } else {
      fips <- NA
      if (printing) {cat(county); cat('\n')}
    }
    return(fips)
  })

  smoking <- mapply(function(x, y) mean(c(x, y), na.rm=TRUE), as.numeric(tobacco.tab$Model.Based.Estimate....),
    as.numeric(tobacco.tab$Model.Based.Estimate.....1))
  smoking[is.nan(smoking)] <- NA

  smoking.tab <- data.frame('fips'=tobacco.fips, 'smoking'=smoking, stringsAsFactors=FALSE)
  smoking.tab <- smoking.tab[order(smoking.tab$fips), ]
  return(smoking.tab)
}

ReadCDCWonder <- function(path) {
  # http://wonder.cdc.gov/
  # Must replace superscripts and save in an encoding R can read for some files.
  tab <- read.delim(path, colClasses='character', na.strings=c('', 'Missing'))
  keep.row <- ! apply(is.na(tab[, -1]), 1, all)
  tab <- tab[keep.row, ]
  return(tab)
}

ReadBRFSS <- function(path){
  # Read SAS .xpt file from the Behavioral Risk Factor Surveillance System
  library('foreign')
  xpt.path <- file.path('Downloads', 'brfss', '2011', 'LLCP2011.XPT')
  brfss <- read.xport(xpt.path)
  dim(brfss) # 504408 x 450
}

ReadHRP <- function(path) {
  # Read a high-radon project file ending in .CT
  # http://energy.lbl.gov/ie/high-radon/data/lbnl-met.html
  fields <- c('longitude', 'latitude', 'state', 'county', 'fips', month.name, 'yearly')
  tab <- read.table(path, col.names=fields, 
    stringsAsFactors=FALSE, colClasses=c('fips'='character'))
  tab$fips <- PadZeros(tab$fips, length.out=5)
  return(tab)
}



