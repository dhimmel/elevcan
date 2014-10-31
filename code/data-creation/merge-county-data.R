# Daniel Himmelstein
# daniel.himmelstein@gmail.com
# Merge all data sources into one dataset where each row is a county
dropbox.dir <- '~/Dropbox/lung'
code.dir <- file.path(dropbox.dir, 'analysis', 'code', 'data-creation')
source(file.path(code.dir, 'source-reader.R'))

options(stringsAsFactors=FALSE)

demography.dir <- file.path(dropbox.dir, 'data', 'county', 'demography')
environment.dir <- file.path(dropbox.dir, 'data', 'county', 'environment')
cancer.dir <- file.path(dropbox.dir, 'data', 'county', 'cancer')
cancer.time.period <- '2005-2009'
risk.dir <- file.path(dropbox.dir, 'data', 'county', 'risk')

for (extras in c(TRUE, FALSE)) {

cat(sprintf('Computing with extras: %s\n', extras))

# Create a table with county name and fips id as a template to variables onto
# check that county tab has every county
data.tab <- data.frame(
  'name'= sprintf('%s, %s', county.tab$county_name, StateName2Abbrev(county.tab$state_name)),
  'state'=StateName2Abbrev(county.tab$state_name),
  'fips'=county.tab$fips)

# Add 2004 USDA typology information
path <- file.path(demography.dir, 'usda', 'typology-2004.csv')
tab <- read.csv(path, colClasses='character')
indices <- match(data.tab$fips, tab$FIPSTXT)
data.tab[, 'metro'] <- tab[indices, 'metro']

# Add immigration
path <- file.path(demography.dir, 'migration', 'gross-migration.csv')
col.names <- c('fips', 'name', 'population', 'nonmovers', 'movers',
  'movers.same.area', 'movers.abroad')
tab <- read.csv(path, header=FALSE, skip=79, nrows=3237 - 79, 
  col.names=col.names, strip.white=TRUE, colClasses='character', na.strings=c('', 'X', '-'))
tab <- tab[rowSums(is.na(tab)) < 2, ] # remove state divider rows
tab$fips <- PadZeros(tab$fips, 5)
for (j in 1:ncol(tab)) {
  tab[, j] <- RemoveCommas(tab[, j])
  if (j > 2) tab[, j] <- as.numeric(tab[, j])
}
tab$immigration <- (tab$movers - tab$movers.same.area) / tab$population
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'immigration'] <- tab[indices, 'immigration'] 
data.tab[, 'immigration'] <- as.numeric(data.tab[, 'immigration']) * 100 # convert to a percentage

# Add Race Percentages

# Native Americans
path <- file.path(demography.dir, 'race', 'native-combo.txt')
tab <- read.delim(path, colClasses='character')
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'native'] <- tab[indices, 'percent.race']

# whites
path <- file.path(demography.dir, 'race', 'white-combo.txt')
tab <- read.delim(path, colClasses='character')
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'white'] <- tab[indices, 'percent.race']

# blacks
path <- file.path(demography.dir, 'race', 'black-combo.txt')
tab <- read.delim(path, colClasses='character')
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'black'] <- tab[indices, 'percent.race']

if (extras) {
# Asians
path <- file.path(demography.dir, 'race', 'asian-combo.txt')
tab <- read.delim(path, colClasses='character')
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'asian'] <- tab[indices, 'percent.race']
}

# male: percent male
path <- file.path(demography.dir, 'national-atlas', 'ce2000t.tdt')
tab <- read.delim(path, colClasses='character', check.names=FALSE)
indices <- match(data.tab$fips, tab[, 'FIPS,C,5'])
males <- as.numeric(tab[indices, 'MALE2000,N,11,0'])
females <- as.numeric(tab[indices, 'FEMALE2000,N,10,0'])
data.tab[, 'male'] <- 100 * males / (males + females)

# Add population and elevation
elev.path <- file.path(environment.dir, 'elevation-blockgroup-averaged-2000.txt')
elevation.tab <- read.delim(elev.path, colClasses='character')
indices <- match(data.tab$fips, elevation.tab$county)
data.tab$population <- elevation.tab$population[indices]
data.tab[, 'elevation'] <- elevation.tab[indices, 'elevation']
data.tab[, 'elevation'] <- as.numeric(data.tab[, 'elevation']) / 1000 # convert to km

# Add barometric pressure (adjusted to sea level and poorly documented)
if (extras) {
path <- file.path(environment.dir, 'high-radon-project', 'PRESAVE.CT')
tab <- ReadHRP(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'pressure'] <- tab[indices, 'yearly']
}

# Add temperature
path <- file.path(environment.dir, 'temperature.txt')
tab <- ReadCDCWonder(path)
indices <- match(data.tab$fips, tab$County.Code)
data.tab[, 'high_temp'] <- tab[indices, 'Avg.Daily.Max.Air.Temperature.C.']
data.tab[, 'low_temp'] <- tab[indices, 'Avg.Daily.Min.Air.Temperature.C.']
data.tab[, 'diurnal_temp'] <- as.numeric(data.tab[, 'high_temp']) - as.numeric(data.tab[, 'low_temp'])
if (extras) {
data.tab[, 'heat_index'] <- tab[indices, 'Avg.Daily.Max.Heat.Index.C.']
}

# Add fine particulate matter
path <- file.path(environment.dir, 'fine-particulate-matter.txt')
tab <- ReadCDCWonder(path)
indices <- match(data.tab$fips, tab$County.Code)
data.tab[, 'particulate'] <- tab[indices, 'Avg.Fine.Particulate.Matter.µg.m3.']

# Add precipitation
path <- file.path(environment.dir, 'precipitation.txt')
tab <- ReadCDCWonder(path)
indices <- match(data.tab$fips, tab$County.Code)
data.tab[, 'precipitation'] <- tab[indices, 'Avg.Daily.Precipitation.mm.']

# Add sunlight
path <- file.path(environment.dir, 'sunlight.txt')
tab <- ReadCDCWonder(path)
indices <- match(data.tab$fips, tab$County.Code)
data.tab[, 'sunlight'] <- tab[indices, 'Avg.Daily.Sunlight.KJ.m3.']

# Add UVB
path <- file.path(environment.dir, 'UV-B_exposure.txt')
tab <- read.delim(path, comment.char='#', col.names=c('fips', 'uvb'), colClasses='character')
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'uvb'] <- tab[indices, 'uvb']

# Add lifetime smoking
tobacco.lifetime.path <- file.path(risk.dir, 'tobacco', 'lifetime', 'all.csv')
lifetime.smoking.tab <- ReadCountyTobaccoSAE(tobacco.lifetime.path)
indices <- match(data.tab$fips, lifetime.smoking.tab$fips)
data.tab[, 'smoking'] <- lifetime.smoking.tab[indices, 'smoking']
data.tab[, 'smoking_lower'] <- lifetime.smoking.tab[indices, 'smoking.lower']
data.tab[, 'smoking_upper'] <- lifetime.smoking.tab[indices, 'smoking.upper']

# Add lifetime smoking for females
path <- file.path(risk.dir, 'tobacco', 'lifetime', 'female.csv')
tab <- ReadCountyTobaccoSAE(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'female_smoking'] <- tab$smoking[indices]

# Add lifetime smoking for males
path <- file.path(risk.dir, 'tobacco', 'lifetime', 'male.csv')
tab <- ReadCountyTobaccoSAE(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'male_smoking'] <- tab$smoking[indices]

# Add current smoking
if (extras) {
tobacco.current.path <- file.path(risk.dir, 'tobacco', 'current', 'all.csv')
current.smoking.tab <- ReadCountyTobaccoSAE(tobacco.current.path)
data.tab[, 'current_smoking'] <- current.smoking.tab$smoking[match(data.tab$fips, current.smoking.tab$fips)]
}

# Read air chek radon data
if (extras) {
path <- file.path(risk.dir, 'radon-airchek', 'parsed-merged-radon.txt')
tab <- read.delim(path, colClasses='character')
tab$name <- paste(tab$county, ', ', tab$state, sep='')
tab$modified.name <- tolower(tab$name)
tab$modified.name <- sub('saint ', 'st ', tab$modified.name, fixed=TRUE)
#tab$modified.name <- sub(' county', '', tab$modified.name, fixed=TRUE)
#data.tab$name[! (tolower(data.tab$name)  %in% tab$modified.name)]
#tab$name[! (tab$lower.name %in% tolower(data.tab$name))]
indices <- match(tolower(data.tab$name), tab$modified.name)
data.tab[, 'aircheck_radon'] <- tab$average_level[indices]
}

# Add cohen radon data
# cohen radon data exists for 1601 counties
if (extras) {
path <- file.path(risk.dir, 'radon-cohen', 'modified.dat')
tab <- read.table(path, colClasses='character')
tab$fips <- PadZeros(tab$V1, 5)
tab$radon <- tab$V72
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'cohen_radon'] <- tab[indices, 'radon']
}

# Read LBL radon data
path <- file.path(risk.dir, 'radon-lbl.txt')
tab <- read.table(path, colClasses='character', header=TRUE)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'radon'] <- tab[indices, 'AMest']

# Read county health reports
path <- file.path(dropbox.dir, 'data', 'county', 'county-health-reports', '2010-extracted.csv')
tab <- read.csv(path, colClasses='character')
indices <- match(data.tab$fips, tab$FIPS)
if (extras) {
data.tab[, 'ozone_days'] <- tab$ozone[indices]
data.tab[, 'particulate_days'] <- tab$particulate[indices]
data.tab[, 'liquor_density'] <- tab$liquor_density[indices]
}
data.tab[, 'drinking'] <- tab[indices, 'drinking']

if (extras) {
path <- file.path(dropbox.dir, 'data', 'county', 'county-health-reports', '2010-mortality-morbidity.txt')
tab <- read.delim(path, colClasses='character')
indices <- match(data.tab$fips, tab$FIPS)
#Age-adjusted years of potential life lost (YPLL) rate per 100,000 (2004–2006)
data.tab[, 'mortality'] <- tab[indices, 'YPLL_Rate']
data.tab[, 'low_birthweight'] <- tab[indices, 'low_birthweight']
data.tab[, 'bad_health'] <- tab[indices, 'bad_health']
data.tab[, 'unhealthy_days'] <- tab[indices, 'unhealthy_days']
data.tab[, 'mental_days'] <- tab[indices, 'mental_days']
}


# Read Food Environment Atlas
path <- file.path(dropbox.dir, 'data', 'county', 'risk', 'food-environment-atlas-2011.txt')
tab <- read.delim(path, colClasses='character', na.strings='-9999')
indices <- match(data.tab$fips, tab[, 'FIPSTXT'])
data.tab[, 'meat'] <- tab[indices, 'PH_MEAT']
if (extras) {
data.tab[, 'produce'] <- tab[indices, 'PH_FRUVEG']
data.tab[, 'fats'] <- tab[indices, 'PH_FATS']
data.tab[, 'soda'] <- tab[indices, 'PH_SODA']
data.tab[, 'snacks'] <- tab[indices, 'PH_SNACKS']
}


# Read Diabetes
path <- file.path(dropbox.dir, 'data', 'county', 'risk', 'diabetes-ade-adjusted.txt')
tab <- read.delim(path, colClasses='character', na.strings='No Data', check.names=FALSE)
indices <- match(data.tab$fips, tab[, 'FIPS Codes'])
data.tab[, 'diabetes'] <- apply(data.matrix(tab[indices, c('2004', '2005', '2006', '2007', '2008')]), 1, mean, na.rm=TRUE)


# Read obesity data
path <- file.path(risk.dir, 'obesity2004.csv')
tab <- read.csv(path, colClasses='character')
tab$fips <- paste(PadZeros(tab$STATE_FIPS, 2), PadZeros(tab$COUNTY_FIPS, 3), sep='')
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'obesity'] <- tab[indices, 'ADJPERCENT']

# Read inactivity data
if (extras) {
path <- file.path(risk.dir, 'inactivity2004.csv')
tab <- read.csv(path, colClasses='character')
tab$fips <- paste(PadZeros(tab$STATE_FIPS, 2), PadZeros(tab$COUNTY_FIPS, 3), sep='')
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'inactivity'] <- tab[indices, 'ADJPERCENT']
}

# Read education - percent bachelors - info
path <- file.path(demography.dir, 'money', 'education-bachelor.csv')
tab <- ReadStateCancerProfileCSV(path)
indices <- match(data.tab$fips, tab$FIPS)
data.tab[, 'education'] <- tab[indices, 'Percent']

# Read education - percent highschool grad - info
if (extras) {
path <- file.path(demography.dir, 'money', 'education-highschool.csv')
tab <- ReadStateCancerProfileCSV(path)
indices <- match(data.tab$fips, tab$FIPS)
data.tab[, 'highschool_education'] <- tab[indices, 'Percent']
}

# Read median household income info
path <- file.path(demography.dir, 'money', 'income.csv')
tab <- ReadStateCancerProfileCSV(path)
indices <- match(data.tab$fips, tab$FIPS)
data.tab[, 'income'] <- tab[indices, 'Percent']
data.tab[, 'income'] <- as.numeric(data.tab[, 'income']) / 1000 # convert to thousands of dollars


# Read poverty percentage
if (extras) {
path <- file.path(demography.dir, 'money', 'poverty.csv')
tab <- ReadStateCancerProfileCSV(path)
indices <- match(data.tab$fips, tab$FIPS)
data.tab[, 'poverty'] <- tab[indices, 'Percent']
}

# Read mammogram info
path <- file.path(risk.dir, 'mammogram.csv')
tab <- ReadStateCancerProfileCSV(path)
indices <- match(data.tab$fips, tab$FIPS)
data.tab[, 'mammogram'] <- tab[indices, 'Model.Based.Percent']


# Add coal mining information
if (FALSE) {
mining.path <- file.path(demography.dir, 'mining', 'mine-number-and-production.csv')
mining.tab <- read.csv(mining.path, stringsAsFactors=FALSE)
name <- paste(mining.tab$county_name, ', ', StateName2Abbrev(mining.tab$state_name), sep='')
indices <- match(tolower(data.tab$name), tolower(name))
stopifnot(all(tolower(name) %in% tolower(data.tab$name)))
data.tab$coal.surface <- mining.tab$surface_production[indices]
data.tab$coal.underground <- mining.tab$underground_production[indices]
data.tab$coal.total <- mining.tab$total_production[indices]
data.tab$coal.mines <- mining.tab$total_mines[indices]
data.tab$coal.surface[is.na(data.tab$coal.surface)] <- 0
data.tab$coal.underground[is.na(data.tab$coal.underground)] <- 0
data.tab$coal.mines[is.na(data.tab$coal.mines)] <- 0
data.tab$coal.total[is.na(data.tab$coal.total)] <- 0
}

# Add total cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'all', 'all.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'all_cancer'] <- tab[indices, 'incidence']

# Add total cancer incidence male
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'all', 'male.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'all_cancer_male'] <- tab[indices, 'incidence']

# Add total cancer incidence female
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'all', 'female.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'all_cancer_female'] <- tab[indices, 'incidence']


# Add colorectal cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'colorectal', 'all.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab$colorectal <- tab$incidence[indices]
data.tab$colorectal_lower <- tab[indices, 'lower.95']
data.tab$colorectal_upper <- tab[indices, 'upper.95']


# Add breast cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'breast', 'female.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab$breast <- tab$incidence[indices]
data.tab$breast_lower <- tab[indices, 'lower.95']
data.tab$breast_upper <- tab[indices, 'upper.95']


# Add breast cancer incidence for 50 years or older
if (extras) {
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'breast', 'female-50-plus.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'over_50_breast'] <- tab$incidence[indices]
}

# Add prostate cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'prostate', 'male.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab$prostate <- tab$incidence[indices]
data.tab$prostate_lower <- tab[indices, 'lower.95']
data.tab$prostate_upper <- tab[indices, 'upper.95']


# Add melanoma cancer incidence
if (extras) {
path <- file.path(cancer.dir, 'incidence', '2006-2010', 'melanoma', 'all.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab$melanoma <- tab$incidence[indices]
data.tab$melanoma_lower <- tab[indices, 'lower.95']
data.tab$melanoma_upper <- tab[indices, 'upper.95']
}

# Add lung cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'lung', 'all.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'lung'] <- tab[indices, 'incidence']
data.tab$lung_lower <- tab[indices, 'lower.95']
data.tab$lung_upper <- tab[indices, 'upper.95']


# Add female lung cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'lung', 'female.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'female_lung'] <- tab[indices, 'incidence']


# Add male lung cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'lung', 'male.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'male_lung'] <- tab[indices, 'incidence']


# Add over65 lung cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'lung', 'over65.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'over_65_lung'] <- tab[indices, 'incidence']


# Add under65 lung cancer incidence
path <- file.path(cancer.dir, 'incidence', cancer.time.period, 'lung', 'under65.csv')
tab <- ReadStateCancerProfileIncidence(path)
indices <- match(data.tab$fips, tab$fips)
data.tab[, 'under_65_lung'] <- tab[indices, 'incidence']


# Write to file

#filename <- sprintf('data-%s%s.txt', cancer.time.period, ifelse(extras, '-extras', ''))
#path <- file.path(dropbox.dir, 'data', 'county', filename)
filename <- sprintf('county-data%s.txt', ifelse(extras, '-extras', ''))
path <- file.path(dropbox.dir, 'analysis', 'data', filename)
write.table(data.tab, path, row.names=FALSE, sep='\t', quote=FALSE, na='')
print(warnings())
}
