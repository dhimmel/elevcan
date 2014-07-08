# daniel.himmelstein@gmail.com
# kamen.simeonov@gmail.com


################################################################################
############################ Parameter Specification ###########################

## File locations
project.dir <- getwd()
code.dir <- file.path(project.dir, 'code')
data.dir <- file.path(project.dir, 'data')
figure.dir <- file.path(project.dir, 'figures')
font.dir <- file.path(project.dir, 'fonts')
table.dir <- file.path(project.dir, 'tables')
output.dir <- file.path(project.dir, 'output')
figdata.dir <- file.path(output.dir, 'figdata')
county.data.path <- file.path(data.dir, 'county-data.txt')

## County filtering parameters
states <- c('AZ', 'CA', 'CO', 'ID', 'MT', 'NV', 'NM', 'OR', 'UT', 'WA', 'WY')
min.population=1e4
max.native=25.0
max.migration=40.0

## Cancers to evaluate (sets plotting order)
cancers <- c('lung', 'breast', 'colorectal', 'prostate')


global.covars <- c('metro', 'white', 'black', 'education', 'income', 'obesity')
specific.covars <- list(
  'lung'       = c('no_lung', 'smoking', 'radon', 'particulate'),
  'breast'     = c('no_breast', 'female_smoking', 'mammogram', 'drinking'),
  'colorectal' = c('no_colorectal', 'smoking', 'drinking', 'diabetes', 'meat'),
  'prostate'   = c('no_prostate', 'meat')
)
all.covars <- sort(unique(c(global.covars, do.call(c, specific.covars))))
envir.covars <- c('elevation', 'uvb', 'sunlight', 'precipitation', 'high_temp', 'diurnal_temp', 'radon', 'particulate')


################################################################################
################################# Execution ####################################

# Read helper functions
source(file.path(code.dir, 'functions.R'))

# Enable logging
sink(file=file.path(output.dir, 'log.txt'))
time.chr <- strftime(Sys.time(), format="%y-%m-%d %H:%M:%S %Z")
cat(sprintf('Initiating Analysis: %s\n', time.chr)); CatDiv()

# Read county data
data.df <- ReadCountyData(county.data.path, 
  states=states, min.population=min.population, 
  max.native=max.native, max.migration=max.migration, verbose=TRUE)

# Calculate models
source(file.path(code.dir, 'create-models.R'))

# Draw Figures
source(file.path(code.dir, 'create-figures.R'))

# Create Tables
source(file.path(code.dir, 'create-tables.R'))

# Save R session
save.image(file=file.path(output.dir, 'session.RData'))

