# Use by writing state fips codes after rscript file.
# example usage:
# Rscript elevation-blockgroup-averaged.R 08 06 04 16 32 35 41 49 53 30
args <- commandArgs(trailingOnly = TRUE)

lung.data.dir <- '/home/dhimmels/Documents/lung-data'
setwd(lung.data.dir)

# Read block group shapefiles
shape.dir <- file.path(lung.data.dir, 'census', '2000', 'blockgroup', 'shapefiles')

shp.info <- lapply(list.dirs(shape.dir, recursive=FALSE), function(directory) {
  file.name <- list.files(directory, pattern="[.]shp$")[1]
  file.name <- sub("[.]shp$", '', file.name)
  state.fips <- substr(file.name, 3, 4) # for 2000 blockgroup filenames
  c(state.fips, directory, file.name)
})
shp.info <- do.call(rbind, shp.info)
shp.info <- as.data.frame(shp.info, stringsAsFactors=FALSE)
colnames(shp.info) <- c('state.fips', 'directory', 'file.name')

# subset shp.info based on args
shp.info <- shp.info[match(args, shp.info$state.fips), ]

library('raster')
library('rgdal')
ReadStateShapefile <- function(fips) {
  shp.row <- shp.info[shp.info$state.fips == fips, ]
  readOGR(shp.row$directory, shp.row$file.name, stringsAsFactors=FALSE)
}

# Read elevation raster files
#elevation.dir <- file.path('raster', 'elevation', 'usgs', '33')
elevation.dir <- file.path('raster', 'elevation', 'worldclim')
elevation.files <- list.files(elevation.dir)

raster.info <- lapply(elevation.files, function(elevation.file) {
  rast <- raster(file.path(elevation.dir, elevation.file))
  raster.extent <- extent(rast)
  extent.vector <- as.numeric(c(xmin(raster.extent), xmax(raster.extent),
    ymin(raster.extent), ymax(raster.extent)))
  c(elevation.file, extent.vector)
})
raster.info <- do.call(rbind, raster.info)
raster.info <- data.frame(raster.info, stringsAsFactors=FALSE)
for (j in 2:4) raster.info[, j] <- as.numeric(raster.info[, j])
colnames(raster.info) <- c('file.name', 'xmin', 'xmax', 'ymin', 'ymax')

# Reversed comparison order.
OverlappingRasters <- function(ext, raster.info) {
  overlap.bool <- apply(raster.info, 1, function(raster.row) {
    suppressWarnings(intersection <- intersect(ext, extent(as.numeric(raster.row[-1]))))
    return(!is.null(intersection))
  })
  file.names <- raster.info[overlap.bool, 'file.name']
  return(file.names)
}


cat('Computing population averaged elevation\n')
#shp.info <- shp.info[shp.info$state.fips != '02', ] # remove alaska
for (state.fips in shp.info$state.fips) {
  shapes <- ReadStateShapefile(state.fips)
  # Save a plot of the blockgroups for the state.
  file.name <- paste(state.fips, '-blockgroup-plot.eps', sep='')
  path <- file.path('processed', 'elevation', 'blockgroup', file.name)
  setEPS()
  pdf(path, width=5, height=5)
  plot(shapes, main=paste('blockgroups for fips ', state.fips, sep=''))
  dev.off()
  #shapes <- shapes[1:5, ] # TESTING ONLY REMOVE
  cat(paste(nrow(shapes), ' blockgroups for fips ', state.fips, '\n', sep=''))
  # Patch together raster elevation files if necessary
  raster.files <- OverlappingRasters(extent(shapes), raster.info)
  rasters <- sapply(raster.files, function(file.name)
    raster(file.path(elevation.dir, file.name)))
  merged.raster <- rasters[[1]]
  for (rast in rasters[-1]) {
    cat('merging raster\n')
    merged.raster <- merge(merged.raster, rast)
  }
  cat('rasters have been merged\n')
  cat('computing mean elevation for each blockgroup\n')
  means <- extract(merged.raster, shapes, weights=TRUE, fun=weighted.mean)
  cat(paste(sum(is.na(means)), ' blockgroups with missing altitude when using weights=TRUE\n', sep=''))
  mean.missing = is.na(means)
  if (any(mean.missing)) {
    means[mean.missing] <- raster::extract(merged.raster, shapes[mean.missing, ], small=TRUE, fun=mean)  
  }
  cat(paste(sum(mean.missing), ' blockgroups with missing altitude after filling in with small=TRUE\n', sep=''))
  means <- round(means, 3)
  #county <- paste(shapes$STATE, shapes$COUNTY, sep='')
  # Fix issue where some tracks have last two zero's omitted.
  tract <- shapes$TRACT
  tract[nchar(tract) == 4] <- paste(tract[nchar(tract) == 4], '00', sep='')
  # Check for other clipped fips
  stopifnot(all(nchar(shapes$STATE) == 2))
  stopifnot(all(nchar(shapes$COUNTY) == 3))
  stopifnot(all(nchar(tract) == 6))
  stopifnot(all(nchar(shapes$BLKGROUP) == 1))

  blockgroup <- paste(shapes$STATE, shapes$COUNTY, tract, shapes$BLKGROUP, sep='')
  shape.data <- data.frame(blockgroup, 'area'=shapes$AREA, 'mean.elev'=means)
  
  shape.data <- shape.data[order(shape.data$blockgroup), ]
  file.name <- paste(state.fips, '-blockgroup-elevation.txt', sep='')
  path <- file.path('processed', 'elevation', 'blockgroup', file.name)
  write.table(shape.data, path, sep='\t', row.names=FALSE, quote=FALSE)
  cat(paste('blockgroup mean elevations calculated for fips ', state.fips, '\n', sep=''))
  cat('---------------------------------------------------------------------\n')
}

