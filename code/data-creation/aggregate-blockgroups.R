
data.dir <- '/home/dhimmels/Documents/lung-data'
dropbox.dir <- '~/Dropbox/lung/'

pop.path <- file.path(data.dir, 'census', '2000', 'blockgroup', 'blkgrp_pop_centroid_withname.txt')
pop.data <- read.csv(pop.path, colClasses='character')

pop.tab <- data.frame('blockgroup'=paste(pop.data$state, pop.data$county, pop.data$tract, pop.data$blkgrp, sep=''),
                      'population'=as.numeric(pop.data$pop))

blockgroup.dir <- file.path(data.dir, 'processed', 'elevation', 'blockgroup')
blockgroup.paths <- list.files(blockgroup.dir, pattern="[.]txt$", full.names=TRUE)

blockgroup.tabs <- lapply(blockgroup.paths, read.delim, colClasses='character')
blockgroup.tab <- do.call(rbind, blockgroup.tabs)
blockgroup.tab$county <- substr(blockgroup.tab$blockgroup, 1, 5)
blockgroup.tab$state <- substr(blockgroup.tab$blockgroup, 1, 2)
blockgroup.tab$area <- as.numeric(blockgroup.tab$area)
blockgroup.tab$elevation <- as.numeric(blockgroup.tab$mean.elev)
blockgroup.tab$population <- pop.tab[match(blockgroup.tab$blockgroup, pop.tab$blockgroup), 'population']

path <- file.path(data.dir, 'processed', 'elevation', 'county', 'blockgroups-2000.txt')
write.table(blockgroup.tab, path, row.names=FALSE, quote=FALSE, sep='\t', na='')
new.path <- file.path(dropbox.dir, 'data', 'blockgroups-2000.txt')
invisible(file.copy(path, new.path, overwrite = TRUE))




AggregateWeightedMean <- function(tab, grouping.var, average.vars, weighting.var) {
  # Like aggregate function but allows weighted mean.
  # Also returns the sum of the weighting variable
  groups <- levels(as.factor(tab[, grouping.var]))
  weighted.mean.list <- lapply(groups, function(group) {
    keep.row <- tab[, grouping.var] == group
    tab.subset <- tab[keep.row, average.vars]
    if (length(average.vars) == 1) {
      tab.subset <- data.frame(tab.subset, stringsAsFactors=FALSE)
      colnames(tab.subset) <- average.vars
    }
    subset.weights <- tab[keep.row, weighting.var]
    weighted.means <- apply(tab.subset, 2, weighted.mean, w=subset.weights, na.rm=TRUE)
    weighted.means <- c(weighted.means, sum(subset.weights))
    names(weighted.means) <- c(average.vars, weighting.var)
    return(weighted.means)
  })
  aggregate.tab <- do.call(rbind, weighted.mean.list)
  aggregate.tab <- data.frame(groups, aggregate.tab, stringsAsFactors=FALSE)
  colnames(aggregate.tab)[1] <- grouping.var
  return(aggregate.tab)
}




county.elev.tab <- AggregateWeightedMean(blockgroup.tab, 'county', c('elevation'), 'population')
county.elev.tab$elevation <- round(county.elev.tab$elevation, 2)
path <- file.path(data.dir, 'processed', 'elevation', 'county', 'blockgroup-averaged-2000.txt')
write.table(county.elev.tab, path, row.names=FALSE, quote=FALSE, sep='\t')
new.path <- file.path(dropbox.dir, 'data', 'county', 'environment', 'elevation-blockgroup-averaged-2000.txt')
invisible(file.copy(path, new.path, overwrite = TRUE))


state.elev.tab <- AggregateWeightedMean(blockgroup.tab, 'state', c('elevation'), 'population')
state.elev.tab$elevation <- round(state.elev.tab$elevation, 2)
path <- file.path(data.dir, 'processed', 'elevation', 'state', 'blockgroup-averaged-2000.txt')
write.table(state.elev.tab, path, row.names=FALSE, quote=FALSE, sep='\t')
new.path <- file.path(dropbox.dir, 'data', 'state', 'environment', 'elevation-blockgroup-averaged-2000.txt')
invisible(file.copy(path, new.path, overwrite = TRUE))

