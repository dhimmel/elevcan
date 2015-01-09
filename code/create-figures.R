# daniel.himmelstein@gmail.com

################################################################################
################################ Package Loading ###############################
library(ggplot2)
library(grid)
library(metafor)
AssertPackages('reshape2')

#options(stringsAsFactors=FALSE)
################################################################################
############################ Parameter Specification ###########################

## PLOS Figure Widths (inches)
full.width  <- 6.83
half.width  <- 3.27
full.height <- 9.19

## Color options
biv.line.col   <- '#990000'
biv.fill.col   <- '#ffb2b2'
par.line.col <- '#000099'
par.fill.col <- '#b2b2ff'

strip.fill <- '#fef2e2'

div.cols <- c('#0000FF', '#f7eff7', '#FF0000')
seq.cols <- c('#ccccff', '#b2b2ff', '#7f7fff', '#6666ff', '#000099', '#00007f')

fill.red <- '#FB7C72'
fill.blue <- '#7284FF'

## Alpha options for point transperency
alpha.range <- c(0.1, 1)
weight.to.alpha.exp <- 1

################################################################################
############################# Function Declaration #############################

OpenPDF <- function(filename, width=4.86, height=4.8) {
  path <- file.path(figure.dir, filename)
  cairo_pdf(path, width=width, height=height, bg='white', family='sans')
}

ClosePDF <- function(filename) {
  dev.off()
  path <- file.path(figure.dir, filename)  
  embedFonts(path)
}

LabelR2 <- function(model) {
  r.squared <- summary(model)$r.squared
  r.squared <- ChrRound(r.squared, 2)
  r.squared <- sprintf("R^2=='%s'", r.squared)
  return(r.squared)
}

LabelCoef <- function(model, coef.var='elevation') {
  slope <- as.numeric(coef(model)[coef.var])
  slope <- ChrRound(slope, 2)
  slope <- sprintf("beta=='%s'", slope)
  return(slope)
}

LabelStandardCoef <- function(model, coef.var='elevation') {
  model.df <- model.frame(model)
  standardizing.denom <- (Hmisc::wtd.var(model.df[, 1], weights=model.df[, '(weights)'], normwt=TRUE) ^ 0.5 /
    Hmisc::wtd.var(model.df[, coef.var], weights=model.df[, '(weights)'], normwt=TRUE) ^ 0.5)
  standard.beta <- as.numeric(coef(model)[coef.var] / standardizing.denom)
  slope <- ChrRound(standard.beta, 2)
  slope <- sprintf("beta[z]=='%s'", slope)
  return(slope)
}

SetGGTheme <- function(gg) {
  gg <- gg + theme_bw()
  gg <- gg + theme(plot.margin=grid::unit(c(2, 2, 2, 2), 'points'))
  return(gg)
}


################################################################################
################################# Execution ####################################

######################################
## Filtering threshold determination 

unfiltered.df <- ReadCountyData(county.data.path, 
  states=states, min.population=min.population, 
  max.native=100, max.immigration=100)

outlier.covars <- c(global.covars, 'smoking', 'male')
outlier.form <- as.formula(paste('lung ~', paste(c('1', outlier.covars), collapse=' + ')))
outlier.model <- lm(outlier.form, weight=weight, unfiltered.df, na.action='na.exclude')

filter.df <- data.frame(
  'fips'=unfiltered.df[, 'fips'],
  'name'=unfiltered.df[, 'name'],
  'weight'=unfiltered.df[, 'weight'],
  'residual'=resid(outlier.model),
  'Immigration'=unfiltered.df[, 'immigration'],
  'Native'=unfiltered.df[, 'native']
)

filter.melt <- reshape2::melt(filter.df, 
  measure.vars=c('Native', 'Immigration'), value.name='percent')
threshold.df <- data.frame('variable'=c('Native', 'Immigration'), 
  'xmin'=c(max.native, max.immigration), xmax=Inf, ymin=-Inf, ymax=Inf)

gg.filter <- ggplot(filter.melt, aes(percent, abs(residual)))
gg.filter <- SetGGTheme(gg.filter) +
  facet_grid(. ~ variable, scales='free_x') +
  geom_rect(data=threshold.df, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, x=NULL, y=NULL),
    fill=biv.fill.col, alpha=0.5) +
  geom_smooth(aes(weight=weight), span=5, color=par.line.col, fill=par.fill.col, method='loess') +
  geom_point(aes(alpha=weight^.5)) + scale_alpha(range=c(0.1, 1)) +
  guides(alpha=FALSE) +
  xlab('Percent') + ylab('Absolute Residual') +
  theme(strip.background=element_rect(fill=strip.fill))

OpenPDF('filter.pdf', width=full.width, height=3)
print(gg.filter); ClosePDF('filter.pdf')


############################
## Variable Correlation Plot

HclustOrder <- function(variable.df) {
  # Returns the colnames, ordered by heirarchical clustering
  cor.mat <- cor(variable.df, use='pairwise.complete.obs')
  clust <- hclust(dist(cor.mat), method='ward')
  clust.order <- colnames(variable.df)[clust$order]
  return(clust.order)
}

SelfCorrelationDF <- function(variable.df) {
  # Returns a data.frame where each row shows the correlation between two variables.
  cor.tri.mat <- cor(variable.df, use='pairwise.complete.obs')
  cor.tri.mat[lower.tri(cor.tri.mat, diag=TRUE)] <- NA
  cor.melt <- reshape2::melt(cor.tri.mat, varnames=c('x', 'y'), value.name='correlation')
  cor.melt <- na.omit(as.data.frame(cor.melt))
}

OutcomePredictorCorPlot <- function(variable.df, outcomes, predictors, 
  sort.predictors=TRUE, barwidth=10, barheight=1, xtext.angle=65,
  legend.position =c(1, 0.75)) {

  # Sort the predictors using heirarchical clustering
  if (sort.predictors) {
    predictors <- HclustOrder(variable.df[, predictors])
  }

  # Compute pairwise correlations
  predictor.cor.df <- SelfCorrelationDF(variable.df[, predictors])
  outcome.cor.df <- reshape2::melt(cor(variable.df[, outcomes, drop=FALSE], 
    variable.df[, predictors], use='pairwise.complete.obs'),
    varnames=c('x', 'y'), value.name='correlation')
  cor.df <- rbind(outcome.cor.df, predictor.cor.df)

  # Plot
  xlimits <- c(outcomes, predictors[-length(predictors)])
  ylimits <- rev(predictors)
  xlabels <- gsub('_', ' ', xlimits)
  ylabels <- gsub('_', ' ', ylimits)
  gg.cor <- ggplot(cor.df, aes(x, y, fill=correlation))
  gg.cor <- SetGGTheme(gg.cor) +
    geom_tile(color='white') + 
    scale_fill_gradientn(name=expression("Pearson's" * ~ rho),
      limits=c(-1, 1), colours=div.cols) + 
    scale_x_discrete(limits=xlimits, labels=xlabels) + 
    scale_y_discrete(limits=ylimits, labels=ylabels) +
    theme(axis.text.x=element_text(angle=xtext.angle, hjust=1)) +
    xlab(NULL) + ylab(NULL) + coord_fixed() +
    theme(plot.background=element_blank(), panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), panel.border=element_blank(),
      axis.line=element_line(color='grey50', size=0.25)) +
    theme(legend.justification=c(1, 0),
      legend.position=legend.position, 
      legend.direction='horizontal') +
    guides(fill=guide_colorbar(barwidth=barwidth, barheight=barheight, title.position='top', title.hjust=0.5)) +
    geom_vline(xintercept=length(outcomes) + 0.5, size=0.25, color='grey50')

  # Return a list with the correlations and ggplot
  cor.list <- list()
  cor.list$data <- cor.df
  cor.list$plot <- gg.cor
  return(cor.list)
}

other.cancers <- c('no_lung', 'no_breast', 'no_prostate', 'no_colorectal')
clust.predictors <- setdiff(c('elevation', 'all_cancer', all.covars), other.cancers)
clust.cancers <- c('lung', 'colorectal', 'breast', 'prostate')
cor.list <- OutcomePredictorCorPlot(data.df, clust.cancers, clust.predictors)

OpenPDF('correlation.pdf', width=4.86, height=4.25)
print(cor.list$plot); ClosePDF('correlation.pdf')

write.table(cor.list$data, file.path(figdata.dir, 'variable-correlation.txt'), 
  row.names=FALSE, sep='\t', quote=FALSE)
#####################################
## Bivariate and Partial Scatterplots

# Create data frames for ggplot
biv.df <- do.call(rbind, lapply(cancers, function(cancer) {
  cancer.df <- cancer.list[[cancer]]$cancer.allvar.df
  data.frame(
    'fips'      = cancer.df[, 'fips'],
    'Elevation' = cancer.df[, 'elevation'],
    'Incidence' = cancer.df[, cancer],
    'weight'    = cancer.df[, 'weight'],
    'Cancer'    = SimpleCap(cancer))
}))

par.df <- do.call(rbind, lapply(cancers, function(cancer) {
  par.cancer.df <- par.df.list[[cancer]]
  data.frame(
    'fips'      = par.cancer.df[, 'fips'],
    'Elevation' = par.cancer.df[, 'elevation_residual'],
    'Incidence' = par.cancer.df[, 'cancer_residual'],
    'weight'    = par.cancer.df[, 'weight'],
    'Cancer'    = SimpleCap(cancer))
}))

# Create data frames with geom_text labels
biv.label.df <- do.call(rbind, lapply(cancers, function(cancer) {
  biv.model <- biv.model.list[[cancer]]
  data.frame('Elevation'=Inf, 'Incidence'=Inf, 'Cancer'=SimpleCap(cancer), 
    'coef'=LabelCoef(biv.model), 'r2'=LabelR2(biv.model))
}))

par.label.df <- do.call(rbind, lapply(cancers, function(cancer) {
  par.model <- par.model.list[[cancer]]
  data.frame('Elevation'=Inf, 'Incidence'=Inf, 'Cancer'=SimpleCap(cancer), 
    'coef'=LabelCoef(par.model, 'elevation_residual'), 'r2'=LabelR2(par.model))
}))


GGScatter <- function(scat.df, line.col, fill.col, x.lab, y.lab, label.df) {
  gg <- ggplot(scat.df, aes(Elevation, Incidence))
  gg <- SetGGTheme(gg)
  gg <- gg + geom_smooth(method='lm', aes(weight=weight), level=0.99,
    color=line.col, fill=fill.col, alpha=1) +
  geom_point(aes(alpha=weight ^ weight.to.alpha.exp)) +
  facet_wrap(~ Cancer, scales='free', ncol=1) +
  scale_alpha(range=alpha.range) + guides(alpha=FALSE) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  xlab(x.lab) + ylab(y.lab) +
  geom_text(aes(label=coef), data=label.df, hjust=1.1, vjust=1.2, parse=TRUE, size=3.5) +
  geom_text(aes(label=r2), data=label.df, hjust=1.1, vjust=2.2, parse=TRUE, size=3.5)
  return(gg)
}

ggbiv <- GGScatter(biv.df, line.col=biv.line.col, fill.col=biv.fill.col,
  x.lab='Elevation (km)', y.lab='Incidence (Age-adjusted cases per 100,000)',
  label.df=biv.label.df)

ggpar <- GGScatter(par.df, line.col=par.line.col, fill.col=par.fill.col,
  x.lab='Elevation Residual', y.lab='Incidence Residual',
  label.df=par.label.df)

OpenPDF('scatterplots.pdf', width=full.width, height=12)
gridExtra::grid.arrange(ggbiv, ggpar, nrow=1)
ClosePDF('scatterplots.pdf')

write.table(biv.df, file.path(figdata.dir, 'scatterplot-bivariate.txt'), 
  row.names=FALSE, sep='\t', quote=FALSE)
write.table(par.df, file.path(figdata.dir, 'scatterplot-partial.txt'), 
  row.names=FALSE, sep='\t', quote=FALSE)

#####################################
## Best Subset Plot

BestSubsetDF <- function(cancer) {
  rsubs <- cancer.list[[cancer]]$rsubs
  max.size <- rsubs$nvmax - 2
  predictors.list <- append(list('elevation'), 
    lapply(1:max.size, function(i) names(coef(rsubs, i)[-1])))
  subset.df <- do.call(rbind, lapply(predictors.list, function(predictors) {
    form <- as.formula(sprintf('%s ~ 1 + %s', cancer, paste(predictors, collapse=' + ')))
    cancer.zdf <- cancer.list[[cancer]]$cancer.zdf
    model <- lm(form, weights=weight, data=cancer.zdf)
    conf.int <- confint(model, level=0.99)['elevation', ]
    size <- length(predictors)
    subset.row <- data.frame(
      'Cancer'=SimpleCap(cancer), 'size'=size, 'bic'=BIC(model),
      'estimate'=coef(model)['elevation'], 'lower'=conf.int[1], 'upper'=conf.int[2],
      'optimal'=length(cancer.list[[cancer]]$predictors) == size)
  }))
  return(subset.df)
}

subset.df <- do.call(rbind, lapply(cancers, BestSubsetDF))
row.names(subset.df) <- NULL
subset.df <- subset(subset.df, size <= 10)


gg.bs <- ggplot(subset.df, aes(x=size, y=estimate, ymin=lower, ymax=upper, color=bic))
gg.bs <- SetGGTheme(gg.bs) +
  geom_hline(yintercept=0, linetype='dashed', color='darkgray') +
  geom_linerange(size=1.1) + facet_grid(~ Cancer) + 
  geom_point(size=2.5, aes(shape=optimal)) +
  scale_shape_manual(values=c(16, 11), guide='none') +
  scale_colour_gradientn(colours=rev(seq.cols)) +
  guides(color=guide_colorbar(title='Model BIC')) + 
  scale_x_discrete(breaks=c(2, 5, 8)) +
  scale_y_continuous(breaks=seq(-1, 1, 0.2)) +
  theme_bw() + xlab('Subset Size') + ylab('Standardized Elevation Coefficient') +
  theme(strip.background=element_rect(fill=strip.fill)) +
  theme(plot.margin=grid::unit(c(2, 0, 2, 2), 'points')) +
  theme(legend.margin=grid::unit(0, 'cm')) +
  theme(legend.title=element_text(face='plain'))

OpenPDF('best-subsets.pdf', width=full.width, height=3.4)
print(gg.bs); ClosePDF('best-subsets.pdf')

write.table(subset.df, file.path(figdata.dir, 'best-subsets.txt'), 
  row.names=FALSE, sep='\t', quote=FALSE)


#########################################
## Lasso Viz

ExtractLassoDF <- function(cancer) {
  lasso <- cancer.list[[cancer]]$lasso 
  coefs <- coef(lasso, s=lasso$lambda.1se)[, 1]
  coefs <- coefs[coefs != 0]
  coefs <- coefs[names(coefs) != '(Intercept)']
  coefs <- coefs[order(abs(coefs), decreasing=TRUE)]
  lasso.df <- data.frame(
    'predictor'=factor(names(coefs), levels=names(coefs)),
    'coefficient'=as.numeric(coefs), 'positive'=sign(coefs) > 0, 'Cancer'=SimpleCap(cancer))
  return(lasso.df)
}


gg.lasso <- ggplot(mapping=aes(predictor, abs(coefficient), fill=positive))
gg.lasso <- SetGGTheme(gg.lasso)
max.coefs <- c()
for (cancer in cancers) {
  lasso.df <- ExtractLassoDF(cancer)
  max.coefs <- append(max.coefs, max(abs(lasso.df$coefficient)))
  gg.lasso <- gg.lasso + geom_bar(data=lasso.df, stat='identity', color='black') +
  geom_text(data=lasso.df, aes(label=ChrRound(coefficient, 2)), vjust=-0.3, size=3.5)
}
gg.lasso <- gg.lasso + facet_grid(. ~ Cancer, space='free_x', scales='free_x') +
  theme(axis.text.x=element_text(angle=38, hjust=1)) +
  theme(strip.background=element_rect(fill=strip.fill)) +
  scale_fill_manual(values=c(par.fill.col, biv.fill.col)) +
  xlab('Predictor') + ylab('Absolute Standardized Coefficient') +
  guides(fill=FALSE) + ylim(c(0, max(max.coefs) + 0.025))

OpenPDF('lasso.pdf', width=full.width, height=4)
print(gg.lasso); ClosePDF('lasso.pdf')

lasso.df <- do.call(rbind, lapply(cancers, ExtractLassoDF))
write.table(lasso.df, file.path(figdata.dir, 'predictors-lasso.txt'), row.names=FALSE, sep='\t', quote=FALSE)



#####################################
## Stratification

# smoking
lung.df$smoking.qauntile <- cut(lung.df$smoking, breaks=quantile(lung.df$smoking, probs=c(0, 1/3, 2/3, 1)), include.lowest=TRUE)
terciles <- rev(levels(lung.df$smoking.qauntile))
cat(AddNewLines(sprintf('Smoking terciles for lung cancer: %s, %s, %s', terciles[1], terciles[2], terciles[3]))); CatDiv()
terc.labels <- c('high', 'mid', 'low')
gg.smoke <- ggplot(lung.df, aes(elevation, lung, linetype=smoking.qauntile, shape=smoking.qauntile))
gg.smoke <- SetGGTheme(gg.smoke) + 
  geom_smooth(method='lm', aes(weight=weight),
    color=NA, fill=biv.fill.col, alpha=1, show_guide=FALSE) +
  geom_smooth(method='lm', aes(weight=weight), color=biv.line.col, fill=NA) +
  geom_point(aes(alpha=weight^.5), fill='black', color=NA) + 
  scale_alpha(range=c(0.25, 1)) +
  scale_linetype_manual(values=c('dotted', 'dashed', 'solid'), breaks=terciles, labels=terc.labels) +
  scale_shape_manual(values=c(25, 21, 24), breaks=terciles, labels=terc.labels) +
  theme(legend.justification=c(1,1), legend.position=c(1, 1),
    legend.background=element_rect(color='#CCCCCC', size=0.2), 
    legend.key.width=grid::unit(2, 'lines'),
    legend.key.height=grid::unit(0.95, 'lines'),
    legend.key=element_rect(linetype='blank')) +
  guides(alpha=FALSE, shape=guide_legend('Smoking\nPrevalence', title.hjust=0.5), 
    linetype=guide_legend('Smoking\nPrevalence', title.hjust=0.5)) +
  theme(legend.title=element_text(face='plain')) +
  xlab('Elevation (km)') + ylab('Lung Cancer Incidence')

# State
lung.allvar.df <- cancer.list$lung$cancer.allvar.df

state.df <- do.call(rbind, lapply(states, function(state) {
  state.df <- lung.allvar.df[lung.allvar.df$state == state, ]
  state.lm <- lm(lung ~ elevation + smoking, weight=weight, data=state.df)
  conf.int <- confint(state.lm, level = 0.95)
  data.frame('state'=state, 
   'state_n'=sprintf('%s (%s)', state, nrow(state.df)),
   'elevation.coef'=coef(state.lm)['elevation'], 
   'elevation.se'=summary(state.lm)$coef['elevation', 'Std. Error'],
   'elevation.lower'=conf.int['elevation', 1], 'elevation.upper'=conf.int['elevation', 2],
   'smoking.coef'=coef(state.lm)['smoking'],
   'smoking.lower'=conf.int['smoking', 1], 'smoking.upper'=conf.int['smoking', 2],
   stringsAsFactors=FALSE)
}))
row.names(state.df) <- state.df$state
state.x.limits <- c('Meta', state.df[order(state.df$elevation.coef), 'state_n'])

# Computation of the I^2
state.random <- metafor::rma.uni(yi=state.df$elevation.coef,
  sei=state.df$elevation.se, level=99,  method='REML')
rma.confint <- metafor::confint.rma.uni(
  state.random, level=95, control=list('tau2.max'=1e5))
i2 <- rma.confint$random['I^2(%)', 'estimate']
i2.lower <- rma.confint$random['I^2(%)', 'ci.lb']
i2.upper <- rma.confint$random['I^2(%)', 'ci.ub']
cat(sprintf('Statewise Meta-Analysis I-squared:\n%.1f%% [95%% CI: %.1f%%, %.1f%%]\n',
  i2, i2.lower, i2.upper))
CatDiv()

# Fixed Effects meta-analysis
state.fixed <- metafor::rma.uni(yi=state.df$elevation.coef,
  sei=state.df$elevation.se, level=99,  method='FE')

state.df <- rbind(state.df, data.frame(
  'state'='Meta', 'state_n'='Meta', 'elevation.coef'=state.fixed$b[1], 
  'elevation.se'=NA, 'elevation.lower'=state.fixed$ci.lb, 'elevation.upper'=state.fixed$ci.ub,
  'smoking.coef'=NA, 'smoking.lower'=NA, 'smoking.upper'=NA))
lung.conf.int <- confint(model.list$lung, level=0.99)['elevation', ]

gg.state <- ggplot(state.df, aes(state_n, elevation.coef, ymin=elevation.lower, ymax=elevation.upper))
gg.state <- SetGGTheme(gg.state) + 
  geom_hline(yintercept=0, linetype='dashed', color='darkgray') + 
  geom_rect(xmin=-Inf, xmax=Inf, fill=par.fill.col,
    ymin=lung.conf.int[1], ymax=lung.conf.int[2]) +
  geom_errorbar(aes(color=state_n == 'Meta'), size=.4, width=.2) +
  geom_point(aes(color=state_n == 'Meta'), size=2, shape=19) + 
  xlab('State') + ylab('Elevation Coefficient') + ylab(expression(beta[elevation])) + 
  scale_x_discrete(limits=state.x.limits) + 
  scale_color_manual(values=c('#444444', 'black'), guide=FALSE) +
  coord_flip() + theme(axis.text.y=element_text(angle=50, hjust=1))

write.table(state.df, file.path(figdata.dir, 'lung-state-strat.txt'), row.names=FALSE, sep='\t', quote=FALSE)

OpenPDF('stratification.pdf', width=full.width, height=4)
gridExtra::grid.arrange(gg.smoke, gg.state, nrow=1, widths=c(2, 1.1))
ClosePDF('stratification.pdf')

################################################################################
# Subgrouping
options(stringsAsFactors=FALSE)
strip.labeller <- function(variable, value) {
  if(variable == 'group') {
    value[value == 'male_lung'] <- 'Male'
    value[value == 'female_lung'] <- 'Female'
    value[value == 'over_65_lung'] <- 'Over 65'
    value[value == 'under_65_lung'] <- 'Under 65'
    return(value)
  }
  if(variable == 'group_type') {
    return(value)
  }
}

SubgroupPlot <- function(group.df, label.df, margin.lines) {
  gg.group <- ggplot(group.df, aes(elevation_residual, incidence_residual)) +
  facet_grid(group ~ group_type, scales='free_y', labeller=strip.labeller) +
  theme_bw() + theme(strip.background=element_rect(fill=strip.fill)) +
  geom_smooth(method='lm', aes(weight=weight),
    color=par.line.col, fill=par.fill.col, alpha=1) +
  geom_point(aes(alpha=weight^.5)) + scale_alpha(range=c(0.1, 1)) +
  xlab(NULL) + ylab(NULL) + guides(alpha=FALSE) +
  theme(plot.margin = grid::unit(margin.lines, 'lines')) +
  geom_text(aes(label=coef), data=label.df, hjust=1.1, vjust=1.2, parse=TRUE, size=3.5) +
  geom_text(aes(label=zcoef), data=label.df, hjust=1.1, vjust=2.3, parse=TRUE, size=3.5) +
  geom_text(aes(label=r2), data=label.df, hjust=1.1, vjust=3.2, parse=TRUE, size=3.5)
  return(gg.group)
}

SubgroupLabelDF <- function(group, group_type) {
  model <- group.par.list[[group]]
  coef.label <- LabelCoef(model, 'elevation_residual')
  zcoef.label <- LabelStandardCoef(model, 'elevation_residual')
  r2.label <- LabelR2(model)
  data.frame('elevation_residual'=Inf, 'incidence_residual'=Inf, 'group'=group,
   'group_type'=group_type, 'coef'=coef.label, 'zcoef'=zcoef.label, 'r2'=r2.label)
}

age.label.df <- rbind(SubgroupLabelDF('over_65_lung', 'Age'), 
  SubgroupLabelDF('under_65_lung', 'Age'))
sex.label.df <- rbind(SubgroupLabelDF('male_lung', 'Sex'), 
  SubgroupLabelDF('female_lung', 'Sex'))


age.df <- rbind(group.df.list$over_65_lung, group.df.list$under_65_lung)
age.df$group_type <- 'Age'
sex.df <- rbind(group.df.list$male_lung, group.df.list$female_lung)
sex.df$group_type <- 'Sex'

write.table(rbind(age.df, sex.df), file.path(figdata.dir, 'lung-subgroup.txt'), row.names=FALSE, sep='\t', quote=FALSE)

gg.age <- SubgroupPlot(age.df, age.label.df, margin.lines=c(0.2, 0.2, 0.6, 0.4))
gg.sex <- SubgroupPlot(sex.df, sex.label.df, margin.lines=c(0.2, 0.0, 0.6, 0.6))

subgroup.xmajor <- ggplot_build(gg.sex)$panel$ranges[[1]]$x.major_source
gg.age <- gg.age + scale_x_continuous(breaks=subgroup.xmajor)


OpenPDF('subroup.pdf', width=full.width, height=5)
gridExtra::grid.arrange(gg.sex, gg.age, ncol=2,
  sub=textGrob('Elevation Residuals', vjust=.1),
  left=textGrob('Lung Cancer Incidence Residuals', vjust=1, rot=90))
ClosePDF('subroup.pdf')
options(stringsAsFactors=TRUE)


################################################################################
## Environmental Replacement - BIC Plot

# Evironmental variable correlation plot
cor.list <- OutcomePredictorCorPlot(data.df, c('lung', 'breast'), predictors=envir.covars, 
barwidth=6.1, barheight=0.75, xtext.angle=35, legend.position=c(1, 0.75))
cor.list$plot <- cor.list$plot + theme(axis.text.y=element_text(angle=35, hjust=1)) 
write.table(cor.list$data, file.path(figdata.dir, 'environment-correlation.txt'), 
  row.names=FALSE, sep='\t', quote=FALSE)

# BIC Plot
envir.lung.df <- subset(envir.df, cancer == 'lung')
env.var.sorted <- envir.lung.df$forced[order(envir.lung.df$bayes_factor, decreasing=TRUE)]
gg.bic <- ggplot(subset(envir.df, cancer %in% c('breast', 'lung')), aes(
  forced, log10(bayes_factor)))
gg.bic <- SetGGTheme(gg.bic) +
  geom_hline(y_intercept=1, linetype='dashed', color='darkgrey') +
  geom_point(aes(fill=cancer, size=abs(coefficient), 
    shape=sign(coefficient) > 0, color=cancer)) +
  scale_shape_manual(values=c(25, 24)) +
  scale_alpha_continuous(range=c(0.3, 0.9)) +
  scale_size_continuous(range=c(2.75, 7)) +
  scale_color_manual(values=c(fill.red, fill.blue)) +
  scale_fill_manual(values=c(fill.red, fill.blue)) +
  theme(axis.text.x=element_text(angle=35, hjust=1)) + 
  theme(legend.position='none') +
  xlab('Elevation Replacement') + ylab(expression('log'[10]*'(Bayes Factor)')) +
  scale_x_discrete(limits=env.var.sorted, labels=gsub('_', ' ', env.var.sorted)) +
  scale_y_continuous(breaks=seq(-20, 20, 2), limits=c(-13.2, 2)) +
  theme(plot.margin=grid::unit(c(5, 2, 2, 2), 'points'))

OpenPDF('environment.pdf', width=full.width, height=half.width)
gridExtra::grid.arrange(cor.list$plot, gg.bic, nrow=1)
ClosePDF('environment.pdf')

