library(methods)
library(leaps)
library(glmnet)

model.list <- list()
zmodel.list <- list()
cancer.list <- list()

######################################
## Best Subset and Lasso Regression

for (cancer in cancers) {
  # Create a list to store results
  rsub.list <- list()

  # Create a cancer specific dataset with the included variables
  covariates <- unique(c(global.covars,  specific.covars[[cancer]]))
  cancer.df <- subset(data.df, select=c(cancer, 'elevation', covariates, 'weight'))

  # Remove any counties with missing data
  keep.row <- rowSums(is.na(cancer.df)) == 0
  cancer.df <- cancer.df[keep.row, ]
  rsub.list$cancer.df <- cancer.df
  rsub.list$cancer.allvar.df <- data.df[keep.row, ]

  # Create a standardized dataset (z-scores with weighted mean and sd)
  cancer.zdf <- cancer.df
  for (variable in c(cancer, 'elevation', covariates)) {
    cancer.zdf[, variable] <- WeightedScale(cancer.zdf[, variable], cancer.zdf[, 'weight'])
  }
  rsub.list$cancer.zdf <- cancer.zdf

  # Evaluate all subsets of covariates (best-subset regression)
  X <- subset(cancer.zdf, select = c('elevation', covariates))
  y <- cancer.zdf[, cancer]
  w <- cancer.zdf[, 'weight']
  rsubs <- leaps::regsubsets(x=X, y=y, weights=w,
    force.in='elevation', nvmax=ncol(X), intercept=TRUE)

  rsub.list$rsubs <- rsubs

  # Calculate optimal models
  optimal.index <- which.min(summary(rsubs)$bic)
  predictors <- names(coef(rsubs, optimal.index))[-1]
  rsub.list$predictors <- predictors
  optimal.form <- as.formula(paste(cancer, '~', paste(c('1', predictors), collapse=' + ')))
  model.list[[cancer]]  <- lm(optimal.form, weights=weight, data=cancer.df)
  zmodel.list[[cancer]] <- lm(optimal.form, weights=weight, data=cancer.zdf)


  # Regularized regression
  glmnet.alpha <- 1
  X.mat <- as.matrix(X)
  set.seed(0)
  cv.lasso <- glmnet::cv.glmnet(X.mat, y, w, alpha=glmnet.alpha, standardize=FALSE)
  index.1se <- cv.lasso$lambda == cv.lasso$lambda.1se
  stopifnot(length(cv.lasso$lambda) == length(cv.lasso$glmnet.fit$dev.ratio))
  rsub.list$lasso.r2 <- cv.lasso$glmnet.fit$dev.ratio[index.1se]
  rsub.list$lasso <- cv.lasso

  # Store results
  cancer.list[[cancer]] <- rsub.list
}

##############################################
# Create a best-subset model summary data frame
summary.df <- do.call(rbind, lapply(cancers, function(cancer) {
  predictors <- c('(Intercept)', cancer.list[[cancer]]$predictors)
  mod <- model.list[[cancer]]
  zmod <- zmodel.list[[cancer]]
  conf <- confint(mod, level=0.95)[predictors, ]
  zconf <- confint(zmod, level=0.95)[predictors, ]
  sum.df <- data.frame(
    'cancer'=cancer, 'predictor'=predictors,
    'coef'=coef(mod)[predictors], 'coef_lower'=conf[, 1], 'coef_upper'=conf[, 2],
    'zcoef'=coef(zmod)[predictors], 'zcoef_lower'=zconf[, 1], 'zcoef_upper'=zconf[, 2],
    'pval'=summary(mod)$coef[predictors, 'Pr(>|t|)'], stringsAsFactors=FALSE)
  return(sum.df)
}))

summary.df$mlog_pval <- -log10(summary.df$pval)
summary.df$half_zci <- summary.df$zcoef_upper - summary.df$zcoef
summary.df$mlog_pval <- -log10(summary.df$pval)

write.table(summary.df, file.path(figdata.dir, 'predictors-best-subset.txt'), row.names=FALSE, sep='\t', quote=FALSE)


######################################
## Bivariate and Partial Models

biv.model.list <- list()
par.model.list <- list()
par.df.list <- list()

for (cancer in cancers) {
  cancer.df <- cancer.list[[cancer]]$cancer.allvar.df

  # Calculate bivariate model
  biv.form <- as.formula(paste(cancer, '~ 1 + elevation'))
  biv.model.list[[cancer]] <- lm(biv.form, weights=weight, data=cancer.df)

  # Calculate partial model
  opt.model <- model.list[[cancer]]
  par.df <- data.frame(
    'fips'               = cancer.df$fips,
    'cancer_residual'    = PartialRegressionY(opt.model),
    'elevation_residual' = PartialRegressionX(opt.model),
    'weight'             = cancer.df$weight)
  par.df.list[[cancer]] <- par.df
  par.form <- as.formula('cancer_residual ~ 1 + elevation_residual')
  par.model.list[[cancer]] <- lm(par.form, weights=weight, data=par.df)
}





################################################################################
## Additional Lung Cancer Models
lung.df <- cancer.list[['lung']]$cancer.df
lung.allvar.df <- cancer.list[['lung']]$cancer.allvar.df
lung.zdf <- cancer.list[['lung']]$cancer.zdf
lung.predictors <- cancer.list$lung$predictors

###########################################################
## Subgroup partial models
lung.subgroups <- c('over_65_lung', 'under_65_lung', 'male_lung', 'female_lung')
subgroup.df <- na.omit(data.df[, c(lung.subgroups, 'male_smoking',
  'female_smoking', lung.predictors, 'weight')])
group.df.list <- list()
group.par.list <- list()

for (group in lung.subgroups) {
  variables <- lung.predictors
  if (group == 'male_lung') {
    variables[variables == 'smoking'] <- 'male_smoking'
    variables <- variables[variables != 'male']
  }
  if (group == 'female_lung') {
    variables[variables == 'smoking'] <- 'female_smoking'
    variables <- variables[variables != 'male']
  }
  form <- as.formula(paste(group, '~', paste(c('1', variables), collapse=' + ')))
  model <- lm(form, weights=weight, data=subgroup.df)
  group.df <- data.frame(
    'elevation_residual'=PartialRegressionX(model),
    'incidence_residual'=PartialRegressionY(model),
    'weight'=subgroup.df$weight, 'group'=group, stringsAsFactors=FALSE)
  group.df.list[[group]] <- group.df
  par.model <- lm(incidence_residual ~ elevation_residual, weight=weight, data=group.df)
  group.par.list[[group]] <- par.model
}


###################################
## Radon Analysis
model.list[['radon1']] <- lm(lung ~ radon + smoking, weights=weight, data=lung.df)
model.list[['radon2']] <- lm(lung ~ radon + smoking + elevation, weights=weight, data=lung.df)
predictors3 <- c('radon', lung.predictors)
formula3 <- as.formula(paste('lung', '~', paste(c('1', predictors3), collapse=' + ')))
model.list[['radon3']] <- lm(formula3, weights=weight, data=lung.df)

###################################
## UVB Analysis
model.list[['uvb1']] <- lm(lung ~ uvb + smoking, weights=weight, data=lung.allvar.df)
model.list[['uvb2']] <- lm(lung ~ uvb + smoking + elevation, weights=weight, data=lung.allvar.df)
predictors3 <- c('uvb', lung.predictors)
formula3 <- as.formula(paste('lung', '~', paste(c('1', predictors3), collapse=' + ')))
model.list[['uvb3']] <- lm(formula3, weights=weight, data=lung.allvar.df)


###################################
## Lung Cancer - Smoking * Elevation interaction
interact.formula <- as.formula(paste('lung', '~', paste(c('1', lung.predictors, 'smoking * elevation'), collapse=' + ')))
model.list[['smokingXelevation']] <- lm(interact.formula, weights=weight, data=lung.zdf)



################################################################################
# Substituting elevation with other environmental variables.

BestSubsetMinBIC <- function(bs.df, response, forced, covariates) {
  X <- subset(bs.df, select = c(forced, covariates))
  y <- bs.df[, response]
  w <- bs.df[, 'weight']
  rsubs <- leaps::regsubsets(x=X, y=y, weights=w,
    force.in=forced, nvmax=ncol(X), intercept=TRUE)

  optimal.index <- which.min(summary(rsubs)$bic)
  predictors <- names(coef(rsubs, optimal.index))[-1]
  optimal.form <- as.formula(paste(response, '~', paste(c('1', predictors), collapse=' + ')))
  opt.mod <- lm(optimal.form, weights=weight, data=bs.df)

  min.bic <- list()
  min.bic$bic <- BIC(opt.mod)
  min.bic$model <- opt.mod
  min.bic$env.coef <- coef(opt.mod)[forced]
  min.bic$env.pval <- summary(opt.mod)$coef[forced, 'Pr(>|t|)']
  return(min.bic)
}


envir.df <- do.call(rbind, lapply(cancers, function(cancer) {
  covariates <- unique(c(global.covars, specific.covars[[cancer]]))
  covariates <- setdiff(covariates, envir.covars)
  bs.df <- subset(data.df, select=unique(
    c(cancer, envir.covars, covariates, 'weight')))
  bs.df <- na.omit(bs.df)
  for (j in 2:ncol(bs.df) - 1) {
    bs.df[, j] <- WeightedScale(bs.df[, j], bs.df$weight)
  }
  envir.df <- do.call(rbind, lapply(envir.covars, function(env.var) {
    min.bic <- BestSubsetMinBIC(bs.df, response=cancer,
      forced=env.var, covariates=setdiff(covariates, env.var))
    data.frame('cancer'=cancer, 'counties'=nrow(bs.df), 'forced'=env.var,
      'coefficient'=min.bic$env.coef, 'coef_pval'=min.bic$env.pval,
      'bic'=min.bic$bic, stringsAsFactors=FALSE)
  }))
  envir.df$delta_bic <- envir.df$bic - envir.df[envir.df$forced == 'elevation', 'bic']
  return(envir.df)
}))

envir.df$bayes_factor <- exp(-0.5 * envir.df$delta_bic)

write.table(envir.df, file.path(figdata.dir, 'environment-replacement.txt'), row.names=FALSE, sep='\t', quote=FALSE)



################################################################################
# Elevation nationwide effect.
us.df <- ReadAllCountyData(county.data.path)
elev.coef <- coef(model.list$lung)['elevation']
elev.ci <- confint(model.list$lung, level=0.99)['elevation', ]
max.elev <- max(us.df$elevation, na.rm=TRUE)
DeltaNewCases <- function(elev.coef) {
  sum(elev.coef * (max.elev - us.df$elevation) * us.df$population / 1e5, na.rm=TRUE)}
cat(AddNewLines(sprintf('Holding everything else constant, were all counties to rise in elevation to that of the highest U.S. county, we estimate %.0f [99%% CI: %.0f to %.0f] fewer new lung cancer cases per year would arise.',
  -DeltaNewCases(elev.coef), -DeltaNewCases(elev.ci[2]), -DeltaNewCases(elev.ci[1]) )))
CatDiv()
