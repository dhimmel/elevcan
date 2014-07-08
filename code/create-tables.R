library(Hmisc)

################################################################################
############################ Function Declaration ###########################
LatexMeanSD <- function(x, w, decimals.mean=2, decimals.sd=2) {
  response.mean <- Hmisc::wtd.mean(x, w)
  response.sd <- sqrt(Hmisc::wtd.var(x, w, normwt=TRUE))
  chr.mean <- ChrRound(response.mean, decimals.mean)
  chr.sd <- ChrRound(response.sd, decimals.sd)
  chr.mean.sd <- sprintf('%s (%s)', chr.mean, chr.sd)
  return(chr.mean.sd)
}

LatexPercent <- function(value, decimals=1) {
  value.percent <- ChrRound(value * 100, decimals)
  percent.tex <- sprintf('%s\\%%', value.percent)
  return(percent.tex)
}

LatexCoefficient <- function(value, decimals=2) {
  return(ChrRound(value, decimals))
}

LatexPvalue <- function(pval) {
  # Latex inline math scientific formation.
  sci.notation <- sprintf("%.2e}$", pval)
  return(sub('e', '$\\\\times10^{', sci.notation))
}

LatexCoefCI <- function(coefficient, lower, upper, decimals=2) {
  coefficient <- LatexCoefficient(coefficient, decimals)
  lower <- LatexCoefficient(lower, decimals)
  upper <- LatexCoefficient(upper, decimals)
  chr <- sprintf('%s [%s, %s]', coefficient, lower, upper)
  return(chr)
}


NegativeCoefTest <- function(model, predictor) {
  tstat <- summary(model)$coef[predictor, 't value']
  pval <- pt(tstat, df=model$df.residual, lower.tail=TRUE)
  #cat(sprintf('Model %s: %s < 0 p-value %.3e\n', name, predictor, pval))
  return(pval)
}

Model2LatexCoefCI <- function(model, predictor, level=0.95, standardize=FALSE, decimals=2) {
  coefficient <- coef(model)[predictor]
  conf.vector <- confint(model, level=level)[predictor, ]

  if (standardize) {
    model.df <- model.frame(model)
    if (! any('(weights)' == colnames(model.df))) 
      {model.df[, '(weights)'] <- 1}
    w <- model.df[, '(weights)']
    response.sd <- Hmisc::wtd.var(model.df[, 1], weights=w, normwt=TRUE) ^ 0.5
    predictor.sd <- Hmisc::wtd.var(model.df[, predictor], weights=w, normwt=TRUE) ^ 0.5

    coefficient <- StandardizeBeta(coefficient, predictor.sd, response.sd)
    conf.vector <- StandardizeBeta(conf.vector, predictor.sd, response.sd)
  }
  chr <- LatexCoefCI(coefficient, conf.vector[1], conf.vector[2], decimals)
  return(chr)
}

################################################################################
################################# Execution ####################################

## Best Subset Summary Table
summary.texdf <- do.call(rbind, lapply(cancers, function(cancer) {
  cancer.df <- cancer.list[[cancer]]$cancer.df
  predictors <- cancer.list[[cancer]]$predictors
  model <- model.list[[cancer]]
  zmodel <- zmodel.list[[cancer]]

  elev.mean <- weighted.mean(cancer.df[, cancer], cancer.df$weight)
  elev.coef <- coef(model)['elevation']
  elev.tstat <- summary(model)$coef['elevation', 't value']
  elev.pval <- pt(elev.tstat, df=model$df.residual, lower.tail=TRUE)
  elev.coef.percent <- elev.coef / elev.mean
  conf.int <- confint(model, level=0.99)['elevation', ]
  zconf.int <- confint(zmodel, level=0.99)['elevation', ]
  pconf.int <- conf.int / elev.mean * 100
  conf.int <- sapply(conf.int, LatexCoefficient)
  zconf.int <- sapply(zconf.int, LatexCoefficient)
  pconf.int <- sapply(pconf.int, LatexCoefficient)

  summary.row <- data.frame('cancer'=cancer,
    'mean (sd)'=LatexMeanSD(cancer.df[, cancer], cancer.df[, 'weight'], 1, 1),
    'n'=nrow(cancer.df), 'size'=length(predictors), 
    '$R^2$'=LatexPercent(summary(model)$r.squared),
    '$p$'=LatexPvalue(elev.pval),
    '$\\beta$'=LatexCoefficient(elev.coef),
    '$\\beta_z$'=LatexCoefficient(coef(zmodel)['elevation']),
    '$\\beta_{\\%}$'=LatexPercent(elev.coef.percent),
    check.names=FALSE, stringsAsFactors=FALSE, row.names=cancer)

  summary.row <- rbind(summary.row, data.frame(
    'cancer'='', 'mean (sd)'='', 'n'='', 'size'='', '$R^2$'='', '$p$'='',
    '$\\beta$'=sprintf('[%s, %s]', conf.int[1], conf.int[2]),
    '$\\beta_z$'=sprintf('[%s, %s]', zconf.int[1], zconf.int[2]),
    '$\\beta_{\\%}$'=sprintf('[%s, %s]\\%%', pconf.int[1], pconf.int[2]),
    check.names=FALSE, stringsAsFactors=FALSE, row.names=cancer))

  return(summary.row)
}))

Hmisc::latex(summary.texdf, file=file.path(table.dir, 'best-subset-summaries.tex'),
  rowname=NULL, label='tab:best-subset', caption='Optimal Best Subset Model Summaries')

## Lasso Summary Table
lasso.texdf <- do.call(rbind, lapply(cancers, function(cancer) {
  cancer.df <- cancer.list[[cancer]]$cancer.df
  w <- cancer.df[, 'weight']
  predictor.sd <- sqrt(Hmisc::wtd.var(cancer.df[, 'elevation'], w, normwt=TRUE))
  response.sd <- sqrt(Hmisc::wtd.var(cancer.df[, cancer], w, normwt=TRUE))

  elev.zcoef <- coef(cancer.list[[cancer]]$lasso)['elevation', ]
  elev.coef <- UnstandardizeBeta(elev.zcoef, predictor.sd, response.sd)
  elev.coef.percent <- elev.coef / weighted.mean(cancer.df[, cancer], cancer.df$weight)

  summary.row <- data.frame(
    'cancer'=cancer,
    'size'=sum(lasso.df$Cancer == SimpleCap(cancer)),
    '$R^2$'=LatexPercent(cancer.list[[cancer]]$lasso.r2),
    '$\\beta$'=LatexCoefficient(elev.coef), 
    '$\\beta_z$'=LatexCoefficient(elev.zcoef), 
    '$\\beta_{\\%}$'=LatexPercent(elev.coef.percent),
    check.names=FALSE, stringsAsFactors=FALSE, row.names=cancer)
  return(summary.row)
}))

Hmisc::latex(lasso.texdf, file=file.path(table.dir, 'lasso-summaries.tex'),
  rowname=NULL, label='tab:lasso', caption='Lasso Model Summaries')


## Predictor Table
all.predictors <- unique(c(all.covars, envir.covars))
predictor.texdf <- data.frame(
  'predictor'=all.predictors,
  'cancers'='CANCERS', stringsAsFactors=FALSE)
predictor.texdf[predictor.texdf$predictor %in% global.covars, 'cancers'] <- 'all'

predictor.texdf[, 'n'] <- sapply(predictor.texdf$predictor, function(predictor) {
  sum(!is.na(data.df[, predictor]))})
predictor.texdf[, 'mean'] <- sapply(predictor.texdf$predictor, function(predictor) {
  format(mean(data.df[, predictor], na.rm=TRUE), digits=2)})
predictor.texdf[, 'sd'] <- sapply(predictor.texdf$predictor, function(predictor) {
  format(sd(data.df[, predictor], na.rm=TRUE), digits=2)})

predictor.texdf[, 'years'] <- 'YEARS'
predictor.texdf[, 'units'] <- 'UNITS'

predictor.texdf$predictor <- gsub('_', ' ', predictor.texdf$predictor)

Hmisc::latex(predictor.texdf, file=file.path(table.dir, 'predictors.tex'),
  rowname=NULL, label='tab:predictors', caption='Predictor Table')


## Best Subset Model Tables
for (cancer in cancers) {
  csum.df <- summary.df[summary.df$cancer == cancer, ]
  csum.df <- csum.df[order(csum.df$predictor  == '(Intercept)', csum.df$pval), ]
  csum.df$predictor <- gsub('_', ' ', csum.df$predictor)
  csum.df$predictor <- factor(csum.df$predictor, levels=csum.df$predictor)

  tex.df <- plyr::ddply(csum.df, 'predictor', plyr::summarize, 
    '$\\beta$'=LatexCoefCI(coef, coef_lower, coef_upper),
    '$\\beta_z$'=LatexCoefCI(zcoef, zcoef_lower, zcoef_upper),
    'p-value'=LatexPvalue(pval))

  path <- file.path(table.dir, sprintf('SI_best-subset-model-%s.tex', cancer))
  cancer.caption <- sprintf('Optimal best subset regression model for %s cancer', cancer)
  table.label <- sprintf('SI_tab:%s', cancer)
  Hmisc::latex(tex.df, file=path, rowname=NULL, 
    caption=cancer.caption, label=table.label)
}



# Radon/UVB Confounding Effects of Elevation Table
radon.list <- model.list[c('radon1', 'radon2', 'radon3')]
uvb.list <- model.list[c('uvb1', 'uvb2', 'uvb3')]
confounding.texdf <- data.frame('model'=1:3,
  #'radon $\\beta$'=sapply(radon.list, Model2LatexCoefCI, predictor='radon'),
  'radon $\\beta_z$'=sapply(radon.list, Model2LatexCoefCI, predictor='radon', standardize=TRUE),
  'radon p-value'=sapply(sapply(radon.list, NegativeCoefTest, predictor='radon'), LatexPvalue),
  #'uvb $\\beta$'=sapply(uvb.list, Model2LatexCoefCI, predictor='uvb'),
  'uvb $\\beta_z$'=sapply(uvb.list, Model2LatexCoefCI, predictor='uvb', standardize=TRUE),
  'uvb p-value'=sapply(sapply(uvb.list, NegativeCoefTest, predictor='uvb'), LatexPvalue),
  check.names=FALSE, stringsAsFactors=FALSE) 

Hmisc::latex(confounding.texdf, file=file.path(table.dir, 'confounding.tex'),
  rowname=NULL, label='tab:confounding', caption='Confounding effect of elevation on radon and UVB lung cancer association.')


