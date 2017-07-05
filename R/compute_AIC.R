compute_AIC           <- function(snp.list,
                                fixed.formula.init = 'grain_yield ~ Env + S + Env:S',
                                Data=d,
                                random.formula = ~ G + G:EC + G:L + G:Y + G:L:Y,
                                rcov.formula   = ~ at(Env):units,
                                iv=NULL,
                                penalty=2) {



  # fixed.formula.init : must be of type character !

  # Data=d;snp.list=long.list; fixed.formula.init = 'grain.yield ~ Env + S + Env:S' ; random.formula = ~ G + G:L + G:Y + G:L:Y ; rcov.formula   = ~ at(Env):units; threshold = 0.05

  # fixed.formula.init = fixed.formula

  snp.part <- paste(snp.list,collapse=' + ')

  fixed.formula <- paste0(fixed.formula.init, ' + Env:(',snp.part,')')

  ft <- call('as.formula',fixed.formula)

  if (is.null(iv)) {

    # 'dummy' call to asreml, to create a data.frame with starting values
    reml.obj <- asreml(fixed  = eval(ft),
                       random = random.formula,
                       rcov   = rcov.formula,
                       data   = Data,
                       maxit  = 500,
                       workspace   = 3e+08,
                       start.values=TRUE)

    iv     <- reml.obj$gammas.table

    # run asreml for the initial model
    reml.obj <- asreml(fixed  = eval(ft),
                       random = random.formula,
                       rcov   = rcov.formula,
                       data   = Data,
                       maxit  = 500,
                       workspace   = 3e+08)

    # asreml should converge ...
    # These starting values are used in all subsequent calls
    stopifnot(reml.obj$converge==TRUE)

      # store the estimated variance components; to be used as starting values in the next asreml call
    iv$Value  <- summary(reml.obj)$varcomp$component

  }
  #else {}

  #iv$Constraint <- 'F'

  AIC_values <- rep(NA,length(snp.list))

  for (snp in snp.list) {

    snp.list.temp <- c(setdiff(snp.list,snp),snp)

    snp.part <- paste(snp.list.temp,collapse=' + ')

    fixed.formula <- paste0(fixed.formula.init, ' + Env:(',snp.part,')')

    ft <- call('as.formula',fixed.formula)

    reml.obj <- asreml(fixed  = eval(ft),
                       random = random.formula,
                       rcov   = rcov.formula,
                       data   = Data,
                       maxit  = 500,
                       workspace   = 3e+08,
                       G.param=iv,
                       R.param=iv)

    w       <- wald(reml.obj)

    AIC_values[which(snp.list==snp)] <-   w[nrow(w)-1, 4]

  }

  result <- data.frame(snp=as.character(snp.list),crit=AIC_values)

  #result <- result[order(result$pval), ]

  row.names(result) <- result$snp

return(result)
}

