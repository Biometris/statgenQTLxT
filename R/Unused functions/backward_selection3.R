backward_selection3 <- function(snp.list,
                                fixed.formula.init = 'grain.yield ~ Env + S + Env:S',
                                Data=d,
                                random.formula = ~ G + G:EC + G:L + G:Y + G:L:Y,
                                rcov.formula   = ~ at(Env):units,
                                threshold = 0.05) {

  # Data=d;snp.list=long.list; fixed.formula.init = 'grain.yield ~ Env + S + Env:S'; random.formula = ~ G + G:EC  + G:L + G:Y + G:L:Y ; rcov.formula   = ~ at(Env):units; threshold = 0.05

  continue <- TRUE

  current.snp.list <- snp.list

  current.length   <- length(current.snp.list)

  snp.part <- paste(current.snp.list,collapse=' + ')

  fixed.formula <- paste0(fixed.formula.init, ' + Env:(',snp.part,')')

  ft <- call('as.formula',fixed.formula)

  ###########################
  # initialize variance comp.'s

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

  # at least the first time, asreml should converge ...
  # In all subsequent calls, we use the starting values from the last fit that converged
  stopifnot(reml.obj$converge==TRUE)

  # store the estimated variance components; to be used as starting values in the next asreml call
  iv$Value  <- summary(reml.obj)$varcomp$component

  a <- sort_snps(Data=Data,
                 snp.list=snp.list,
                 fixed.formula.init = fixed.formula.init,
                 random.formula = random.formula,
                 rcov.formula   = rcov.formula,
                 iv=iv)

  ###########

  if (all(a$pval <= threshold)) {
    continue <- FALSE
  }
  #else {
  #  #snp.number.to.be.removed <- which.max(pvalues)
  #  snp.name.to.be.removed   <- as.character(a[nrow(a),1]) #current.snp.list[snp.number.to.be.removed]
  #}

  while (continue) {

    current.snp.list <- as.character(a[1:(nrow(a) - 1),1])

    current.length   <- length(current.snp.list)

    snp.part <- paste(current.snp.list,collapse=' + ')

    fixed.formula <- paste0(fixed.formula.init, ' + Env:(',snp.part,')')

    ft <- call('as.formula',fixed.formula)

    ######################
    # Fit the mixed model

    # first try with the starting values from the previous fit
    reml.obj <- asreml(fixed  = eval(ft),
                       random = random.formula,
                       rcov   = rcov.formula,
                       data   = Data,
                       maxit  = 500,
                       workspace   = 3e+08,
                       G.param=iv,
                       R.param=iv)

    # In case of non-convergence, try with the default starting values
    if (reml.obj$converge==FALSE) {

      reml.obj <- asreml(fixed  = eval(ft),
                         random = random.formula,
                         rcov   = rcov.formula,
                         data   = Data,
                         maxit  = 500,
                         workspace   = 3e+08)
    }

    #############
    # Store the variance components

    if (reml.obj$converge==TRUE) {

      # store the estimated variance components; to be used as starting values in the next asreml call
      iv$Value  <- summary(reml.obj)$varcomp$component

      #w <- wald(reml.obj)

    } else {

      iv$Constraint <- rep('F',nrow(iv))

      reml.obj <- asreml(fixed  = eval(ft),
                         random = random.formula,
                         rcov   = rcov.formula,
                         data   = Data,
                         maxit  = 500,
                         workspace   = 3e+08,
                         G.param=iv,
                         R.param=iv)

      #w <- wald(reml.obj)

      iv$Constraint <- rep('P',nrow(iv))

    }

    ###########################

    a <- sort_snps(Data=Data,
                   snp.list=current.snp.list,
                   fixed.formula.init = fixed.formula.init,
                   random.formula = random.formula,
                   rcov.formula   = rcov.formula,
                   iv=iv)

    ###########

    if (all(a$pval <= threshold)) {
      continue <- FALSE
    }
    #else {
    #  #snp.number.to.be.removed <- which.max(pvalues)
    #  snp.name.to.be.removed   <- as.character(a[nrow(a),1]) #current.snp.list[snp.number.to.be.removed]
    #}


  }

return(list(snps=a, random.effects=iv))
}

