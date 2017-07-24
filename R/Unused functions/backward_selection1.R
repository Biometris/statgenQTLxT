backward_selection1 <- function(snp.list,
                                fixed.formula.init = 'grain_yield ~ Env + proba_Oh43 + proba_B14a + proba_Mo17 + proba_W117 + proba_PH207 + Env:(proba_Oh43 + proba_B14a + proba_Mo17 + proba_W117 + proba_PH207)',
                                Data=d,
                                random.formula = ~ G:soilPsi_contrast + G:L + G:Y + G:L:Y,
                                rcov.formula   = ~ at(Env):units,
                                threshold = 0.05) {

  # Data=d;snp.list=long.list; fixed.formula.init = 'grain_yield ~ Env + proba_Oh43 + proba_B14a + proba_Mo17 + proba_W117 + proba_PH207 + Env:(proba_Oh43 + proba_B14a + proba_Mo17 + proba_W117 + proba_PH207)';random.formula = ~ G:soilPsi_contrast + G:L + G:Y + G:L:Y ; rcov.formula   = ~ at(Env):units; threshold = 0.05

  continue <- TRUE

  current.snp.list <- snp.list

  current.length   <- length(current.snp.list)

  snp.part <- paste(current.snp.list,collapse=' + ')

  fixed.formula <- paste0(fixed.formula.init, ' + Env:(',snp.part,')')

  ft <- call('as.formula',fixed.formula)

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

  # Number of iterations
  n.iter <- 1
  #summary(reml.obj)$varcomp ; eff <- coef(reml.obj)$fixed

  w       <- wald(reml.obj)

  row.ind <- (nrow(w)-current.length):(nrow(w)-1)

  rn      <- rownames(w)[row.ind]

  pvalues <- w[row.ind, 4]

  if (all(pvalues <= threshold)) {
    continue <- FALSE
  } else {
    snp.number.to.be.removed <- which.max(pvalues) #nrow(w) - current.length -1 + which.max(pvalues)
    snp.name.to.be.removed   <- current.snp.list[snp.number.to.be.removed]
  }

  while (continue) {

    current.snp.list <- setdiff(current.snp.list,snp.name.to.be.removed)

    current.length   <- length(current.snp.list)

    snp.part <- paste(current.snp.list,collapse=' + ')

    fixed.formula <- paste0(fixed.formula.init, ' + Env:(',snp.part,')')

    ft <- call('as.formula',fixed.formula)

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

    n.iter <- n.iter + 1

    if (reml.obj$converge==TRUE) {

      # store the estimated variance components; to be used as starting values in the next asreml call
      iv$Value  <- summary(reml.obj)$varcomp$component

      w <- wald(reml.obj)

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

      w <- wald(reml.obj)

      iv$Constraint <- rep('P',nrow(iv))

    }

    cat(n.iter,'\n\n')

    # summary(reml.obj)$varcomp
    # eff <- coef(reml.obj)$fixed

    row.ind <- (nrow(w)-current.length):(nrow(w)-1)

    rn <- rownames(w)[row.ind]

    pvalues <- w[row.ind, 4]

    if (all(pvalues <= threshold)) {
      continue <- FALSE
    } else {
      snp.number.to.be.removed <- which.max(pvalues)
      snp.name.to.be.removed   <- current.snp.list[snp.number.to.be.removed]
    }
    #snp.name.to.be.removed
  }

return(list(fixed.effects=w, random.effects=iv, final.snp.list=current.snp.list))
}

##########################################################

