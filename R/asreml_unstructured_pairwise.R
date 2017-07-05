


asreml_unstructured_pairwise <- function(d, K, fix.diag=FALSE,
                                         correlation.matrix=TRUE,
                                         vE.diag=TRUE,
                                         genotype.column,
                                         traitname.column,
                                         phenotype.column,
                                         covariates=integer()) {

#d = Y.long; K=K; fix.diag=FALSE;correlation.matrix=TRUE;vE.diag=vE.diag;genotype.column=1;traitname.column=2;phenotype.column=3;covariates=integer()

# d=long.format.q1; K=K; fix.diag=FALSE; correlation.matrix=TRUE; vE.diag=T; genotype.column=1; traitname.column=3; phenotype.column=4; covariates=integer()

#
#

# TO DO: univariate G-BLUPs; + correlations in case of non-convergence
#        covariates
#        diagonal Ve
#        Checks on data-structure
#        p-values for correlations

### INPUT :

# d : multi-trait or environment-data in long-format
#     If your data are in a wide format, use the melt function in the reshape package

# genotype.column (integer) : a column number, corresponding to the column in d
#                             containing the genotype labels

# traitname.column (integer): a column number, corresponding to the column in d
#                             containing the trait or environment names

# phenotype.column (integer): a column number, corresponding to the column in d
#                             containing the phenotype

# covariates : vector of integers; column numbers,
#              corresponding to the columns in d; no intercept

# K : a kinship matrix ...


### OUTPUT :



################################################################################

  require(asreml)
  require(MASS)

  names(d)[genotype.column]    <- 'genotype'
  names(d)[phenotype.column]   <- 'pheno'
  names(d)[traitname.column]   <- 'trait.name'

  d$trait.name <- as.character(d$trait.name)

  d$genotype <- as.character(d$genotype)

  d <- d[order(d$trait.name),]

  K     <- K[unique(d$genotype), unique(d$genotype)]

  K_inv <- ginv(K) #solve(K)

  rownames(K_inv) <- colnames(K_inv) <- rownames(K)

  trait.names <- unique(d$trait.name)

  n.trait     <- length(trait.names)

  vG.vector   <- rep(0, n.trait)

  vE.vector   <- rep(0, n.trait)

  convergence.vector  <- rep(FALSE, n.trait)

  vG.matrix   <- matrix(0, n.trait, n.trait)

  vE.matrix   <- matrix(0, n.trait, n.trait)

  convergence.matrix <- matrix(FALSE, n.trait, n.trait)

  names(vG.vector) <- names(vE.vector) <- names(convergence.vector) <- trait.names

  colnames(vG.matrix) <- rownames(vG.matrix) <- trait.names

  colnames(vE.matrix) <- rownames(vE.matrix) <- trait.names

  colnames(convergence.matrix) <- rownames(convergence.matrix) <- trait.names

  # First, univariate analyses; estimation of genetic- and residual variances

  if (length(covariates) > 0) {

    fixed.formula <- paste0('pheno ~ ', paste(names(x)[covariates], collapse ='+'))

  } else {

    fixed.formula <- 'pheno ~ 1'

  }

  ft <- call('as.formula',fixed.formula)

  #gi_f <- function(M){M_inv <- ginv(M); rownames(M_inv) <- colnames(M_inv) <- colnames(M); return(M_inv)}
  #gi   <- call('gi_f',K)

  for (i in 1:n.trait) {

    d.red <- d[d$trait.name==trait.names[i], ]

    reml.obj <- asreml(fixed    = eval(ft),
                       random   = ~ giv(genotype),
                       data     = d.red,
                       maxit    =30,
                       workspace=2e+08,
                       ginverse = list(genotype = solve(K)))
                       #ginverse = list(genotype = K_inv))
                       #                        ginverse = list(genotype = eval(gi))

    vG.vector[i] <- summary(reml.obj)$varcomp$component[1]

    vE.vector[i] <- summary(reml.obj)$varcomp$component[2]

    convergence.vector[i] <- reml.obj$converge

  }

  if (correlation.matrix==TRUE) {

    diag(vG.matrix) <- rep(1,nrow(vG.matrix))

    diag(vE.matrix) <- rep(1,nrow(vE.matrix))

  } else {

    diag(vG.matrix) <- vG.vector

    diag(vE.matrix) <- vE.vector

  }


  if (length(covariates) > 0) {

    fixed.formula <- paste0('pheno ~ trait.name + ', paste(names(x)[covariates],collapse ='+'))

  } else {

    fixed.formula <- 'pheno ~ trait.name'

  }

  ft <- call('as.formula',fixed.formula)

  for (i in 1:(n.trait-1)) {

    for (j in (i+1):n.trait) {

      #i=1; j=2
      cat(i,'   ',j,'\n')

      d.red <- d[d$trait.name %in% trait.names[c(i,j)], ]

      if (vE.diag==FALSE) {

        reml.obj <- asreml(fixed = eval(ft),
                           random = ~ us(trait.name):giv(genotype),
                           ginverse = list(genotype = solve(K)),
                           rcov = ~ us(trait.name):id(genotype),
                           data = d.red,
                           maxit=30,
                           workspace =2e+08,
                           start.values=T)


        iv <- reml.obj$gammas

        if (fix.diag) {

          if ( 1 %in% grep(trait.names[i],iv$Gamma)) {
            iv$Value[c(1,3,5,7)] <- c(vG.vector[c(i,j)], vE.vector[c(i,j)])
          } else {
            iv$Value[c(1,3,5,7)] <- c(vG.vector[c(j,i)], vE.vector[c(j,i)])
          }

          iv$Constraint[c(1,3,5,7)] <- 'F'

        }

        iv$Constraint[c(2,6)] <- 'U'

        reml.obj <- asreml(fixed = pheno ~ trait.name,
                           random = ~ us(trait.name):giv(genotype),
                           ginverse = list(genotype = solve(K)),
                           rcov = ~ us(trait.name):id(genotype),
                           data = d.red,
                           maxit = 30,
                           workspace = 2e+08,
                           R.param = iv,
                           G.param = iv)

        summary(reml.obj)$varcomp

        varcov <- matrix(data=c(summary(reml.obj)$varcomp$component[1:2],summary(reml.obj)$varcomp$component[2:3]),nrow=2,ncol=2)

        ve     <- matrix(data=c(summary(reml.obj)$varcomp$component[5:6],summary(reml.obj)$varcomp$component[6:7]),nrow=2,ncol=2)

        if (correlation.matrix==TRUE) {

          varcov <- cov2cor(varcov)

          ve <- cov2cor(ve)

        }

        vG.matrix[i, j] <- vG.matrix[j, i] <- varcov[1,2]

        vE.matrix[i, j] <- vE.matrix[j, i] <- ve[1,2]

      }



      if (vE.diag==TRUE) {

        reml.obj <- asreml(fixed = eval(ft),
                           random = ~ us(trait.name):giv(genotype),
                           ginverse = list(genotype = solve(K)),
                           rcov = ~ at(trait.name):units,
                           data = d.red,
                           maxit=30,
                           workspace =2e+08,
                           start.values=T)

        iv <- reml.obj$gammas

        if (fix.diag) {

          if ( 1 %in% grep(trait.names[i],iv$Gamma)) {
            iv$Value[c(1,3,4,5)] <- c(vG.vector[c(i,j)], vE.vector[c(i,j)])
          } else {
            iv$Value[c(1,3,4,5)] <- c(vG.vector[c(j,i)], vE.vector[c(j,i)])
          }

          iv$Constraint[c(1,3,4,5)] <- 'F'

        }

        iv$Constraint[2] <- 'U'

        reml.obj <- asreml(fixed = pheno ~ trait.name,
                           random = ~ us(trait.name):giv(genotype),
                           ginverse = list(genotype = solve(K)),
                           rcov = ~ at(trait.name):units,
                           data = d.red,
                           maxit = 30,
                           workspace = 2e+08,
                           R.param = iv,
                           G.param = iv)

        varcov <- matrix(data=c(summary(reml.obj)$varcomp$component[1:2],summary(reml.obj)$varcomp$component[2:3]),nrow=2,ncol=2)

        ve <- matrix(data=c(summary(reml.obj)$varcomp$component[4],rep(0,2),summary(reml.obj)$varcomp$component[5]),nrow=2,ncol=2)

        if (correlation.matrix==TRUE) {

          varcov <- cov2cor(varcov)

          ve <- cov2cor(ve)

        }

        vG.matrix[i, j] <- vG.matrix[j, i] <- varcov[1,2]

        vE.matrix[i, j] <- vE.matrix[j, i] <- ve[1,2]

      }

      convergence.matrix[i, j] <- convergence.matrix[j, i] <- reml.obj$converge

    }

  }

  out <- list(vG.matrix = vG.matrix, vE.matrix = vE.matrix,
              convergence.matrix=convergence.matrix,
              vG.vector=vG.vector, vE.vector=vE.vector,
              convergence.vector=convergence.vector
              )

return(out)
}