# the same as test_contrasts , but ...
test_contrasts2      <- function(snp.list,
                                fixed.formula = 'grain.yield ~ Env + S + Env:S',
                                Data=d,
                                random.formula = ~ G + G:EC + G:L + G:Y + G:L:Y,
                                rcov.formula   = ~ at(Env):units,
                                iv,
                                threshold = 0.05,
                                indices=NULL) {

  # Data=d;snp.list=long.list; fixed.formula.init = 'grain.yield ~ Env + S + Env:S'; random.formula = ~ G + G:EC  + G:L + G:Y + G:L:Y ; rcov.formula   = ~ at(Env):units; threshold = 0.05

  # library(asreml); library(MASS); setwd("D:/willem/Dropbox/drops gxe emilie"); load(file='drops_qtl_E_develop2.RData'); iv=g$random.effects; snp.list=as.character(g$snps$snp); fixed.formula='grain.yield ~ Env + S + Env:S'; Data=d; threshold = 0.05; random.formula = ~ G + G:EC  + G:L + G:Y + G:L:Y ; rcov.formula   = ~ at(Env):units

  # iv=g$random.effects; snp.list=as.character(g$snps$snp); Data=d   ; threshold = 0.05; random.formula = ~ G + G:EC + G:L + G:Y + G:L:Y;

  #iv=g$random.effects,
  #                  snp.list=as.character(g$snps$snp),
  #                  fixed.formula=fixed.formula,
  #                  Data=d

  iv$Constraint <- 'F'

  n.snps   <- length(snp.list)

  results  <- data.frame(snp  = as.character(snp.list),
                         p.EC = rep(NA,n.snps),
                         p.EC.Env = rep(NA,n.snps),
                         p.hyd = rep(NA,n.snps),
                         p.hyd.temp = rep(NA,n.snps))

  if (!is.null(indices)) {
    extra.results <- matrix(NA,nrow(results),length(indices) * 2)
    colnames(extra.results) <- paste0('V',1:(length(indices) * 2))
    colnames(extra.results)[2*(1:length(indices)) - 1] <- paste0('p.',indices)
    colnames(extra.results)[2*(1:length(indices))] <- paste0('p.',indices,'RES')
    extra.results <- as.data.frame(extra.results)
    rownames(extra.results) <- rownames(results)
    results <- data.frame(results,extra.results)
    n.indices <- length(indices)
  } else {
    n.indices <- 0
  }

  for (test.snp in snp.list) {

    # test.snp=snp.list[2]

    snp.index <- which(snp.list == test.snp)

    snp.list.temp <- c(snp.list[-snp.index],snp.list[snp.index])

    snp.part <- paste(snp.list.temp[1:(length(snp.list.temp)-1)],collapse=' + ')

    test.snp <- snp.list.temp[length(snp.list.temp)]

    ##########

    snp.contrast1 <- paste0(test.snp,' + ',test.snp,':EC + ',test.snp,':Env')

    fixed.formula1 <- paste0(fixed.formula, ' + Env:(',snp.part,')', ' + ', snp.contrast1)

    ft1 <- call('as.formula',fixed.formula1)

    reml.obj1 <- asreml(fixed  = eval(ft1),
                        random = ~ G + G:EC + G:L + G:Y + G:L:Y,
                        rcov   = ~ at(Env):units,
                        data   = Data,
                        maxit  = 500,
                        workspace   = workspace.size,
                        G.param=iv,
                        R.param=iv)

    w1 <- wald(reml.obj1)

    p1 <- w1[(nrow(w1)-2):(nrow(w1)-1),4]

    results[snp.index,2:3] <- p1

    ##########

    snp.contrast2 <- paste0(test.snp,' + ',test.snp,':water.grouped.scenarios_Gsize',' + ',test.snp,':water.grouped.scenarios_Gsize:temp.water.scenarios_Gsize')

    fixed.formula2 <- paste0(fixed.formula, ' + Env:(',snp.part,')', ' + ', snp.contrast2)

    ft2 <- call('as.formula',fixed.formula2)

    reml.obj2 <- asreml(fixed  = eval(ft2),
                        random = ~ G + G:EC + G:L + G:Y + G:L:Y,
                        rcov   = ~ at(Env):units,
                        data   = Data,
                        maxit  = 500,
                        workspace   = workspace.size,
                        na.method.X="include",
                        G.param=iv,
                        R.param=iv)

    w2 <- wald(reml.obj2)

    p2 <- w2[(nrow(w2)-2):(nrow(w2)-1),4]

    results[snp.index,4:5] <- p2

    if (!is.null(indices)) {

      for (ind in indices) {

       # ind = indices[1]
        snp.contrast.ind <- paste0(test.snp,' + ',test.snp,':',ind, ' + ',test.snp,':Env')

        fixed.formula.ind <- paste0(fixed.formula, ' + Env:(',snp.part,')', ' + ', snp.contrast.ind)

        ft.ind <- call('as.formula',fixed.formula.ind)

        reml.obj.ind <- asreml(fixed  = eval(ft.ind),
                            random = ~ G + G:EC + G:L + G:Y + G:L:Y,
                            rcov   = ~ at(Env):units,
                            data   = Data,
                            maxit  = 500,
                            workspace   = workspace.size,
                            na.method.X="include",
                            G.param=iv,
                            R.param=iv)

        w.ind <- wald(reml.obj.ind)

        p.ind <- w.ind[(nrow(w.ind)-2):(nrow(w.ind)-1),4]

        results[snp.index, 5 + 2*(which(ind==indices) - 1) + 1:2] <- p.ind

      }





    }

  }

return(results)
}

