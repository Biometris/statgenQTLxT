EM_on_2D_grid <- function(Y,K,X,pen.Cm.grid,pen.Dm.grid,tol.em,max.iter.em) {

# pen.Cm.grid and pen.Dm.grid: vectors defining the values of the penalties for
# the genetic (Cm) and environmental precision matrices
# important : ORDER OF THE VALUES IN THE GRIDS SHOULD BE DESCENDING; IF NOT THEY ARE MADE DESCENDING!

# THE OTHER PARAMETERS ARE THE SAME AS FOR EM_function (WHICH IS CALLED FOR ALL VALUES
# OF THE PENALTIES). IN PARTICULAR:
# Y            : the n x p matrix of phenotypic observations (n individuals, p traits)
#                without missing values; NOT transformed
# K            : the n x n kinship matrix; NOT transformed
# X            : the n x c design matrix  (n individuals, c covariates), can be data.frame()
#                NOT transformed

# For each value in pen.Cm.grid (from large to small), we loop over pen.Dm.grid
# (also from large to small)

# output:

# to do : do not print output to screen

  p1 <- length(pen.Cm.grid)
  p2 <- length(pen.Dm.grid)

  pen.Cm.grid <- sort(pen.Cm.grid,decr=T)
  pen.Dm.grid <- sort(pen.Dm.grid,decr=T)

  convergence.matrix <- matrix(FALSE,p1,p2)

  p <- ncol(Y)

  result.array.Cm       <- array(dim=c(p1,p2,p,p))

  result.array.Dm       <- array(dim=c(p1,p2,p,p))

  for (penC in pen.Cm.grid) {

    # penC <- pen.Cm.grid[1]
    penC.index <- which(penC==pen.Cm.grid)

    has.converged <- T

    penD.index    <- 1

    while (has.converged &   penD.index <= p2) {

      cat('PenC.index = ',penC.index,'\t','PenC.index = ',penD.index,'\n')

      if (penD.index == 1) {
        t <- EM_function(Y=Y,K=K,X=X,Dm.is.diagonal=F,
                            tol.em=tol.em,max.iter.em=max.iter.em,
                            pen.Cm=penC,pen.Dm=pen.Dm.grid[penD.index])
      } else {
        t <- EM_function(Y=Y,K=K,X=X,Dm.is.diagonal=F,
                            tol.em=tol.em,max.iter.em=max.iter.em,
                            pen.Cm=penC,pen.Dm=pen.Dm.grid[penD.index],
                            Cm.start=Cm.start,Dm.start=Dm.start)
      }

      Cm.start <- t$Cm
      Dm.start <- t$Dm

      convergence.matrix[penC.index,penD.index] <- t$converged

      if (t$converged) {
        result.array.Cm[penC.index,penD.index,,] <- t$Cm
        result.array.Dm[penC.index,penD.index,,] <- t$Dm
      } else {
        has.converged <- F
      }

      penD.index    <- penD.index + 1

    }

  }

return(list(convergence.matrix=convergence.matrix,
            result.array.Cm=result.array.Cm,result.array.Dm=result.array.Dm,
            pen.Cm.grid=pen.Cm.grid,pen.Dm.grid=pen.Cm.grid))
}

