EM_on_diag_grid <- function(Y,K,X,pen.Cm.grid,pen.Dm.grid,tol.em,max.iter.em) {

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

# pen.Cm.grid and pen.Dm.grid need to have the same length, say m, and we call EM_function
#   for each of the m pairs. Again, it is assumed that the values in pen.Cm.grid
#   are in descending order

# output:

# to do : do not print output to screen

  p1 <- length(pen.Cm.grid)
  p2 <- length(pen.Dm.grid)

  stopifnot(p1==p2)
  grid.length <- p1

  pen.Cm.grid <- sort(pen.Cm.grid,decr=T)
  pen.Dm.grid <- sort(pen.Dm.grid,decr=T)

  convergence.vector <- rep(FALSE,grid.length)

  p <- ncol(Y)

  result.array.Cm       <- array(dim=c(grid.length,p,p))

  result.array.Dm       <- array(dim=c(grid.length,p,p))

  has.converged <- T

  pen.index    <- 1

  while (has.converged &   pen.index <= grid.length) {

      cat('pen.C = ',pen.Cm.grid[pen.index],'\t','pen.D = ',pen.Cm.grid[pen.index],'\n')

      if (pen.index == 1) {
        t <- EM_function(Y=Y,K=K,X=X,Dm.is.diagonal=F,
                            tol.em=tol.em,max.iter.em=max.iter.em,
                            pen.Cm=pen.Cm.grid[pen.index],
                            pen.Dm=pen.Dm.grid[pen.index])
      } else {
        t <- EM_function(Y=Y,K=K,X=X,Dm.is.diagonal=F,
                            tol.em=tol.em,max.iter.em=max.iter.em,
                            pen.Cm=pen.Cm.grid[pen.index],
                            pen.Dm=pen.Dm.grid[pen.index],
                            Cm.start=Cm.start,Dm.start=Dm.start)
      }

      Cm.start <- t$Cm
      Dm.start <- t$Dm

      convergence.vector[pen.index] <- t$converged

      if (t$converged) {
        result.array.Cm[pen.index,,] <- t$Cm
        result.array.Dm[pen.index,,] <- t$Dm
      } else {
        has.converged <- F
      }

      pen.index    <- pen.index + 1

  }

return(list(convergence.vector=convergence.vector,
            result.array.Cm=result.array.Cm,result.array.Dm=result.array.Dm,
            pen.Cm.grid=pen.Cm.grid,pen.Dm.grid=pen.Cm.grid))
}

