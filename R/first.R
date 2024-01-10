twin.sample <- function(X, n) {
  N <- nrow(X)
  r <- floor(N/n)
  if (r > 1) {
    rind <- sample(1:N, (N-r*n))
    kind <- setdiff(1:N, rind)
    ind <- kind[twinning::twin(X[kind,,drop=FALSE],r)]
  } else {
    ind <- sample(1:N, n)
  }
  return (ind)
}

preproess.input <- function(X) {
  # one hot encoding of categorical features
  p <- ncol(X)
  ind.cat <- which(sapply(X,class)=="factor")
  if (length(ind.cat) > 0) {
    ind.con <- setdiff(1:p, ind.cat)
    X.tmp <- as.matrix(X[,ind.con,drop=FALSE])
    for (i in 1:length(ind.cat)) {
      # normalization to have unit euclidean distance in difference
      ftr.tmp <- stats::model.matrix(~.-1,X[,ind.cat[i],drop=FALSE])
      X.tmp <- cbind(X.tmp, ftr.tmp)
    }
    X <- X.tmp
  } else {
    X <- as.matrix(X,ncol=p)
  }
  
  return (X)
}

exp.var.knn <- function(X, y, subset, n.knn=2, n.mc=nrow(X), mc.option="random") {
  # preprocess
  if (is.logical(subset[1])) {
    p <- ncol(X)
    if (length(subset)==p) {
      subset <- c(1:p)[subset]
    } else {
      stop(sprintf("length of subset(TRUE/FALSE) is not %d", p))
    }
  }
  X <- preproess.input(X[,subset,drop=FALSE])
  
  # compute D.ind
  if (n.mc < nrow(X)){
    if (mc.option == "random"){
      D.ind <- sample(1:nrow(X), n.mc)
    } else if (mc.option == "twinning") {
      D.ind <- twin.sample(X, n.mc)
    } else {
      stop(sprintf("no mc.option %s. It can either be random or twinning.", mc.option))
    }
  } else {
    D.ind <- c(1:nrow(X))
  }
  
  # find the n.knn cloeset points
  if (nrow(unique(X)) < nrow(X)) {
    # compute exact distance
    D <- X[D.ind,,drop=FALSE]
    suppressWarnings(D.dist <- as.matrix(pdist::pdist(D,X)))
    c.var <- numeric(length(D.ind))
    for (i in 1:length(D.ind)){
      d.dist <- D.dist[i,]
      rad <- d.dist[order(d.dist)[n.knn]]
      id.sort <- c(1:length(d.dist))[d.dist<=rad]
      c.var[i] <- stats::var(y[id.sort])
    }
    ev <- mean(c.var)
  } else {
    # use fast knn
    D <- X[D.ind,,drop=FALSE]
    id.sort <- FNN::knnx.index(data=X, query=D, k=n.knn)
    ev <- mean(apply(matrix(y[id.sort],ncol=n.knn),1,stats::var))
  }
  
  return (ev)
}

#' @title
#' Estimating Total Sobol' Indices from Data 
#' 
#' @description 
#' \code{totalsobol.knn} implements the estimation of the total Sobol' indices directly from scattered data.
#' Parallel computation is available to accelerate the estimation.  
#' For categorical inputs, please convert them to factor before calling this function. 
#' For large datasets, we support the use of subsample to reduce the computational cost.
#'   
#' @importFrom FNN knnx.index
#' @importFrom pdist pdist
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#' @importFrom stats var model.matrix
#' @importFrom twinning twin
#' 
#' @param X a matrix or data frame containing the inputs.
#' @param y a vector containing the outputs.
#' @param noise a logical indicating whether the outputs are noisy. 
#' @param n.knn number of nearest-neighbors for the inner loop conditional variance estimation.
#' @param n.mc number of Monte Carlo samples for the outer loop expectation estimation.
#' @param mc.option option for the Monte Carlo samples. The options are random(default)/twinning. See below for details.
#' @param rescale a logical indicating whether continuous inputs is rescaled to zero mean unit variance.
#' @param parl number of cores on which to parallelize the computation. If \code{NULL}, then no parallelization is done.
#' 
#' @details  
#' \code{totalsobol.knn} provides consistent estimation of the total Sobol' indices (Sobol, 1993) from scattered data.
#' When the output is clean/noiseless (\code{noise=FALSE}), \code{totalsobol.knn} implements the nearest-neighbor estimator proposed in Broto et al. (2020).
#' See \code{shapleysobol_knn} from the \pkg{sensitivity} package for another implementation of the nearest-neighbor estimator. 
#' When the output is noisy (\code{noise=TRUE}), \code{totalsobol.knn} implements the Noise-Adjusted Nearest-Neighbor (NANNE) estimator (Huang and Joseph, 2024). 
#' NANNE estimator can correct the estimation bias of the nearest-neighbor estimator caused by the random noise. 
#' Please see Huang and Joseph (2024) for a more detailed discussion and comparison. 
#' 
#' For integer/numeric output, \code{n.knn=2} nearest-neighbors is recommended for the noisy data (Huang and Joseph, 2024),
#' and \code{n.knn=3} nearest-neighbors is suggested for the clean/noiseless data (Broto et al., 2020). 
#' For numeric inputs, it is recommended to standardize them via setting the argument \code{rescale=TRUE}.
#' Categorical inputs are transformed via one-hot encoding for the nearest-neighbor search. 
#' To speed up the nearest-neighbor search, k-d tree from the \pkg{FNN} package is used. 
#' Also, parallel computation is also supported via the \pkg{parallel} package.
#' 
#' Last, for large datasets, we support the use of subsamples for further acceleration.
#' Use argument \code{n.mc} to specify the number of subsamples. 
#' Two options (argument \code{mc.option}) are available for finding the subsamples: random and twinning (Vakayil and Joseph, 2022). 
#' Twinning is able to find subsamples that better represent the big data, i.e., 
#' providing a more accurate estimation, but at a slightly higher computational cost. 
#' For more details, please see the \pkg{twinning} package.
#' 
#' @return 
#' A numeric vector of the total Sobol' indices estimation.
#' 
#' @author 
#' Chaofan Huang \email{chaofan.huang@@gatech.edu} and V. Roshan Joseph \email{roshan@@gatech.edu}
#' 
#' @references
#' Broto, B., Bachoc, F., & Depecker, M. (2020). Variance reduction for estimation of Shapley effects and adaptation to unknown input distribution. SIAM/ASA Journal on Uncertainty Quantification, 8(2), 693-716. 
#' 
#' Huang, C., & Joseph, V. R. (2024). Factor Importance Ranking and Selection using Total Indices. arXiv preprint arXiv:2401.00800.
#' 
#' Sobol', I. Y. M. (1990). On sensitivity estimation for nonlinear mathematical models. Matematicheskoe modelirovanie, 2(1), 112-118.
#' 
#' Vakayil, A., & Joseph, V. R. (2022). Data twinning. Statistical Analysis and Data Mining: The ASA Data Science Journal, 15(5), 598-610.
#' 
#' @export
#' 
#' @examples
#' ishigami <- function(x) {
#'   x <- -pi + 2*pi*x
#'   y <- sin(x[1]) + 7*sin(x[2])^2 + 0.1*x[3]^4*sin(x[1])
#'   return (y)
#' }
#' 
#' set.seed(123)
#' n <- 10000
#' p <- 3
#' X <- matrix(runif(n*p), ncol=p)
#' y <- apply(X,1,ishigami) + rnorm(n)
#' tsi <- totalsobol.knn(X, y, noise=TRUE, n.knn=2, rescale=FALSE)
#' print(round(tsi,3)) # Analytical Total Sobol' Indices: 0.558, 0.442, 0.244
#' 

totalsobol.knn <- function(X, y, noise=TRUE, n.knn=2, n.mc=nrow(X), mc.option="random", 
                           rescale=TRUE, parl=NULL) {
  # arguments check
  if(!(inherits(X, "matrix") | inherits(X, "data.frame")))
    stop("X must be either a matrix of a data frame.")
  if(n.mc <= 0) 
    stop("n.mc must be a positive integer.")
  
  # preprocess
  n <- nrow(X)
  p <- ncol(X)
  n.mc <- min(n.mc, n)
  if (inherits(X, "matrix")) X <- data.frame(X)
  # verify column type
  col.class.error <- F
  col.class.error.msg <- "\n"
  col.class <- sapply(X,class)
  for (i in 1:p) {
    if (!(col.class[i] %in% c("integer","numeric","factor"))) {
      col.class.error.msg <- paste(
        col.class.error.msg, 
        sprintf("Column %d (%s) is type %s.\n", i, colnames(X)[i], col.class[i]),
        sep="")
      col.class.error <- T
    }
  }
  if (col.class.error) {
    col.class.error.msg <- paste(
      col.class.error.msg, 
      "Please convert above columns to the supported column types (integer/numeric/factor).\n",
      sep="")
    stop(col.class.error.msg)
  }
  # rescale continuous inputs 
  if (rescale) {
    for (i in 1:p) {
      if (col.class[i]=="numeric")
        X[,i] <- c(scale(X[,i]))
    }
  }
  y <- c(y)
  
  # compute the variance for the noise
  noise.var <- if(noise) exp.var.knn(X=X,y=y,subset=rep(TRUE,p),n.knn=n.knn,n.mc=n.mc,mc.option=mc.option) else 0
  y.var <- pmax(stats::var(y)-noise.var, 0)
  if (y.var==0) return (rep(0,p))
  if (p==1) return (c(1))
  
  # compute total sobol' indices, with last term being the noise
  indices <- (diag(p)==0)
  if (!is.null(parl)) {
    clus <- parallel::makeCluster(parl)
    parallel::clusterExport(clus, list("exp.var.knn","preproess.input"))
    tsi <- parallel::parApply(clus, indices, 1, 
                              function(ss) exp.var.knn(X=X,y=y,subset=ss,n.knn=n.knn,n.mc=n.mc,mc.option=mc.option))
    parallel::stopCluster(clus)
  } else {
    tsi <- apply(indices, 1, function(ss) exp.var.knn(X=X,y=y,subset=ss,n.knn=n.knn,n.mc=n.mc,mc.option=mc.option))
  }
  tsi <- pmax(tsi-noise.var, 0) / y.var
  
  return (tsi)
}

#' @title NANNE estimator with Backward Elimination
#' @importFrom FNN knnx.index
#' @importFrom pdist pdist
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#' @importFrom stats var model.matrix
#' @importFrom twinning twin
#' @noRd
nanne.be <- function(X, y, n.knn=2, n.mc=nrow(X), mc.option="random", 
                     rescale=TRUE, parl=NULL, verbose=FALSE) {
  # arguments check
  if(!(inherits(X, "matrix") | inherits(X, "data.frame")))
    stop("X must be either a matrix of a data frame.")
  if(n.mc <= 0) 
    stop("n.mc must be a positive integer.")
  
  # preprocess
  n <- nrow(X)
  p <- ncol(X)
  n.mc <- min(n.mc, n)
  if (inherits(X, "matrix")) X <- data.frame(X)
  # verify column type
  col.class.error <- FALSE
  col.class.error.msg <- "\n"
  col.class <- sapply(X,class)
  for (i in 1:p) {
    if (!(col.class[i] %in% c("integer","numeric","factor"))) {
      col.class.error.msg <- paste(
        col.class.error.msg, 
        sprintf("Column %d (%s) is type %s.\n", i, colnames(X)[i], col.class[i]),
        sep="")
      col.class.error <- TRUE
    }
  }
  if (col.class.error) {
    col.class.error.msg <- paste(
      col.class.error.msg, 
      "Please convert above columns to the supported column types (integer/numeric/factor).\n",
      sep="")
    stop(col.class.error.msg)
  }
  # rescale continuous inputs 
  if (rescale) {
    for (i in 1:p) {
      if (col.class[i]=="numeric")
        X[,i] <- c(scale(X[,i]))
    }
  }
  y <- c(y)
  
  # iterate until convergence
  subset <- c(1:p)
  while (TRUE) {
    subset.old <- subset
    subset.imp <- totalsobol.knn(X[,subset,drop=FALSE],y,noise=TRUE,
                                 n.knn=n.knn,n.mc=n.mc,mc.option=mc.option,
                                 rescale=FALSE,parl=parl)
    if (verbose) {
      print(noquote(paste(c("Subset:",subset), collapse=" ")))
      print(noquote(paste(c("Total Sobol' Indices:",round(subset.imp,3)), collapse=" ")))
    }
    subset <- subset[subset.imp>0]
    subset.imp <- subset.imp[subset.imp>0]
    
    if (length(subset)==0) break
    if (length(subset.old)==length(subset)) break
  }
  
  imp <- rep(0,p)
  if (length(subset)>0) imp[subset] <- subset.imp
  return (imp)
}

#' @title 
#' Factor Importance Ranking and Selection using Total (Sobol') indices
#' 
#' @description 
#' \code{first} implements the \emph{\strong{model-independent}} factor importance and selection procedure proposed in Huang and Joseph (2024).
#' The importance measure is based on total Sobol' indices from global sensitivity analysis. 
#' Factor importance computation and selection are performed directly from the noisy data.
#' Parallel computations are available to accelerate the estimation. 
#' For categorical data inputs, please convert them to factor type before calling the function.
#' 
#' @importFrom FNN knnx.index
#' @importFrom pdist pdist
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#' @importFrom stats var model.matrix
#' @importFrom twinning twin
#' 
#' @param X a matrix or data frame containing the inputs.
#' @param y a vector containing the outputs.
#' @param n.knn number of nearest-neighbors for the inner loop conditional variance estimation. See \code{\link{totalsobol.knn}} for details.
#' @param n.mc number of Monte Carlo samples for the outer loop expectation estimation. See \code{\link{totalsobol.knn}} for details.
#' @param mc.option option for the Monte Carlo samples. The options are random(default)/twinning. See \code{\link{totalsobol.knn}} for details.
#' @param rescale a logical indicating whether continuous inputs is rescaled to zero mean unit variance.
#' @param sparse a logical indicating whether to run the fast version of forward selection that integrates with the effect sparsity principle.
#' @param parl number of cores on which to parallelize the computation. If \code{NULL}, then no parallelization is done.
#' @param verbose a logical indicating whether to display intermediate results of running \code{first}.
#' 
#' @details
#' \code{first} provides factor importance ranking and selection directly from scattered data without any model fitting. 
#' Factor importance is based on total Sobol' indices (Sobol, 1993), 
#' and its connection to the approximation error introduced by excluding the factor of interest is shown in Huang and Joseph (2024).
#' The estimation of the total Sobol' indices is carried out using \code{\link{totalsobol.knn}}. 
#' Integrating it with forward selection and backward elimination allows for factor selection.     
#' Please see Huang and Joseph (2024) for the details of the algorithm. 
#' 
#' \code{n.knn=2} nearest-neighbors is recommended for integer/numeric output, and \code{n.knn=3} is suggested for binary output.
#' For numeric inputs, it is recommended to standardize them via setting the argument \code{rescale=TRUE}.
#' Categorical inputs are transformed via one-hot encoding for the nearest-neighbor search. 
#' To speed up the nearest-neighbor search, k-d tree from the \pkg{FNN} package is used. 
#' Also, parallel computation is also supported via the \pkg{parallel} package.
#' 
#' For very large dimensional problems, we also provide the fast forward selection method (\code{sparse=TRUE}) 
#' that is based on the effect sparsity principle, but it sacrifices some accuracy.
#' See Appendix of Huang and Joseph (2024) for more details. 
#' For large datasets, we support the use of subsamples for further acceleration.
#' Use argument \code{n.mc} to specify the number of subsamples. 
#' Two options (argument \code{mc.option}) are available for finding the subsamples: random and twinning (Vakayil and Joseph, 2022). 
#' Twinning is able to find subsamples that better represent the big data, i.e., 
#' providing a more accurate estimation, but at a slightly higher computational cost. 
#' For more details, please see the \pkg{twinning} package.
#' 
#' @return
#' A numeric vector of the factor importance, with zero indicating that the factor is not important to the prediction of the response.
#' 
#' @author 
#' Chaofan Huang \email{chaofan.huang@@gatech.edu} and V. Roshan Joseph \email{roshan@@gatech.edu}
#' 
#' @references
#' Huang, C., & Joseph, V. R. (2024). Factor Importance Ranking and Selection using Total Indices. arXiv preprint arXiv:2401.00800.
#' 
#' Sobol', I. Y. M. (1990). On sensitivity estimation for nonlinear mathematical models. Matematicheskoe modelirovanie, 2(1), 112-118.
#' 
#' Vakayil, A., & Joseph, V. R. (2022). Data twinning. Statistical Analysis and Data Mining: The ASA Data Science Journal, 15(5), 598-610.
#' 
#' @export
#'
#' @examples
#' ishigami <- function(x) {
#'   x <- -pi + 2*pi*x
#'   y <- sin(x[1]) + 7*sin(x[2])^2 + 0.1*x[3]^4*sin(x[1])
#'   return (y)
#' } 
#‘
#' set.seed(123)
#' n <- 1000
#' p <- 6
#' X <- matrix(runif(n*p), ncol=p)
#' y <- apply(X,1,ishigami) + rnorm(n)
#' imp <- first(X, y, n.knn=2, rescale=FALSE, verbose=TRUE)
#' print(round(imp,3)) # Only first 3 factors are important
#' 

first <- function(X, y, n.knn=2, n.mc=nrow(X), mc.option="random", 
                  rescale=TRUE, sparse=FALSE, parl=NULL, verbose=FALSE) {
  # arguments check
  if(!(inherits(X, "matrix") | inherits(X, "data.frame")))
    stop("X must be either a matrix of a data frame.")
  if(n.mc <= 0) 
    stop("n.mc must be a positive integer.")
  
  # preprocess
  n <- nrow(X)
  p <- ncol(X)
  n.mc <- min(n.mc, n)
  if (inherits(X, "matrix")) X <- data.frame(X)
  # verify column type
  col.class.error <- FALSE
  col.class.error.msg <- "\n"
  col.class <- sapply(X,class)
  for (i in 1:p) {
    if (!(col.class[i] %in% c("integer","numeric","factor"))) {
      col.class.error.msg <- paste(
        col.class.error.msg, 
        sprintf("Column %d (%s) is type %s.\n", i, colnames(X)[i], col.class[i]),
        sep="")
      col.class.error <- TRUE
    }
  }
  if (col.class.error) {
    col.class.error.msg <- paste(
      col.class.error.msg, 
      "Please convert above columns to the supported column types (integer/numeric/factor).\n",
      sep="")
    stop(col.class.error.msg)
  }
  # rescale continuous inputs 
  if (rescale) {
    for (i in 1:p) {
      if (col.class[i]=="numeric")
        X[,i] <- c(scale(X[,i]))
    }
  }
  y <- c(y)
  
  # iterate forward selection
  if (verbose) print(noquote("Starting forward selection..."))
  y.var <- stats::var(y)
  subset <- c()
  candidate <- c(1:p)
  x.isi.max <- 0
  while (TRUE) {
    if (length(candidate)==0) break
    # compute total sobol' index for noise
    if (!is.null(parl)) {
      clus <- parallel::makeCluster(parl)
      parallel::clusterExport(clus, list("exp.var.knn","preproess.input"))
      e.tsi <- parallel::parSapply(clus, candidate, function(i) {
        exp.var.knn(X=X,y=y,subset=c(subset,i),n.knn=n.knn,n.mc=n.mc,mc.option=mc.option)
      })
      parallel::stopCluster(clus)
    } else {
      e.tsi <- sapply(candidate, function(i) {
        exp.var.knn(X=X,y=y,subset=c(subset,i),n.knn=n.knn,n.mc=n.mc,mc.option=mc.option)
      })
    }
    # compute inclusive sobol' index x in the model
    x.isi <- pmax(y.var-e.tsi, 0)
    if (verbose) {
      print(noquote(paste(c("Candidate:",candidate), collapse=" ")))
      print(noquote(paste(c("Variance Explained:",round(x.isi,3)), collapse=" ")))
    }
    if (max(x.isi) > x.isi.max) {
      # find the one that maximize the inclusive total sobol' index
      max.ind <- candidate[which.max(x.isi)]
      subset <- c(subset, max.ind)
      if(sparse) candidate <- candidate[x.isi>x.isi.max]
      candidate <- setdiff(candidate, max.ind)
      x.isi.max <- max(x.isi)
    } else {
      break
    }
  }
  
  imp <- rep(0,p)
  if (length(subset)==1) imp[subset] <- 1
  if (length(subset)>1) {
    # backward elimination
    if (verbose) print(noquote("Starting backward elimination..."))
    subset.imp <- nanne.be(X=X[,subset,drop=FALSE],y=y,n.knn=n.knn,n.mc=n.mc,mc.option=mc.option,
                           rescale=FALSE,parl=parl,verbose=verbose)
    imp[subset] <- subset.imp
  }
  
  return (imp)
}