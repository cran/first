preproess.input <- function(X) {
  # one hot encoding of categorical features
  p <- ncol(X)
  cat.col.ind <- which(sapply(X,class)=="factor")
  num.col.ind <- setdiff(1:p, cat.col.ind)
  Xp <- as.matrix(X[,num.col.ind,drop=FALSE])
  for (ci in cat.col.ind)
    Xp <- cbind(Xp, stats::model.matrix(~.-1,X[,ci,drop=FALSE]))
  Xp <- as.matrix(Xp, nrow=nrow(Xp))
  return (Xp)
}

exp.var.knn <- function(
    X, 
    y, 
    subset, 
    duplicate = NULL,
    n.knn = 2,
    n.mc = nrow(X),
    twin.mc = FALSE,
    random.seed = NULL
) {
  
  # check for random seed
  if (!is.null(random.seed)) set.seed(random.seed)
  
  # arguments check for X and y
  if(!(inherits(X, "matrix") | inherits(X, "data.frame")))
    stop("X must be either a matrix or a data frame.")
  if (inherits(X, "matrix")) 
    X <- data.frame(X)
  if (inherits(y, "matrix"))
    y <- c(y) 
  if (nrow(X) != length(y))
    stop(sprintf("Size of X (%d) must match with size of y (%d)", nrow(X), length(y)))
  n <- nrow(X)
  p <- ncol(X)
  
  # argument check of using subsample
  n.mc <- as.integer(n.mc)
  n.mc <- if (n.mc > 0) min(n.mc, n) else n 
  twin.mc <- if (n %/% n.mc < 2) FALSE else twin.mc
  
  # argument check for subset
  if (length(subset) == 0) return (stats::var(y))
  
  # check if any row of X duplicate and preprocess categorical inputs
  if (is.null(duplicate)) {
    duplicate <- rep(FALSE, length(subset))
    for (i in 1:length(subset)) duplicate[i] <- (length(unique(X[,subset[i]])) < n)
  } else {
    if (length(duplicate) != p)
      stop(sprintf("Size of duplicate (%d) must be %d", length(duplicate), p))
    duplicate <- duplicate[subset]
  }
  X <- X[,subset,drop=FALSE]
  row.duplicated <- FALSE
  if (all(duplicate)) {
    if (length(subset) > 1) {
      if (nrow(unique(X)) < n)
        row.duplicated <- TRUE
    } else {
      row.duplicated <- TRUE
    }
  }
  X <- preproess.input(X)
  if (row.duplicated) {
    # add a column for random jittering
    X <- cbind(X, 1e-3*runif(n,min=-1,max=1))
  }
  
  # nearest-neighbor search and compute expected conditional variance
  if (n.mc < n) {
    if (twin.mc) {
      r <- n %/% n.mc
      keep.ind <- sample(1:n, n.mc*r)
      twin.ind <- twinning::twin(data=X[keep.ind,,drop=F], r=r, u1=1)
      query.ind <- keep.ind[twin.ind]
    } else {
      query.ind <- sample(1:n, n.mc, replace=FALSE)
    }
  } else {
    query.ind <- c(1:n)
  }
  nn.index <- FNN::knnx.index(data=X, query=X[query.ind,], k=n.knn)
  ev <- mean(apply(matrix(y[nn.index],ncol=n.knn),1,stats::var))
  
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
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#' @importFrom stats var model.matrix runif
#' @importFrom twinning twin
#' 
#' @param X a matrix or data frame for the factors / predictors.
#' @param y a vector for the responses.
#' @param noise a logical indicating whether the responses are noisy. 
#' @param n.knn number of nearest-neighbors for the inner loop conditional variance estimation. \code{n.knn=2} is recommended for regression, and \code{n.knn=3} for binary classification.
#' @param n.mc number of Monte Carlo samples for the outer loop expectation estimation.
#' @param twin.mc a logical indicating whether to use twinning subsamples, otherwise random subsamples are used. It is supported when the reduction ratio is at least 2. 
#' @param rescale a logical logical indicating whether to standardize the factors / predictors.
#' @param parl number of cores on which to parallelize the computation. If \code{NULL}, then no parallelization is done.
#' 
#' @details  
#' \code{totalsobol.knn} provides consistent estimation of the total Sobol' indices (Sobol, 1993) from scattered data.
#' When the output is clean/noiseless (\code{noise=FALSE}), \code{totalsobol.knn} implements the Nearest-Neighbor estimator proposed in Broto et al. (2020).
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
#' Two options are available for finding the subsamples: random and twinning (Vakayil and Joseph, 2022). 
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
#' Huang, C., & Joseph, V. R. (2024). Factor Importance Ranking and Selection using Total Indices. arXiv preprint arXiv:2401.00800.
#' 
#' Sobol', I. M. (2001). Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates. Mathematics and computers in simulation, 55(1-3), 271-280.
#' 
#' Broto, B., Bachoc, F., & Depecker, M. (2020). Variance reduction for estimation of Shapley effects and adaptation to unknown input distribution. SIAM/ASA Journal on Uncertainty Quantification, 8(2), 693-716.
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

totalsobol.knn <- function(
    X,
    y,
    noise,
    n.knn = NULL,
    n.mc = nrow(X),
    twin.mc = FALSE,
    rescale = TRUE, 
    parl = NULL
) {
  
  # arguments check for X and y
  if(!(inherits(X, "matrix") | inherits(X, "data.frame")))
    stop("X must be either a matrix of a data frame.")
  if (inherits(X, "matrix")) 
    X <- data.frame(X)
  if (inherits(y, "matrix"))
    y <- c(y) 
  if (nrow(X) != length(y))
    stop(sprintf("Size of X (%d) must match with size of y (%d)", nrow(X), length(y)))
  n <- nrow(X)
  p <- ncol(X)
  
  # argument check of using subsample
  n.mc <- as.integer(n.mc)
  n.mc <- if (n.mc > 0) min(n.mc, n) else n 
  twin.mc <- if (n %/% n.mc < 2) FALSE else twin.mc
  
  # preprocess for features
  col.class.error <- FALSE
  col.class.error.msg <- "\n"
  col.class <- sapply(X,class)
  for (i in 1:p) {
    if (!(col.class[i] %in% c("logical","integer","numeric","factor"))) {
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
      "Please convert above columns to the supported column types (logical/integer/numeric/factor).\n",
      sep="")
    stop(col.class.error.msg)
  }
  # find if any column has duplicate value
  duplicate <- rep(FALSE, p)
  for (i in 1:p) duplicate[i] <- (length(unique(X[,i])) < n)
  # rescale logical and continuous inputs
  if (rescale) {
    for (i in 1:p) {
      if (col.class[i]=="logical")
        X[,i] <- 2 * (X[,i] - 0.5)
      if (col.class[i] %in% c("integer","numeric"))
        X[,i] <- c(scale(X[,i]))
    }
  }
  
  # argument check for n.knn
  if (length(unique(y)) == 1) 
    stop("y must have more than one unique value.")
  if (is.null(n.knn)) {
    if (length(unique(y)) > 2) {
      n.knn <- 2
    } else {
      n.knn <- 3
    }
  }
  
  # compute total Sobol' indices
  if (noise) {
    noise.var <- exp.var.knn(
      X = X,
      y = y,
      subset = c(1:p), 
      duplicate = duplicate,
      n.knn = n.knn, 
      n.mc = n.mc,
      twin.mc = twin.mc
    )
  } else {
    noise.var <- 0
  }
  y.var <- max(stats::var(y) - noise.var, 0)
  if (y.var == 0) {
    tsi <- rep(0, p)
  } else {
    if (is.null(parl)) {
      xe.var <- sapply(1:p, function(i) exp.var.knn(
        X = X,
        y = y,
        subset = setdiff(1:p, i),
        n.knn = n.knn, 
        duplicate = duplicate,
        n.mc = n.mc,
        twin.mc = twin.mc
      ))
    } else {
      seeds <- sample(1:1e9, p)
      clus <- parallel::makeCluster(parl)
      parallel::clusterExport(clus, list("preproess.input","exp.var.knn"), envir=environment())
      xe.var <- parallel::parSapply(clus, 1:p, function(i) exp.var.knn(
        X = X,
        y = y,
        subset = setdiff(1:p, i),
        n.knn = n.knn, 
        duplicate = duplicate,
        n.mc = n.mc,
        twin.mc = twin.mc,
        random.seed = seeds[i]
      ))
      parallel::stopCluster(clus)
    }
    xe.var <- pmax(xe.var - noise.var, 0)
    tsi <- xe.var / y.var
  }
  
  return (tsi)
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
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#' @importFrom stats var model.matrix runif
#' @importFrom twinning twin
#' 
#' @param X a matrix or data frame for the factors / predictors.
#' @param y a vector for the responses.
#' @param n.knn number of nearest-neighbors for the inner loop conditional variance estimation. \code{n.knn=2} is recommended for regression, and \code{n.knn=3} for binary classification.
#' @param n.mc number of Monte Carlo samples for the outer loop expectation estimation. 
#' @param twin.mc a logical indicating whether to use twinning subsamples, otherwise random subsamples are used. It is supported when the reduction ratio is at least 2. 
#' @param rescale a logical logical indicating whether to standardize the factors / predictors.
#' @param n.forward number of times to run the forward selection phase to tradeoff between efficiency and accuracy. \code{n_forward=2} is recommended. To run the complete forward selection, please set \code{n_forward} to the number of factors / predictors. 
#' @param parl number of cores on which to parallelize the computation. If \code{NULL}, then no parallelization is done.
#' @param verbose a logical indicating whether to display intermediate results.
#' 
#' @details
#' \code{first} provides factor importance ranking and selection directly from scattered data without any model fitting. 
#' Factor importance is computed based on total Sobol' indices (Sobol', 2001), 
#' which is connected to the approximation error introduced by excluding the factor of interest (Huang and Joseph, 2024). 
#' The estimation procedure adapts from the Nearest-Neighbor estimator in Broto et al. (2020) to account for the noisy data. 
#' Integrating it with forward selection and backward elimination allows for factor selection.
#' 
#' \code{first} belongs to the class of forward-backward selection with early dropping algorithm (Borboudakis and Tsamardinos, 2019). 
#' In forward selection, each time we find the candidate that maximizes the output variance that can be explained. 
#' For candidates that cannot improve the variance explained conditional on the selected factors, they are removed from the candidate set. 
#' This forward selection step is run `n_forward` times to tradeoff between accuracy and efficiency. \code{n_forward=2} is recommended in Yu et al. (2020). 
#' To run the complete forward selection, please set `n_forward` to the number of factors / predictors. 
#' In backward elimination, we again remove one factor at a time, starting with the factor that can improve the explained variance most, till no factor can further improve. 
#' 
#' \code{n.knn=2} nearest-neighbors is recommended for integer/numeric output, and \code{n.knn=3} is suggested for binary output.
#' For numeric inputs, it is recommended to standardize them via setting the argument \code{rescale=TRUE}.
#' Categorical inputs are transformed via one-hot encoding for the nearest-neighbor search. 
#' To speed up the nearest-neighbor search, k-d tree from the \pkg{FNN} package is used. 
#' Also, parallel computation is also supported via the \pkg{parallel} package.
#' 
#' For large datasets, we support the use of subsamples for further acceleration.
#' Use argument \code{n.mc} to specify the number of subsamples. 
#' Two options are available for finding the subsamples: random and twinning (Vakayil and Joseph, 2022). 
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
#' Sobol', I. M. (2001). Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates. Mathematics and computers in simulation, 55(1-3), 271-280.
#' 
#' Broto, B., Bachoc, F., & Depecker, M. (2020). Variance reduction for estimation of Shapley effects and adaptation to unknown input distribution. SIAM/ASA Journal on Uncertainty Quantification, 8(2), 693-716.
#' 
#' Borboudakis, G., & Tsamardinos, I. (2019). Forward-backward selection with early dropping. The Journal of Machine Learning Research, 20(1), 276-314.
#' 
#' Yu, K., Guo, X., Liu, L., Li, J., Wang, H., Ling, Z., & Wu, X. (2020). Causality-based feature selection: Methods and evaluations. ACM Computing Surveys (CSUR), 53(5), 1-36.
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

first <- function(
    X,
    y,
    n.knn = NULL,
    n.mc = nrow(X),
    twin.mc = FALSE,
    rescale = TRUE,
    n.forward = 2,
    parl = NULL,
    verbose = FALSE
) {
  
  # arguments check for X and y
  if(!(inherits(X, "matrix") | inherits(X, "data.frame")))
    stop("X must be either a matrix of a data frame.")
  if (inherits(X, "matrix")) 
    X <- data.frame(X)
  if (inherits(y, "matrix"))
    y <- c(y) 
  if (nrow(X) != length(y))
    stop(sprintf("Size of X (%d) must match with size of y (%d)", nrow(X), length(y)))
  n <- nrow(X)
  p <- ncol(X)
  
  # argument check of using subsample
  n.mc <- as.integer(n.mc)
  n.mc <- if (n.mc > 0) min(n.mc, n) else n 
  twin.mc <- if (n %/% n.mc < 2) FALSE else twin.mc
  
  # preprocess for features
  col.class.error <- FALSE
  col.class.error.msg <- "\n"
  col.class <- sapply(X,class)
  for (i in 1:p) {
    if (!(col.class[i] %in% c("logical","integer","numeric","factor"))) {
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
      "Please convert above columns to the supported column types (logical/integer/numeric/factor).\n",
      sep="")
    stop(col.class.error.msg)
  }
  # find if any column has duplicate value
  duplicate <- rep(FALSE, p)
  for (i in 1:p) duplicate[i] <- (length(unique(X[,i])) < n)
  # rescale logical and continuous inputs
  if (rescale) {
    for (i in 1:p) {
      if (col.class[i]=="logical")
        X[,i] <- 2 * (X[,i] - 0.5)
      if (col.class[i] %in% c("integer","numeric"))
        X[,i] <- c(scale(X[,i]))
    }
  }
  
  # argument check for n.knn
  if (length(unique(y)) == 1) 
    stop("y must have more than one unique value.")
  if (is.null(n.knn)) {
    if (length(unique(y)) > 2) {
      n.knn <- 2
      if (verbose) print("y has more than 2 unique values, setting it to regression problem with suggested n.knn = 2.")
    } else {
      n.knn <- 3
      if (verbose) print("y has only 2 unique values, setting it to binary classification problem with suggested n.knn = 3.")
    }
  }
  
  # forward selection 
  if (verbose) print("Starting forward selection...")
  y.var <- stats::var(y)
  subset <- c()
  x.var.max <- 0
  for (l in 1:n.forward) {
    if (verbose) print(sprintf("Phase-%d forward selection...", l))
    none.added.to.subset <- TRUE
    candidate <- setdiff(1:p, subset)
    while (length(candidate) > 0) {
      # compute total Sobol' effect for -x (x for current subset)
      if (is.null(parl)) {
        nx.var <- sapply(candidate, function(i) exp.var.knn(
          X = X,
          y = y,
          subset = c(subset,i),
          duplicate = duplicate,
          n.knn = n.knn, 
          n.mc = n.mc,
          twin.mc = twin.mc
        ))
      } else {
        seeds <- sample(1:1e9, p)
        clus <- parallel::makeCluster(parl)
        parallel::clusterExport(clus, list("preproess.input","exp.var.knn"), envir=environment())
        nx.var <- parallel::parSapply(clus, candidate, function(i) exp.var.knn(
          X = X,
          y = y,
          subset = c(subset,i),
          duplicate = duplicate,
          n.knn = n.knn, 
          n.mc = n.mc,
          twin.mc = twin.mc,
          random.seed = seeds[i]
        ))
        parallel::stopCluster(clus)
      }
      x.var <- pmax(y.var - nx.var, 0)
      if (verbose) {
        print(paste(c("current selection:", subset), collapse=" "))
        print(sprintf("current variance explained: %.3f", x.var.max))
        print(paste(c("candidate to add:", paste(paste(candidate,round(x.var,3),sep="("),")",sep="")), collapse=" "))
      }
      if (max(x.var) > x.var.max) {
        # find the index to add such that the variance explained is maximized
        add.ind <- candidate[which.max(x.var)]
        candidate <- candidate[x.var > x.var.max]
        candidate <- setdiff(candidate, add.ind)
        subset <- c(subset, add.ind)
        none.added.to.subset <- FALSE
        x.var.max <- max(x.var)
        if (verbose) print(sprintf("add candidate %d (%.3f).", add.ind, x.var.max))
      } else {
        break
      }
    }
    if (none.added.to.subset) {
      if (verbose) print("early termination since none of the candidates can be added in this phase. ")
      break
    }
  }
  
  # backward elimination 
  if (verbose) print("Starting backward elimination...")
  subset <- sort(subset)
  while (length(subset) > 0) {
    # compute total sobol' effect for -(x/{i}) (x for current subset)
    if (is.null(parl)) {
      nx.var <- sapply(subset, function(i) exp.var.knn(
        X = X,
        y = y,
        subset = setdiff(subset, i),
        duplicate = duplicate,
        n.knn = n.knn, 
        n.mc = n.mc,
        twin.mc = twin.mc
      ))
    } else {
      seeds <- sample(1:1e9, p)
      clus <- parallel::makeCluster(parl)
      parallel::clusterExport(clus, list("preproess.input","exp.var.knn"), envir=environment())
      nx.var <- parallel::parSapply(clus, subset, function(i) exp.var.knn(
        X = X,
        y = y,
        subset = setdiff(subset, i),
        duplicate = duplicate,
        n.knn = n.knn, 
        n.mc = n.mc,
        twin.mc = twin.mc,
        random.seed = seeds[i]
      ))
      parallel::stopCluster(clus)
    }
    x.var <- pmax(y.var - nx.var, 0)
    if (verbose) {
      print(paste(c("current selection:", subset), collapse=" "))
      print(sprintf("current variance explained: %.3f", x.var.max))
      print(paste(c("candidate to remove:", paste(paste(subset,round(x.var,3),sep="("),")",sep="")), collapse=" "))
    }
    if (max(x.var) > x.var.max) {
      # find the index to remove such that the variance explained is maximized
      remove.ind <- subset[which.max(x.var)]
      subset <- setdiff(subset, remove.ind)
      x.var.max <- max(x.var)
      if (verbose) print(sprintf("remove candidate %d (%.3f).", remove.ind, x.var.max))
    } else {
      break
    }
  }
  
  # compute importance via total Sobol' indices
  imp <- rep(0, p)
  if (length(subset) > 0) {
    noise.var <- y.var - x.var.max
    imp[subset] <- (nx.var - noise.var) / x.var.max
  }
  
  return (imp)
  
}
