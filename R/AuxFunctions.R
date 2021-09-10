#' Rediscover Internal Functions
#'
#' Internal functions used by EventPointer in the different
#' steps of the algorithm
#'
#' @keywords internal
#' @name InternalFunctions
#' @return Internal outputs
#'
#'
NULL


#' @rdname InternalFunctions
speedglm.wfit2 <- function (y, X, intercept = TRUE, weights = NULL, row.chunk = NULL,
                            family = gaussian(), start = NULL, etastart = NULL, mustart = NULL,
                            offset = NULL, acc = 1e-08, maxit = 25, k = 2, sparselim = 0.9,
                            camp = 0.01, eigendec = TRUE, tol.values = 1e-07, tol.vectors = 1e-07,
                            tol.solve = .Machine$double.eps, sparse = NULL, method = c("eigen",
                                                                                       "Cholesky", "qr"), trace = FALSE, ...)
{
  nobs <- NROW(y)
  nvar <- ncol(X)
  if (missing(y))
    stop("Argument y is missing")
  if (missing(X))
    stop("Argument X is missing")
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  if (is.null(weights))
    weights <- rep(1, nobs)
  col.names <- dimnames(X)[[2]]
  method <- match.arg(method)
  fam <- family$family
  link <- family$link
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  if (is.null(sparse))
    sparse <- is.sparse(X = X, sparselim, camp)
  if (is.null(start)) {
    if (is.null(mustart))
      eval(family$initialize)
    eta <- if (is.null(etastart))
      family$linkfun(mustart)
    else etastart
    mu <- mustart
    start <- rep(0, nvar)
  }
  else {
    eta <- offset + as.vector(if (nvar == 1)
      X * start
      else {
        if (sparse)
          X %*% start
        else tcrossprod(X, t(start))
      })
    mu <- linkinv(eta)
  }
  iter <- 0
  dev <- sum(dev.resids(y, mu, weights))
  tol <- 1
  if ((fam == "gaussian") & (link == "identity"))
    maxit <- 1
  C_Cdqrls <- getNativeSymbolInfo("Cdqrls", PACKAGE = getLoadedDLLs()$stats)
  while ((tol > acc) & (iter < maxit)) {
    iter <- iter + 1
    beta <- start
    dev0 <- dev
    varmu <- variance(mu)
    mu.eta.val <- mu.eta(eta)
    z <- (eta - offset) + (y - mu)/mu.eta.val
    W <- (weights * mu.eta.val * mu.eta.val)/varmu
    X1 <- sqrt(W) * X
    XTX <- crossprod(X1)
    XTz <- t(crossprod((W * z), X))
    if (iter == 1 & method != "qr") {
      variable <- colnames(X)
      ris <- if (eigendec)
        control(XTX, , tol.values, tol.vectors, , method)
      else list(rank = nvar, pivot = 1:nvar)
      ok <- ris$pivot[1:ris$rank]
      if (eigendec) {
        XTX <- ris$XTX
        X <- X[, ok]
        XTz <- XTz[ok]
        start <- start[ok]
      }
      beta <- start
    }
    if (method == "qr") {
      ris <- .Call(C_Cdqrls, XTX, XTz, tol.values, FALSE)
      start <- if (ris$rank < nvar)
        ris$coefficients[ris$pivot]
      else ris$coefficients
    }
    else {
      start <- solve(XTX, XTz, tol = tol.solve)
    }
    eta <-  drop(X %*% start)
    mu <- linkinv(eta <- eta + offset)
    dev <- sum(dev.resids(y, mu, weights))
    tol <- max(abs(dev0 - dev)/(abs(dev) + 0.1))
    if (trace)
      cat("iter", iter, "tol", tol, "\n")
  }
  wt <- sum(weights)
  wtdmu <- if (intercept)
    sum(weights * y)/wt
  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- ris$rank
  dfr <- nobs - rank - sum(weights == 0)
  aic.model <- aic(y, nobs, mu, weights, dev) + k * rank
  #ll.nuovo <- speedglm:::ll.speedglm(fam, aic.model, rank)
  res <- (y - mu)/mu.eta(eta)
  resdf <- n.ok - rank
  RSS <- sum(W * res * res)
  var_res <- RSS/dfr
  dispersion <- if (fam %in% c("poisson", "binomial"))
    1
  else var_res
  if (method == "qr") {
    coefficients <- start
    coefficients[coefficients == 0] = NA
    ok <- ris$pivot[1:rank]
  }
  else {
    coefficients <- rep(NA, nvar)
    start <- as(start, "numeric")
    coefficients[ok] <- start
  }
  names(coefficients) <- col.names
  rval <- list(coefficients = coefficients,
               iter = iter, tol = tol, family = family, link = link,
               df = dfr, XTX = XTX, dispersion = dispersion, ok = ok,
               rank = rank, RSS = RSS, method = method, aic = aic.model,
               sparse = sparse, deviance = dev, nulldf = nulldf, nulldev = nulldev,
               ngoodobs = n.ok, n = nobs, intercept = intercept, convergence = (!(tol >
                                                                                    acc)))
  class(rval) <- "speedglm"
  rval
}

#' @rdname InternalFunctions 
expand.grid_fast <- function(a,b){
  cbind(
    rep(a,each=length(b))
    ,
    # rep(b,length(a))
    b
  )
}


# taken from https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
#' @rdname InternalFunctions 
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}



