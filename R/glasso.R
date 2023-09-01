glasso <- function(S, rho, maxit = 100, thr = 1e-4, zero = NULL, trace = FALSE, penalize.diagonal = FALSE) {
  L = rho
  #abs2 <- abs3 <- abs4 <- abs

  BIG=10e9
  
  if (!is.matrix(L)) {
    L = matrix(L, nrow(S), ncol(S))
  }

  if(!is.null(zero)){
    if(!is.matrix(zero)){ zero=matrix(zero,nrow=TRUE)}
    for(k in 1:nrow(zero)){
      i=zero[k,1]
      j=zero[k,2]
      L[i,j]=BIG
      L[j,i]=BIG
    }
  }
  
  shr <- sum(abs(S)) - sum(diag(S))
  n <- ncol(S)
  
  if (shr == 0) {
    W <- matrix(0, nrow(S), ncol(S))
    X <- matrix(0, nrow(S), ncol(S))
    diag(W) <- diag(L)
    for (i in seq_len(nrow(W))) {
      X[i,i] <- 1.0/max(W[i,i], .Machine$double.eps)
    }
  }
  shr = thr*shr/(n-1)
  
    
  thrLasso = max(shr/n, 2 * .Machine$double.eps)
  
  WXj <- rep(0, n)
  W <- S + 0
  X <- S * 0
  
  Wd <- diag(S) + diag(L)
  diag(W) <- Wd
  
  for (iter in seq_len(maxit)) {
    
    #print(iter)
    
    # dw <- 0
    # for (j in 1:n) {
    #   # WXj[] <- 0
    #   # for (i in 1:n) {
    #   #   if (X[i, j] != 0) {
    #   #     WXj = WXj + W[, i] * X[i, j]
    #   #   }
    #   # }
    #   ind <- which(X[, j] != 0)
    #   WXj <- W[, ind, drop = FALSE] %*% X[ind, j]
    #
    #   # repeat {
    #   #
    #   #   dlx = 0.0
    #   #   for (i in 1:n) {
    #   #     if (j != i) {
    #   #       a = S[i, j] - WXj[i] + Wd[i] * X[i, j]
    #   #       b = abs3(a) - L
    #   #       c = `if`(b > 0, sign(a) * b / Wd[i], 0)
    #   #       delta = c - X[i,j]
    #   #       if (delta != 0) {
    #   #         X[i,j] = c
    #   #         # WXj = WXj + W[, i] * delta
    #   #         update_WXj(WXj, W, i - 1L, delta)
    #   #         dlx = max(dlx, abs4(delta))
    #   #       }
    #   #     }
    #   #   }
    #   #
    #   #   if (dlx < thrLasso) break
    #   # }
    #   inner_loop(S, W, X, WXj, Wd, n, j - 1L, L, thrLasso)
    #
    #   WXj[j] = Wd[j]
    #   # dw = max(dw, sum(abs(WXj - W[, j])))
    #   dw = max(dw, abs_dist(WXj, W[, j]))
    #   W[, j] = W[j, ] <- WXj
    # }
    dw <- inner_loop(S, W, X, WXj, Wd, n, L, thrLasso)
    
    if (dw < shr) break
  }
  
  for (i in 1:n) {
    tmp <- 1 / drop(Wd[i] - crossprod(X[, i], W[, i]))
    X[, i] = -tmp * X[, i]
    X[i, i] = tmp
  }
  for (i in seq_len(n - 1)) {
    X[(i+1):n,i] = (X[(i+1):n,i] + X[i,(i+1):n])/2
    X[i,(i+1):n] = X[(i+1):n,i]
  }
  
  list(w = W, wi = X, niter = iter)
}
