require(MASS)

#------------------------------#
#      Projection matrices     #
#------------------------------#

P <- function(A){
  return(A %*% ginv(t(A) %*% A) %*% t(A))
}

Q <- function(A) {
  n = nrow(A)
  return(diag(n) - P(A))
}

#------------------------------#
#       Response SIMPLS        #
#------------------------------#

resSIMPLS <- function(X, Y, u) {
  n = nrow(X)
  r = ncol(Y)
  res = Q(X) %*% Y
  S = t(res) %*% res / n
  U = t(Y) %*% P(X) %*% Y
  W_k = matrix(rep(0, r))
  for (k in 0:(u-1)) {
    Q_e = Q(S %*% W_k)
    eigensystem = eigen(Q_e %*% U %*% Q_e)
    w = Re(eigensystem$vectors[, 1])
    W_k = cbind(W_k, w)
  }
  Gamma = W_k[, 2:(u+1)]
  return(ginv(t(X) %*% X) %*% t(X) %*% Y %*% P(Gamma))
}

#------------------------------#
#  K-fold CV for response PLS  #
#------------------------------#
cv.resSIMPLS = function(X, Y, k = 5) {
  n = nrow(X)
  r = ncol(Y)
  index = sample(1:n, floor(n/k), replace = FALSE)
  X_test = X[index, ]
  Y_test = Y[index, ]
  X_train = X[-index, ]
  Y_train = Y[-index, ]
  err = numeric(r)
  for (i in 1:r) {
    tmp_beta = resSIMPLS(X_train, Y_train, i)
    hat_Y = X_test %*% tmp_beta
    err[i] = sum((hat_Y - Y_test)^2)
  }
  u = which(err == min(err))
  return(u)
}

