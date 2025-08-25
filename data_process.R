get_Z = function(n,k){
  return(matrix(rnorm(n*k),n,k))
}

get_alpha = function(n,alpha_l,alpha_u){
  return(runif(n,alpha_l,alpha_u))
}

get_A = function(alpha,Z){
  n = length(alpha)
  logit_A = Z %*% t(Z)+matrix(1,n,1)%*%matrix(alpha,1,n)+matrix(alpha,n,1)%*%matrix(1,1,n)
  prob_A = plogis(logit_A)
  U_A = matrix(runif(n*n),n,n)
  A = matrix(0,n,n)
  A[which(U_A < prob_A)] = 1
  for (i in 1:n){
    for (j in 1:(i - 1)) {
      A[j,i] = A[i,j]
    }
    A[i,i] = 1
  }
  return(A)
}

get_B = function(q,k,l,u,sample = F){
  if(sample){
    B = matrix(0,k,q)
    len = floor(q / k)
    ll = len * (k-1)
    for(i in 1:(k-1)){
      B[i,((i-1)*len+1):(i*len)] = runif(len,l,u)
    }
    B[k,(ll+1):q] = runif(q - ll,l,u)
    return(B)
  }
  else{
    return(matrix(runif(q*K),k,k))
  }
}

get_gamma = function(q){
  return(rnorm(q))
}

get_Y = function(gamma_Y,Z,B,continous = T){
  n = dim(Z)[1]
  q = dim(B)[2]
  if(continous){
    return(Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q)) + matrix(rnorm(n * q),n,q))
  }
  else{
    logit_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
    prob_Y = plogis(logit_Y)
    U_Y = matrix(runif(n*q),n,q)
    Y = matrix(0,n,q)
    Y[which(U_Y < prob_Y)] = 1
    return(Y)
  }
}

getmode <- function(v) {
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

get_AIC_binary = function(A,Y,alpha,gamma_Y,B,Z){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z)[2]
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - log(1 + exp(theta_Y))
  return(2 * (k+1)*(n+q) - 2 * (sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y)))
}

get_AIC_normal = function(A,Y,alpha,gamma_Y,B,Z){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z)[2]
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - theta_Y **2 / 2
  return(2 * (k+1)*(n+q) - 2 * (sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y)))
}

get_BIC_binary = function(A,Y,alpha,gamma_Y,B,Z){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z)[2]
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - log(1 + exp(theta_Y))
  return((k+1)*(n+q) * log(n) - 2 * (sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y)))
}

get_BIC_normal = function(A,Y,alpha,gamma_Y,B,Z){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z)[2]
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - theta_Y **2 / 2
  return((k+1)*(n+q) * log(n) - 2 * (sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y)))
}

get_DIC_binary = function(A,Y,alpha_mcmc, gamma_Y_mcmc,B_mcmc,Z_mcmc){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z_mcmc)[3]
  len = dim(alpha_mcmc)[1]
  alpha = apply(alpha_mcmc, 2, mean)
  gamma_Y = apply(gamma_Y_mcmc, 2, mean)
  B = apply(B_mcmc,c(2,3),mean)
  Z = apply(Z_mcmc,c(2,3),mean)
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - log(1 + exp(theta_Y))
  D_hat = - 2 * (sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y))
  bar_D = 1:len
  for (i in 1:len) {
    alphai = alpha_mcmc[i,]
    gamma_Yi = gamma_Y_mcmc[i,]
    Zi = matrix(Z_mcmc[i, , ], nrow = n, ncol = k) 
    Bi = matrix(B_mcmc[i, , ], nrow = k, ncol = q) 
    theta_Ai = Zi %*% t(Zi) + matrix(1,n,1)%*%matrix(alphai,1,n) + matrix(alphai,n,1) %*%matrix(1,1,n)
    theta_Yi = Zi %*% Bi + (matrix(1,n,1) %*% matrix(gamma_Yi,1,q))
    L_Ai = A * theta_Ai - log(1 + exp(theta_Ai))
    L_Yi = Y * theta_Yi - log(1 + exp(theta_Yi))
    bar_D[i] = - 2 * (sum(L_Ai[upper.tri(L_Ai, diag = FALSE)]) + sum(L_Yi))
  }
  return(2 * mean(bar_D) - D_hat)
}


get_DIC_normal = function(A,Y,alpha_mcmc, gamma_Y_mcmc,B_mcmc,Z_mcmc){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z_mcmc)[3]
  len = dim(alpha_mcmc)[1]
  alpha = apply(alpha_mcmc, 2, mean)
  gamma_Y = apply(gamma_Y_mcmc, 2, mean)
  B = apply(B_mcmc,c(2,3),mean)
  Z = apply(Z_mcmc,c(2,3),mean)
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - theta_Y **2 / 2
  D_hat = - 2 * (sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y))
  bar_D = 1:len
  for (i in 1:len) {
    alphai = alpha_mcmc[i,]
    gamma_Yi = gamma_Y_mcmc[i,]
    Zi = matrix(Z_mcmc[i, , ], nrow = n, ncol = k) 
    Bi = matrix(B_mcmc[i, , ], nrow = k, ncol = q) 
    theta_Ai = Zi %*% t(Zi) + matrix(1,n,1)%*%matrix(alphai,1,n) + matrix(alphai,n,1) %*%matrix(1,1,n)
    theta_Yi = Zi %*% Bi + (matrix(1,n,1) %*% matrix(gamma_Yi,1,q))
    L_Ai = A * theta_Ai - log(1 + exp(theta_Ai))
    L_Yi = Y * theta_Yi - theta_Yi **2 / 2
    bar_D[i] = - 2 * (sum(L_Ai[upper.tri(L_Ai, diag = FALSE)]) + sum(L_Yi))
  }
  return(2 * mean(bar_D) - D_hat)
}

get_WAIC_binary = function(A,Y,alpha_mcmc, gamma_Y_mcmc,B_mcmc,Z_mcmc){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z_mcmc)[3]
  len = dim(alpha_mcmc)[1]
  L_A = array(0,c(len,n,n))
  L_Y = array(0,c(len,n,q))
  for (i in 1:len) {
    alphai = alpha_mcmc[i,]
    gamma_Yi = gamma_Y_mcmc[i,]
    Zi = matrix(Z_mcmc[i, , ], nrow = n, ncol = k) 
    Bi = matrix(B_mcmc[i, , ], nrow = k, ncol = q) 
    theta_Ai = Zi %*% t(Zi) + matrix(1,n,1)%*%matrix(alphai,1,n) + matrix(alphai,n,1) %*%matrix(1,1,n)
    theta_Yi = Zi %*% Bi + (matrix(1,n,1) %*% matrix(gamma_Yi,1,q))
    L_A[i,,] = A * theta_Ai - log(1 + exp(theta_Ai))
    L_Y[i,,] = Y * theta_Yi - log(1 + exp(theta_Yi))
  }
  lppd_A = log(apply(exp(L_A), c(2,3),mean))
  lppd_Y = log(apply(exp(L_Y), c(2,3),mean))
  pwaic_A = apply(L_A, c(2,3),var)
  pwaic_Y = apply(L_Y, c(2,3),var)
  return(-2 * (sum(lppd_A[upper.tri(lppd_A, diag = FALSE)]) - sum(pwaic_A[upper.tri(pwaic_A, diag = FALSE)]) + sum(lppd_Y) - sum(pwaic_Y)))
}

get_WAIC_normal = function(A,Y,alpha_mcmc, gamma_Y_mcmc,B_mcmc,Z_mcmc){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z_mcmc)[3]
  len = dim(alpha_mcmc)[1]
  L_A = array(0,c(len,n,n))
  L_Y = array(0,c(len,n,q))
  for (i in 1:len) {
    alphai = alpha_mcmc[i,]
    gamma_Yi = gamma_Y_mcmc[i,]
    Zi = matrix(Z_mcmc[i, , ], nrow = n, ncol = k) 
    Bi = matrix(B_mcmc[i, , ], nrow = k, ncol = q) 
    theta_Ai = Zi %*% t(Zi) + matrix(1,n,1)%*%matrix(alphai,1,n) + matrix(alphai,n,1) %*%matrix(1,1,n)
    theta_Yi = Zi %*% Bi + (matrix(1,n,1) %*% matrix(gamma_Yi,1,q))
    L_A[i,,] = A * theta_Ai - log(1 + exp(theta_Ai))
    L_Y[i,,] = Y * theta_Yi - theta_Yi **2 / 2
  }
  lppd_A = log(apply(exp(L_A), c(2,3),mean))
  lppd_Y = log(apply(exp(L_Y), c(2,3),mean))
  pwaic_A = apply(L_A, c(2,3),var)
  pwaic_Y = apply(L_Y, c(2,3),var)
  return(-2 * (sum(lppd_A[upper.tri(lppd_A, diag = FALSE)]) - sum(pwaic_A[upper.tri(pwaic_A, diag = FALSE)]) + sum(lppd_Y) - sum(pwaic_Y)))
}

create_folds <- function(n, k = 5) {
  shuffled_indices <- sample(1:n)
  fold_indices <- split(shuffled_indices, ceiling(seq_along(shuffled_indices) / (length(shuffled_indices) / k)))
  folds_list <- list()
  for (i in 1:k) {
    test_indices <- fold_indices[[i]]
    train_indices <- unlist(fold_indices[-i])
    folds_list[[i]] <- list(
      train = train_indices,
      test = test_indices
    )
  }
  return(folds_list)
}

get_loglikehood_binary = function(A,Y,alpha,gamma_Y,B,Z){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z)[2]
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - log(1 + exp(theta_Y))
  return(sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y))
}

get_loglikehood_normal = function(A,Y,alpha,gamma_Y,B,Z){
  n = dim(A)[1]
  q = dim(Y)[2]
  k = dim(Z)[2]
  theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_A = A * theta_A - log(1 + exp(theta_A))
  L_Y = Y * theta_Y - theta_Y **2 / 2
  return(sum(L_A[upper.tri(L_A, diag = FALSE)]) + sum(L_Y))
}

get_rate = function(H_hat,H){
  re = dim(H_hat)[1]
  q = dim(H_hat)[2]
  H_count = matrix(0,re,q)
  H_count[which(H_hat == H)] = 1
  rate = apply(H_count, 2, mean)
  return(rate)
}
get_bias = function(H_hat,H){
  re = dim(H_hat)[1]
  q = dim(H_hat)[2]
  H_count = matrix(1,re,q)
  H_count[which(H_hat == H)] = 0
  bias = apply(abs(H_hat-H), 2, sum)/apply(H_count,2,sum)
  return(bias)
}