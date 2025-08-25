library(pgdraw)
library(LaplacesDemon)
library(MASS)
network_briny_kf = function(A,Y,H,nrun,burn,thin){
  n = dim(A)[1]
  q = dim(Y)[2]
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  theta_inv = rep(1,H)
  gamma_Y = rep(1,q)
  B = matrix(rnorm(H * q),nrow = H, ncol = q)
  if(H > 1){
    for (i in 2:H) {
      for (j in 1:(i-1)) {
        B[i,j] = 0
      }
    } 
  }
  N_sample = ceiling((nrun - burn)/thin)
  alpha_hat = matrix(0,N_sample,n)
  gamma_hat = matrix(0,N_sample,q)
  
  B_hat = array(0,dim = c(N_sample,H,q))
  Z_hat = array(0,dim = c(N_sample,n,H))
  m = 1
  for (run in 1:nrun){
    theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
    D_A = matrix(pgdraw(1, theta_A), nrow = n,ncol = n)
    D_A = (D_A + t(D_A)) / 2
    theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
    D_Y = matrix(pgdraw(1,theta_Y),nrow = n,ncol = q)
    #update alpha
    Z_temp = Z %*% t(Z)
    for(i in 1:n){
      sigma_alphai = 1 / (sum(D_A[i,]) - D_A[i,i] + 1/100)#此处默认alpha先验的方差为1，后续可能会改动
      u_temp = A[i,] - 0.5 - D_A[i,] * (alpha + Z_temp[i,])
      u_alphai = sigma_alphai * (sum(u_temp) - u_temp[i])
      alpha[i] = rnorm(1,u_alphai,sqrt(sigma_alphai))
    }
    #update Z
    for(i in 1:n){
      Z_i = Z[-i,]
      D_Ai = (D_A[i,])[-i]
      alp_cons = alpha[i] + alpha[-i]
      kappa_Ai = (A[i,])[-i] - 0.5
      D_Yi = D_Y[i,]
      kappa_Yi = Y[i,] - 0.5
      Sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i + B %*% diag(D_Yi,nrow = q) %*% t(B)))
      u_Zi = Sigma_Zi %*% ((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) + (B %*%(kappa_Yi - diag(D_Yi,nrow = q) %*% gamma_Y)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)      
    }
    #update gamma
    si_g1 = apply(D_Y, 2, sum) + 1/100
    sigma_gamma = diag(1 / si_g1)
    u_gamma_temp = Y - 0.5 - D_Y * (Z %*% B)
    u_gamma = sigma_gamma %*% (apply(u_gamma_temp, 2, sum))
    gamma_Y = mvrnorm(1, u_gamma, sigma_gamma)
    for(j in 1:q){
      D_Yj = D_Y[,j]
      if(j<H){
        sigma_Bj = chol2inv(chol(diag(j, nrow = j) + t(Z[,1:j]) %*% diag(D_Yj, nrow = n) %*% Z[,1:j]))
        u_Bj = sigma_Bj %*% t(Z[,1:j]) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
        B[1:j,j] = mvrnorm(1, u_Bj, sigma_Bj)
      }else{
        sigma_Bj = chol2inv(chol(diag(H, nrow = H) + t(Z) %*% diag(D_Yj, nrow = n) %*% Z))
        u_Bj = sigma_Bj %*% t(Z) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
        B[,j] = mvrnorm(1, u_Bj, sigma_Bj)
      }
    }
    if((run > burn) &((run-burn) %% thin == 0)){
      gamma_hat[m,] = gamma_Y
      alpha_hat[m,] = alpha
      B_hat[m,,] = B
      Z_hat[m,,] = Z
      m = m+1
    }
  }
  output = list("alpha" = alpha_hat,"B" = B_hat,
                "gamma" = gamma_hat,"Z" = Z_hat)
  return(output)
}


network_briny_kf_getZ = function(A,Y,gamma_Y,B,nrun,burn,thin){
  n = dim(A)[1]
  q = dim(Y)[2]
  H = dim(B)[1]
  theta_inv = rep(1/100,H)
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  N_sample = ceiling((nrun - burn)/thin)
  alpha_hat = matrix(0,N_sample,n)
  Z_hat = array(0,dim = c(N_sample,n,H))
  m = 1
  for (run in 1:nrun){
    theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
    D_A = matrix(pgdraw(1, theta_A), nrow = n,ncol = n)
    D_A = (D_A + t(D_A)) / 2
    theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
    D_Y = matrix(pgdraw(1,theta_Y),nrow = n,ncol = q)
    Z_temp = Z %*% t(Z)
    for(i in 1:n){
      sigma_alphai = 1 / (sum(D_A[i,]) - D_A[i,i] + 1/100)#此处默认alpha先验的方差为1，后续可能会改动
      u_temp = A[i,] - 0.5 - D_A[i,] * (alpha + Z_temp[i,])
      u_alphai = sigma_alphai * (sum(u_temp) - u_temp[i])
      alpha[i] = rnorm(1,u_alphai,sqrt(sigma_alphai))
    }
    for(i in 1:n){
      Z_i = Z[-i,]
      D_Ai = (D_A[i,])[-i]
      alp_cons = alpha[i] + alpha[-i]
      kappa_Ai = (A[i,])[-i] - 0.5
      D_Yi = D_Y[i,]
      kappa_Yi = Y[i,] - 0.5
      Sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i + B %*% diag(D_Yi,nrow = q) %*% t(B)))
      u_Zi = Sigma_Zi %*% ((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) + (B %*%(kappa_Yi - diag(D_Yi,nrow = q) %*% gamma_Y)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)      
    }
    if((run > burn) &((run-burn) %% thin == 0)){
      alpha_hat[m,] = alpha
      Z_hat[m,,] = Z
      m = m+1
    }
  }
  output = list("alpha" = alpha_hat,"Z" = Z_hat)
  return(output)
}