library(pgdraw)
library(LaplacesDemon)
library(MASS)
network_normal_kf = function(A,Y,H,a_sig,b_sig,nrun,burn,thin){
  n = dim(A)[1]
  q = dim(Y)[2]
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  theta_inv = rep(1,H)
  gamma_Y = rep(1,q)
  inv_sigma = rep(1,q)
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
  sigma_hat = matrix(0,N_sample,q)
  B_hat = array(0,dim = c(N_sample,H,q))
  Z_hat = array(0,dim = c(N_sample,n,H))
  m = 1
  for (run in 1:nrun){
    theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
    D_A = matrix(pgdraw(1, theta_A), nrow = n,ncol = n)
    D_A = (D_A + t(D_A)) / 2
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
      Sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i + B %*% diag(inv_sigma,nrow = q) %*% t(B)))
      u_Zi = Sigma_Zi %*% ((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) + (B %*% diag(inv_sigma,nrow = q) %*%(Y[i,] - gamma_Y)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)
    }
    si_g = 1 / (n * inv_sigma + 1/100)
    sigma_gamma = diag(si_g)##此处默认gamma先验的方差为100，后续可能会改动
    #u_gamma_temp = Y - Z %*% B
    u_gamma = inv_sigma * si_g * (apply(Y - Z %*% B, 2, sum))
    gamma_Y = mvrnorm(1, u_gamma, sigma_gamma)
    #print(gamma_Y)
    #update B 一个下三角矩阵的估计
    for(j in 1:q){
      if(j<H){
        sigma_Bj = chol2inv(chol(diag(j) + inv_sigma[j] * t(Z[,1:j]) %*% Z[,1:j]))
        u_Bj = inv_sigma[j] * sigma_Bj %*% t(Z[,1:j]) %*% (Y[,j] - gamma_Y[j])
        B[1:j,j] = mvrnorm(1, u_Bj, sigma_Bj)      
      }else{
        sigma_Bj = chol2inv(chol(diag(H) + inv_sigma[j] * t(Z) %*% Z))
        u_Bj = inv_sigma[j] * sigma_Bj %*% t(Z) %*% (Y[,j] - gamma_Y[j])
        B[,j] = mvrnorm(1, u_Bj, sigma_Bj)    
      }
    }
    #update inv_theta_B
    # for(h in 1:H){
    #   theta_inv_B[h] = rgamma(1,a_theta_B + 0.5 * q, b_theta_B + 0.5 * t(B[h,]) %*% B[h,])
    # } 
    #update inv_sigma
    diff_Y = apply((Y - (matrix(1,n,1) %*% matrix(gamma_Y,1,q)) - (Z %*% B))**2, 2, sum)
    for (j in 1:q){
      inv_sigma[j] = rgamma(n=1,shape=a_sig + 0.5 * n,rate=b_sig + 0.5 * diff_Y[j])
    }
    if((run > burn) &((run-burn) %% thin == 0)){
      sigma_hat[m,] = 1 / inv_sigma
      gamma_hat[m,] = gamma_Y
      alpha_hat[m,] = alpha
      B_hat[m,,] = B
      Z_hat[m,,] = Z
      m = m+1
    }
  }
  output = list("alpha" = alpha_hat,"B" = B_hat,
                "gamma" = gamma_hat,"Z" = Z_hat,"sigma" = sigma_hat)
  return(output)  
}


network_normal_kf_getZ = function(A,Y,gamma_Y,inv_sigma,B,nrun,burn,thin){
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
    #update alpha
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
      Sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i + B %*% diag(inv_sigma,nrow = q) %*% t(B)))
      u_Zi = Sigma_Zi %*% ((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) + (B %*% diag(inv_sigma,nrow = q) %*%(Y[i,] - gamma_Y)))
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