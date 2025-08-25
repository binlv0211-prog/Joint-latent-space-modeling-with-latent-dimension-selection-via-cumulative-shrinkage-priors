library(pgdraw)
library(LaplacesDemon)
library(MASS)
network_briny_miss = function(A,Y,M,nrun,burn,thin,delta_n,a_theta,b_theta,
                              theta_inf,start_adapt,Hmax,a,alpha0,alpha1){
  #set.seed(my_seed)
  n = dim(A)[1]
  q = dim(Y)[2]
  u<-runif(nrun)
  #初始化
  H = Hmax + 1
  Hstar = Hmax
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  theta_inv = rep(1,H)
  w = rep(1,H)
  zta = rep(1,H)
  gamma_Y = rep(1,q)
  B = matrix(rnorm(H * q),nrow = H, ncol = q)
  # for (i in 2:H) {
  #   for (j in 1:(i-1)) {
  #     B[i,j] = 0
  #   }
  # }
  logit_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  prob_Y = plogis(logit_Y)
  U_Y = matrix(runif(n*q),n,q)
  Y_hat = matrix(0,n,q)
  Y_hat[which(U_Y < prob_Y)] = 1
  Y[M == 1] = Y_hat[M == 1]
  #一些结果,后续可能会根据需求改变
  N_sample = ceiling((nrun - burn)/thin)
  H_hat = rep(NA,N_sample)
  alpha_hat = matrix(0,N_sample,n)
  gamma_hat = matrix(0,N_sample,q)
  B_hat = list()
  Z_hat = list()
  active_hat = list()
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
    # sample zta
    lhd_spike<-rep(0,H)
    lhd_slab<-rep(0,H)
    for(h in 1:H){
      lhd_spike[h] = exp(sum(log(dnorm(Z[,h], mean = 0, sd = theta_inf^(1/2), log = FALSE))))
      #lhd_spike[h] = min(G, lhd_spike[h])
      lhd_slab[h] = dmvt(x = Z[,h], mu=rep(0,n), S=(b_theta/a_theta)*diag(n), df=2*a_theta)
      prob_h = w*c(rep(lhd_spike[h],h),rep(lhd_slab[h],H - h))
      if (sum(prob_h) == 0){
        prob_h = c(rep(0,H-1),1)
      }
      else{
        prob_h = prob_h/sum(prob_h)
      }
      zta[h] = c(1:H)%*%rmultinom(n=1, size=1, prob=prob_h)
    }
    #sample v and update w
    v = rep(NA,H)
    for (h in 1:(H - 1)){
      if (h == 1){
        v[h] = rbeta(1, shape1 = (Hmax + 1)**(delta_n) + sum(zta == h), shape2 = 1 + sum(zta > h))
      }else{
        v[h] = rbeta(1, shape1 = a + sum(zta == h), shape2 = 1 + sum(zta > h))
      }
    }
    v[H] = 1
    w[1] = v[1]
    for (h in 2:H){
      w[h] = v[h]*prod(1-v[1:(h-1)])  
    }
    # 6) sample theta^{-1}
    for (h in 1:H){
      if (zta[h] <= h){
        theta_inv[h] = theta_inf^(-1)
      }
      else{
        theta_inv[h] = rgamma(n=1,shape = a_theta + 0.5 * n,rate=b_theta + 0.5 * t(Z[,h]) %*% Z[,h])
      }
    }
    #update gamma
    si_g1 = apply(D_Y, 2, sum) + 1/100
    sigma_gamma = diag(1 / si_g1)
    u_gamma_temp = Y - 0.5 - D_Y * (Z %*% B)
    u_gamma = sigma_gamma %*% (apply(u_gamma_temp, 2, sum))
    gamma_Y = mvrnorm(1, u_gamma, sigma_gamma)
    #update B
    for(j in 1:q){
      D_Yj = D_Y[,j]
      sigma_Bj = chol2inv(chol(diag(H) + t(Z) %*% diag(D_Yj, nrow = n) %*% Z))
      u_Bj = sigma_Bj %*% t(Z) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
      B[,j] = mvrnorm(1, u_Bj, sigma_Bj)
    }
    # for(j in 1:q){
    #   D_Yj = D_Y[,j]
    #   if(j<H){
    #     sigma_Bj = chol2inv(chol(diag(j, nrow = j) + t(Z[,1:j]) %*% diag(D_Yj, nrow = n) %*% Z[,1:j]))
    #     u_Bj = sigma_Bj %*% t(Z[,1:j]) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
    #     B[1:j,j] = mvrnorm(1, u_Bj, sigma_Bj)
    #   }else{
    #     sigma_Bj = chol2inv(chol(diag(H, nrow = H) + t(Z) %*% diag(D_Yj, nrow = n) %*% Z))
    #     u_Bj = sigma_Bj %*% t(Z) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
    #     B[,j] = mvrnorm(1, u_Bj, sigma_Bj)
    #   }
    # }
    #update H[t]
    active = which(zta > c(1:H))
    Hstar = length(active)
    if (run >= start_adapt & u[run] <= exp(alpha0 + alpha1 * run)){
      if (Hstar < H - 1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H = Hstar + 1
        theta_inv = c(theta_inv[active],theta_inf^(-1))
        w = c(w[active],1-sum(w[active]))
        Z = cbind(Z[,active],rnorm(n,mean=0,sd=sqrt(theta_inf)))
        B = rbind(B[active,],rnorm(q))
        zta = c(zta[active],H-1)
      } else if (H < Hmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        H = H + 1
        v[H - 1] = rbeta(1,shape1=a,shape2=1)
        v = c(v,1)
        w = rep(NA,H)
        w[1] = v[1]
        for (h in 2:H){
          w[h] = v[h]*prod(1-v[1:(h-1)])
        }
        theta_inv = c(theta_inv,theta_inf^(-1))
        Z = cbind(Z,rnorm(n,mean=0,sd=sqrt(theta_inf)))
        B = rbind(B,rnorm(q))
        zta = c(zta,H-1)
      }
    }
    # 对缺失值进行插补
    logit_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
    prob_Y = plogis(logit_Y)
    U_Y = matrix(runif(n*q),n,q)
    Y_hat = matrix(0,n,q)
    Y_hat[which(U_Y < prob_Y)] = 1
    Y[M == 1] = Y_hat[M == 1]
    
    
    if((run > burn) &((run-burn) %% thin == 0)){
      H_hat[m] = Hstar
      gamma_hat[m,] = gamma_Y
      B_hat[[m]] = B
      Z_hat[[m]] = Z
      active_hat[[m]] = which(zta > c(1:H))
      #sigma_hat[m,] = 1 / inv_sigma
      alpha_hat[m,] = alpha
      m = m+1
    }  
  }
  output = list("H" = H_hat,"active" = active_hat,"alpha" = alpha_hat,"B" = B_hat,
                "gamma" = gamma_hat,"Z" = Z_hat)
  return(output)
}
