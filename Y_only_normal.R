library(pgdraw)
library(LaplacesDemon)
library(MASS)
normal_Y_only_1212 = function(Y,nrun,burn,thin,delta_n,a_sig,b_sig,a_theta_B,b_theta_B,
                              theta_inf,start_adapt,Hmax,a,alpha0,alpha1){
  #  set.seed(my_seed)
  n = dim(Y)[1]
  q = dim(Y)[2]
  u<-runif(nrun)
  #初始化
  H = Hmax + 1
  Hstar = Hmax
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  theta_inv = rep(1,H)
  theta_inv_B = rep(1,H)
  w = rep(1,H)
  zta = rep(1,H)
  gamma_Y = rep(0,q)
  B = matrix(rnorm(H * q),nrow = H, ncol = q)
  for (i in 2:H) {
    for (j in 1:(i-1)) {
      B[i,j] = 0
    }
  }
  inv_sigma = rep(1,q)
  #一些结果,后续可能会根据需求改变
  H_hat = rep(NA,nrun)
  N_sample = ceiling((nrun - burn)/thin)
  gamma_hat = matrix(0,N_sample,q)
  #  Y_hat = matrix(0,n,q)
  Z_hat = list()
  B_hat = list()
  active_hat = list()
  sigma_hat = matrix(0,N_sample,q)
  m = 1
  for (run in 1:nrun){
    for(i in 1:n){
      Sigma_Zi = chol2inv(chol(diag(H) + B %*% diag(inv_sigma,nrow = q) %*% t(B)))
      u_Zi = Sigma_Zi %*% ((B %*% diag(inv_sigma,nrow = q) %*%(Y[i,] - gamma_Y)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)
    }
    si_g1 = 1 / (n * inv_sigma + 1/100)
    sigma_gamma = diag(si_g1)##此处默认gamma先验的方差为n，后续可能会改动
    u_gamma_temp = Y - Z %*% B
    u_gamma = inv_sigma * si_g1 * (apply(u_gamma_temp, 2, sum))
    gamma_Y = mvrnorm(1, u_gamma, sigma_gamma)
    #update B
    for(j in 1:q){
      if(j<H){
        # print(dim(diag(theta_inv_B[1:j],nrow = j)))
        # print(dim(inv_sigma[j] * t(Z[,1:j]) %*% Z[,1:j]))
        sigma_Bj = chol2inv(chol(diag(theta_inv_B[1:j],nrow = j) + inv_sigma[j] * t(Z[,1:j]) %*% Z[,1:j]))
        u_Bj = inv_sigma[j] * sigma_Bj %*% t(Z[,1:j]) %*% (Y[,j] - gamma_Y[j])
        B[1:j,j] = mvrnorm(1, u_Bj, sigma_Bj)      
      }else{
        sigma_Bj = chol2inv(chol(diag(theta_inv_B,nrow = H) + inv_sigma[j] * t(Z) %*% Z))
        u_Bj = inv_sigma[j] * sigma_Bj %*% t(Z) %*% (Y[,j] - gamma_Y[j])
        B[,j] = mvrnorm(1, u_Bj, sigma_Bj)    
      }
    }
    #sample zta
    lhd_spike<-rep(0,H)
    lhd_slab<-rep(0,H)
    for(h in 1:H){
      lhd_spike[h] = exp(sum(log(dnorm(B[h,], mean = 0, sd = theta_inf^(1/2), log = FALSE))))
      lhd_slab[h] = dmvt(x = B[h,], mu=rep(0,q), S=(b_theta_B/a_theta_B)*diag(q), df=2*a_theta_B)
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
    #update inv_theta_B
    for (h in 1:H){
      if (zta[h] <= h){
        theta_inv_B[h] = theta_inf^(-1)
      }
      else{
        theta_inv_B[h] = rgamma(n=1,shape = a_theta_B + 0.5 * q,rate=b_theta_B + 0.5 * t(B[h,]) %*% B[h,])
      }
    }
    #update inv_sigma
    diff_Y = apply((Y - (matrix(1,n,1) %*% matrix(gamma_Y,1,q)) - (Z %*% B))**2, 2, sum)
    for (j in 1:q){
      inv_sigma[j] = rgamma(n=1,shape=a_sig + 0.5 * n,rate=b_sig + 0.5 * diff_Y[j])
    }	
    #update H[t]
    active = which(zta > c(1:H))
    Hstar = length(active)
    if (run >= start_adapt & u[run] <= exp(alpha0 + alpha1 * run)){
      if (Hstar < H - 1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H = Hstar + 1
        w = c(w[active],1-sum(w[active]))
        Z = cbind(Z[,active],rnorm(n))
        B = rbind(B[active,],c(rep(0,H-1),rnorm(q-H+1,0,theta_inf^(-1))))
        theta_inv_B = c(theta_inv_B[active],theta_inf^(-1))
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
        Z = cbind(Z,rnorm(n))
        B = rbind(B,c(rep(0,H-1),rnorm(q-H+1,0,theta_inf^(-1))))
        theta_inv_B = c(theta_inv_B,theta_inf^(-1))
        zta = c(zta,H-1)
      }
    }
    H_hat[run] = Hstar
    if((run > burn) &((run-burn) %% thin == 0)){
      gamma_hat[m,] = gamma_Y
      B_hat[[m]] = B
      Z_hat[[m]] = Z
      active_hat[[m]] = which(zta > c(1:H))
      sigma_hat[m,] = 1 / inv_sigma
      m = m + 1
    }
  }
  output = list("H" = H_hat,"Z" = Z_hat,"B" = B_hat,"active" = active_hat,
                "gamma" = gamma_hat,"sigma" = sigma_hat)
  return(output)
}
