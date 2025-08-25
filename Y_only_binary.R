library(pgdraw)
library(LaplacesDemon)
library(MASS)
briny_Y_only_1212 = function(Y,nrun,burn,thin,delta_n,a_theta_B,b_theta_B,theta_inf,start_adapt,Hmax,a,alpha0,alpha1){
  n = dim(Y)[1]
  q = dim(Y)[2]
  u = runif(nrun)
  H = Hmax + 1
  Hstar = Hmax
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
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
  ####
  H_hat = rep(NA,nrun)
  N_sample = ceiling((nrun - burn)/thin)
  gamma_hat = matrix(0,N_sample,q)
  #  Y_hat = matrix(0,n,q)
  Z_hat = list()
  B_hat = list()
  active_hat = list()
  sigma_hat = matrix(0,nrun,q)
  m = 1
  for (run in 1:nrun){
    theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
    D_Y = matrix(pgdraw(1,theta_Y),nrow = n,ncol = q)       
    for(i in 1:n){
      D_Yi = D_Y[i,]
      kappa_Yi = Y[i,] - 0.5      
      Sigma_Zi = chol2inv(chol(diag(H) +  B %*% diag(D_Yi,nrow = q) %*% t(B)))
      u_Zi = Sigma_Zi %*% ((B %*%(kappa_Yi - diag(D_Yi,nrow = q) %*% gamma_Y)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)       
    }
    #update gamma
    si_g1 = apply(D_Y, 2, sum) + 1/100
    sigma_gamma = diag(1 / si_g1)
    u_gamma_temp = Y - 0.5 - D_Y * (Z %*% B)
    u_gamma = sigma_gamma %*% (apply(u_gamma_temp, 2, sum))
    gamma_Y = mvrnorm(1, u_gamma, sigma_gamma)
    #update B
    # for(j in 1:q){
    #   D_Yj = D_Y[,j]
    #   sigma_Bj = chol2inv(chol(diag(theta_inv_B) + t(Z) %*% diag(D_Yj, nrow = n) %*% Z))
    #   u_Bj = sigma_Bj %*% t(Z) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
    #   B[,j] = mvrnorm(1, u_Bj, sigma_Bj)
    # }
    for(j in 1:q){
      D_Yj = D_Y[,j]
      if(j<H){
        sigma_Bj = chol2inv(chol(diag(theta_inv_B[1:j], nrow = j) + t(Z[,1:j]) %*% diag(D_Yj, nrow = n) %*% Z[,1:j]))
        u_Bj = sigma_Bj %*% t(Z[,1:j]) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
        B[1:j,j] = mvrnorm(1, u_Bj, sigma_Bj)
      }else{
        sigma_Bj = chol2inv(chol(diag(theta_inv_B, nrow = H) + t(Z) %*% diag(D_Yj, nrow = n) %*% Z))
        u_Bj = sigma_Bj %*% t(Z) %*% (Y[,j] - 0.5 - gamma_Y[j] * D_Yj)
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
        Z = cbind(Z,  rnorm(n))
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
      m = m + 1
    }
  }
  output = list("H" = H_hat,"Z" = Z_hat,"B" = B_hat,"active" = active_hat,
                "gamma" = gamma_hat)
  return(output)
}
# n = 400
# k = 4
# q = 50
# nrun = 15000
# burn = 5000
# thin = 5
# N_sample = ceiling((nrun - burn) / thin)
# gamma_ture = 0#rnorm(q)
# Z_ture = matrix(rnorm(n * k),nrow = n,ncol = k)
# B_ture = matrix(runif(k*q,0.25,1.25),nrow = k, ncol = q)
# logit_Y = Z_ture %*% B_ture + (matrix(1,n,1) %*% matrix(gamma_ture,1,q))
# prob_Y = plogis(logit_Y)
# U_Y = matrix(runif(n*q),n,q)
# Y = matrix(0,n,q)
# Y[which(U_Y < prob_Y)] = 1
# out = fix_H_B(Y, nrun, 4, 1 , 1, 2, 2,2,2)
# gamma_Y = rep(0,q)
# B = matrix(0,q,q)
# Z = matrix(0,n,n)
# j = 1
# for(i in burn:nrun){
#   if(i %% thin == 0){
#     gamma_Y = out$gamma[i,] / N_sample + gamma_Y
#     B = t(out$B[[i]]) %*% out$B[[i]] / N_sample + B
#     Z = out$Z[[i]] %*% t(out$Z[[i]]) / N_sample + Z
#     pY = out$Z[[i]]%*%out$B[[i]] 
#     j = j+1
#   }
# }
# print(mean((B - t(B_ture) %*% B_ture)**2))
# print(mean((Z - Z_ture %*% t(Z_ture))**2))
# plot(plogis(pY + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))),prob_Y)