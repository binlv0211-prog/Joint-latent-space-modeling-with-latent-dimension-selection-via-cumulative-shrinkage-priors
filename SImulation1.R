library(pgdraw)
library(LaplacesDemon)
library(MASS)
library(pROC)
library(foreach)
library(doParallel)
source("C:/Users/33701/Desktop/code/network_only.R")
source("C:/Users/33701/Desktop/code/network_binary.R")
source("C:/Users/33701/Desktop/code/network_normal.R")
source("C:/Users/33701/Desktop/code/Y_only_binary.R")
source("C:/Users/33701/Desktop/code/Y_only_normal.R")
source("C:/Users/33701/Desktop/code/network_normal_kl.R")
source("C:/Users/33701/Desktop/code/network_binary_kl.R")
source("C:/Users/33701/Desktop/code/data_process.R")
nrun = 15000
burn = 10000
thin = 5
cn = c(50,100,150,300)
cq = c(10,20,30,60)
k = 3
a = 8
a_sig = 1
b_sig = 1
a_theta = 3
b_theta = 3
a_theta_B = 3
b_theta_B = 3
theta_inf = 0.10
start_adapt = 500
Hmax = 7
alpha0 = -1
alpha1 = -5*10^(-4)
alpha_l = -0.5
alpha_u = 0.5#保证网络的密度大于0.5
B_l = 0.25
B_u = 1.25
pll = 1
replation = 100
K = 5


cl = makeCluster(pll)      
registerDoParallel(cl)

res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon"), .combine = list, .multicombine = TRUE) %dopar% {
  out_net = list()
  out_net_Y = list()
  out_Y = list()
  data_ture = list()
  res_criterion_total = list()
  res_likehood_total = list()
  
  
  for (i in 1:length(cn)) {
    n = cn[i]
    q = cq[i]
    delta_n = 1 + 0.001 * n 
    alpha_ture = get_alpha(n,alpha_l,alpha_u)
    Z_ture = get_Z(n,k)
    A = get_A(alpha_ture,Z_ture)
    gamma_ture = get_gamma(q)
    B_ture = get_B(q,k,B_l,B_u,T)
    Y = get_Y(gamma_ture,Z_ture,B_ture,continous = T)
    data_ture[[i]] = list("alpha" = alpha_ture,"Z" = Z_ture,"gamma" = gamma_ture,"B" = B_ture,"A" = A,"Y" = Y)
    net_out = network_only_1123(A, 1314, nrun,burn,thin, delta_n, a_theta, b_theta, 
                                theta_inf, start_adapt, Hmax,a, alpha0,alpha1)
    net_Y_out = network_nomal(A, Y,nrun,burn,thin, delta_n,a_sig,b_sig
                              ,a_theta,b_theta, #a_theta_B,b_theta_B,
                              theta_inf,start_adapt, Hmax, a,alpha0,alpha1)
    Y_out = normal_Y_only_1212(Y,nrun,burn,thin,delta_n,a_sig,b_sig,
                               a_theta_B,b_theta_B,
                               theta_inf = 0.05,start_adapt,Hmax,a,alpha0,alpha1)
    out_net[[i]] = net_out
    out_net_Y[[i]] = net_Y_out
    out_Y[[i]] = Y_out
    
    res_criterion = list()
    for (j in 1:Hmax) {
      res_criterion[[j]] = network_normal_kf(A,Y,j,a_sig,b_sig,nrun,burn,thin)
    }
    res_criterion_total[[i]] = res_criterion
    loglikehood = matrix(0,Hmax,K)
    folds <- create_folds(n, k = K) 
    
    for (jj in 1:K) {
      train_ids <- folds[[jj]]$train
      test_ids <- folds[[jj]]$test
      n_train = length(train_ids)
      n_test = length(test_ids)
      A_train = A[train_ids,train_ids]
      Y_train = Y[train_ids,]
      A_test = A[test_ids,test_ids]
      Y_test = Y[test_ids,]
      for (j in 1:Hmax) {
        kf_out = network_normal_kf(A_train,Y_train,j,a_sig,b_sig,nrun,burn,thin) 
        inv_sigma =  1 / apply(kf_out$sigma, 2, mean)
        gamma_Y = apply(kf_out$gamma, 2, mean)
        B = apply(kf_out$B,c(2,3),mean)
        kf_outz = network_normal_kf_getZ(A_test,Y_test,gamma_Y,inv_sigma,B,nrun,burn,thin)
        alpha_test = apply(kf_outz$alpha,2,mean)
        Z_test = apply(kf_outz$Z,c(2,3),mean)
        loglikehood[j,jj] = get_loglikehood_binary(A_test,Y_test,alpha_test,gamma_Y,B,Z_test)      
      }
    }
    res_likehood_total[[i]] = loglikehood
  }
  re = list("out_net" = out_net,"out_net_Y" = out_net_Y,"out_Y" = out_Y,"data_ture" = data_ture, "res_criterion_total" = res_criterion_total,"res_likehood_total" = res_likehood_total)
  return(re)
}
stopCluster(cl)


cl = makeCluster(pll)      
registerDoParallel(cl)

res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon"), .combine = list, .multicombine = TRUE) %dopar% {
  out_net = list()
  out_net_Y = list()
  out_Y = list()
  data_ture = list()
  res_criterion_total = list()
  res_likehood_total = list()
  
  for (i in 1:length(cn)) {
    n = cn[i]
    q = cq[i]
    delta_n = 1 + 0.001 * n 
    alpha_ture = get_alpha(n,alpha_l,alpha_u)
    Z_ture = get_Z(n,k)
    A = get_A(alpha_ture,Z_ture)
    gamma_ture = get_gamma(q)
    B_ture = get_B(q,k,B_l,B_u,T)
    Y = get_Y(gamma_ture,Z_ture,B_ture,continous = F)
    data_ture[[i]] = list("alpha" = alpha_ture,"Z" = Z_ture,"gamma" = gamma_ture,"B" = B_ture,"A" = A,"Y" = Y)
    net_out = network_only_1123(A, 1313, nrun,burn,thin, delta_n, a_theta, b_theta, 
                                theta_inf, start_adapt, Hmax, a, alpha0,alpha1)
    net_Y_out = network_briny_1130(A,Y,nrun,burn,thin,delta_n,a_theta,b_theta,
                                   theta_inf,start_adapt,Hmax,a,alpha0,alpha1)
    Y_out = briny_Y_only_1212(Y,nrun,burn,thin, delta_n ,
                              a_theta_B,b_theta_B, theta_inf = 0.05,
                              start_adapt = 2 * nrun, Hmax, a,alpha0,alpha1)
    out_net[[i]] = net_out
    out_net_Y[[i]] = net_Y_out
    out_Y[[i]] = Y_out
    
    res_criterion = list()
    for (j in 1:Hmax) {
      res_criterion[[j]] = network_briny_kf(A,Y,j,nrun,burn,thin)
    }
    res_criterion_total[[i]] = res_criterion
    loglikehood = matrix(0,Hmax,K)
    folds <- create_folds(n, k = K) 
    
    for (jj in 1:K) {
      train_ids <- folds[[jj]]$train
      test_ids <- folds[[jj]]$test
      n_train = length(train_ids)
      n_test = length(test_ids)
      A_train = A[train_ids,train_ids]
      Y_train = Y[train_ids,]
      A_test = A[test_ids,test_ids]
      Y_test = Y[test_ids,]
      for (j in 1:Hmax) {
        kf_out = network_briny_kf(A_train,Y_train,j,nrun,burn,thin) 
        gamma_Y = apply(kf_out$gamma, 2, mean)
        B = apply(kf_out$B,c(2,3),mean)
        kf_outz = network_briny_kf_getZ(A_test,Y_test,gamma_Y,B,nrun,burn,thin)
        alpha_test = apply(kf_outz$alpha,2,mean)
        Z_test = apply(kf_outz$Z,c(2,3),mean)
        loglikehood[j,jj] = get_loglikehood_binary(A_test,Y_test,alpha_test,gamma_Y,B,Z_test)   
      }
    }
    res_likehood_total[[i]] = loglikehood
  }
  re = list("out_net" = out_net,"out_net_Y" = out_net_Y,"out_Y" = out_Y,"data_ture" = data_ture, "res_criterion_total" = res_criterion_total,"res_likehood_total" = res_likehood_total)
  return(re)
}
stopCluster(cl)