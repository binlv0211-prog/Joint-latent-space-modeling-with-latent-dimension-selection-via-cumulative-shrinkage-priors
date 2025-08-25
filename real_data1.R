library(pgdraw)
library(LaplacesDemon)
library(MASS)
source("C:/Users/33701/Desktop/code/network_only.R")
source("C:/Users/33701/Desktop/code/network_binary.R")
source("C:/Users/33701/Desktop/code/network_normal.R")
source("C:/Users/33701/Desktop/code/Y_only_binary.R")
source("C:/Users/33701/Desktop/code/Y_only_normal.R")
source("C:/Users/33701/Desktop/code/network_normal_kl.R")
source("C:/Users/33701/Desktop/code/network_binary_kl.R")
source("C:/Users/33701/Desktop/code/data_process.R")
load("C:/Users/33701/Desktop/毕业论文/实际数据/french.RData")
A = french_data$A
Y = french_data$Y
n = dim(Y)[1]
q = dim(Y)[2]
nrun = 15000
burn = 10000
thin = 5
a = 5
delta_n = 1 + 0.001 * n 
a_sig = 1
b_sig = 1
a_theta = 3
b_theta = 3
a_theta_B = 3
b_theta_B = 3
theta_inf = 0.10
start_adapt = 500
Hmax = 4
alpha0 = -1
alpha1 = -5*10^(-4)
alpha_l = -0.5
alpha_u = 0.5#保证网络的密度大于0.5
B_l = 0.5
B_u = 1.5
K = 5
replation = 100

res = list()
res_criterion_total = list()
res_likehood_total = list()
for (i in 1:replation) {
  res[[i]] = network_briny_1130(A,Y,nrun,burn,thin,delta_n,a_theta,b_theta,
                                theta_inf,start_adapt,Hmax,a,alpha0,alpha1)
}

for (i in 1:replation) {
  res_criterion = list()
  for (j in 1:Hmax) {
    res_criterion[[j]] = network_normal_kf(A,Y,j,a_sig,b_sig,nrun,burn,thin)
  }
  res_criterion_total[[i]] = res_criterion
}

for (i in 1:replation) {
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
