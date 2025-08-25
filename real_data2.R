library(pgdraw)
library(LaplacesDemon)
library(MASS)
library(pROC)
library(foreach)
library(doParallel)
source("C:/Users/33701/Desktop/code/mice.R")
source("C:/Users/33701/Desktop/code/network_miss.R")
source("C:/Users/33701/Desktop/code/network_binary_miss.R")
source("C:/Users/33701/Desktop/code/data_process.R")
load("C:/Users/33701/Desktop/毕业论文/实际数据/real_data.RData8")
A = real_data$A
Y = real_data$Y
n = dim(Y)[1]
q = dim(Y)[2]
nrun = 15
burn = 10
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
Hmax = 15
alpha0 = -1
alpha1 = -5*10^(-4)
alpha_l = -0.5
alpha_u = 0.5#保证网络的密度大于0.5
B_l = 0.5
B_u = 1.5
K = 5
pll = 1
replation = 30
miss_rate = 0.05

t = 15
s = 10


cl = makeCluster(pll)      
registerDoParallel(cl)

res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon"), .combine = list, .multicombine = TRUE) %dopar% {
  M = miss(Y,miss_rate)
  prob_Y_mice_hat = mice_prob(Y,M,s,t)
  prob_Y_mice = apply(prob_Y_mice_hat, c(1,2), mean)
  
  net_out = network_only_1123(A,1234,nrun,burn,thin,delta_n,a_theta,b_theta,
                              theta_inf,start_adapt,Hmax,a,alpha0,alpha1)
  Z_post_hat = get_Z_post(net_out)
  Y_out = miss_Y_binary(Y,M,Z_post_hat,t,s)
  
  prob_Y_net_hat = array(0,dim = c(n,q,t-s))
  for (i in 1:(t-s)) {
    prob_Y_net_hat[,,i] = plogis(Z_post_hat %*% Y_out$B[i,,] + (matrix(1,n,1) %*% matrix(Y_out$gamma[i,],1,q)))
  }
  prob_Y_net = apply(prob_Y_net_hat, c(1,2), mean)
  
  net_binary_out = network_briny_miss(A,Y,M,nrun,burn,thin,delta_n,a_theta,b_theta,
                                      theta_inf,start_adapt,Hmax,a,alpha0,alpha1)
  N_sample = ceiling((nrun - burn)/thin)
  prob_Y_joint_hat = array(0,dim = c(n,q,N_sample))
  for (i in 1:N_sample) {
    prob_Y_joint_hat[,,i] = plogis(net_binary_out$Z[[i]][,net_binary_out$active[[i]]] %*% net_binary_out$B[[i]][net_binary_out$active[[i]],] + (matrix(1,n,1) %*% matrix(net_binary_out$gamma[i,],1,q)))
  }
  prob_Y_joint = apply(prob_Y_joint_hat, c(1,2), mean)
  re = list("net" = prob_Y_net, "joint" = prob_Y_joint, "mice" = prob_Y_mice,"Miss" = M, "net_out" = net_out, "joint_out" = net_binary_out)
  return(re)
}
stopCluster(cl)
