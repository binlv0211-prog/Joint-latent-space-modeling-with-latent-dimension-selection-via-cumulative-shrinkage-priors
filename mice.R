library(mice) 
mice_prob = function(Y,M,s,t){
  n = dim(Y)[1]
  p = dim(Y)[2]
  Y_miss = Y
  Y_miss[M == 1] = 0
  prob_Y = matrix(0,n,p)
  prob_Y_hat = array(0,c(n,p,t-s))
  for (run in 1:t) {
    for (i in 1:p) {
      y_p = Y_miss[,i]
      y_cov = Y_miss[,-i]
      data_df = data.frame(y = y_p, X = y_cov)
      logit_model_df = glm(y ~ ., family = binomial(), data = data_df)
      prob_Y[,i] = predict(logit_model_df, type = "response")      
    }
    U_Y = matrix(runif(n*q),n,q)
    Y_hat = matrix(0,n,q)
    Y_hat[which(U_Y < prob_Y)] = 1   
    Y_miss[M == 1] = Y_hat[M == 1]
    if(run > s){
      prob_Y_hat[,,run-s] = prob_Y
    }
  }
  return(prob_Y_hat)
}

miss = function(Y,rate){
  n = dim(Y)[1]
  q = dim(Y)[2]
  M = matrix(0,n,q)
  U = matrix(runif(n * q),n,q)
  M[U < rate] = 1
  return(M)
}
# load("real_data1.RData")
# A = real_data$A
# Y = real_data$Y
# miss_rate = 0.05
# M = miss(Y,miss_rate)
# prob_Y_hat = mice_prob(Y,M,20,30)
# prob_Y = apply(prob_Y_hat, c(1,2), mean)
# library(pROC)
# print(auc(roc(Y[M==1],prob_Y [M==1])))