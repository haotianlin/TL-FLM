require(extraDistr)
require(gss)
require(pracma)
require(RandomFieldsUtils)
require(MASS)
require(truncnorm)
require(caTools)

source("data_generation.R")
source("Transfer_FLR_functions.R")
set.seed(2022)

grid.size = 50;  sigma = .5; K = 100 ;

grid = gauss.quad(size = grid.size, interval = c(0,1))
######### pick the covariance kernel for X ###############
C = function(s,t){
  Matern_kernel(x1 = s, x2 = t, kernel = "Exp", rho = 1/15)
}

r1 = 1; # eigenvalue decay rate of RKHS kernel
r2 = 1; # eigenvalue decay rate of covariance kernel of X
rate = 2*(r1+r2)/( 2*(r1+r2) + 1 )

######### pick the covariance kernel and mean function for non-transferable beta ###############
nontran_kernel = function(s,t){return( min(s,t) )}
# nontran_kernel = function(s,t){return( exp(-4*abs(s-t)) )}
nontrans_mu = function(x){ return( 1*cos(2*pi*x) )}
# nontrans_mu = function(x){return( rep(0,length(x)) )}

######### pick diff beta, A size and h ###############
beta_ture_type = 5; ## this type represents three different beta functions.

A_candi = seq(0,20,by=2)
h_candi = seq(4,12,by = 4)
II = 5; JJ = 1 # II and JJ decide the proportion of transferable source and their h

##### generate training data including target and sources ##### 
beta_0 = beta_true(grid$pt)
total_size = 20
target_size = 150
source_size = 100
A_size = A_candi[II]
trans_index =  sample(1:total_size, size = A_size, replace = F)
DATA = gene_transdata(N_target = target_size, N_source = source_size, total_size = total_size, trans_index = trans_index, h = h_candi[JJ],
                       nontrans_mu = nontrans_mu, nontran_kernel = nontran_kernel, grid = grid, sigma = sigma )

##### generate test data ##### 
N_test = 1000
X_test = gene_target_data(N_target = N_test, grid = grid, sigma = sigma)

##### parameter tuning for lambda 
##### all the transfer learning algorithms are using the same factor 
##### here, for all beta0, we set the same scaling factor,
##### but one can also set the lambda_factor to other value based on the complexity of beta0. 
lambda_factor = 1/4


ER = rep(0,5)
#### osflr on target data only 
nontrans_res = sflr(xmat = DATA$X_target, ymat = DATA$Y_target, grid = grid, lambda = lambda_factor*length(DATA$Y_target)^{-rate}, alpha = 1.4, ifadap = T)
ER[1] = ER_func(X = X_test, beta0 = beta_0, beta_hat = nontrans_res$beta_hat, grid = grid, alpha_hat = nontrans_res$alpha_hat)


#### Orcale Trans-Osflr
if(A_size == 0){
  oracle_res = nontrans_res
}else{
  DATA1 = list( X_target = DATA$X_target,
                 Y_target = DATA$Y_target,
                 X_auxiliary = DATA$X_auxiliary[trans_index],
                 Y_auxiliary = DATA$Y_auxiliary[trans_index] )
  oracle_res =  oracle_trans_sflr(Data = DATA1, grid = grid,
                                   lambda1 = lambda_factor*(source_size*length(trans_index))^{-rate},
                                   lambda2 = lambda_factor*(target_size)^{-rate}, 
                                   alpha = 1.4, ifadap = T)
}
ER[2] = ER_func(X = X_test, beta0 = beta_0, beta_hat = oracle_res$beta_hat, grid = grid, alpha_hat = oracle_res$alpha_hat)


#### Detect-TL
detect_res = trans_detect(Data = DATA, grid = grid, alpha = 1.4, ifadap = T, lambda_factor = lambda_factor)
detect_index2 = detect_res$index2
if(is.null(detect_index2)){
  trans_detect_res =  sflr(xmat = DATA$X_target, ymat = DATA$Y_target,
                            grid = grid, lambda = lambda_factor*(target_size)^{-rate}, alpha = 1.4, ifadap = T)
}else{
  DATA2 = list( X_target = DATA$X_target,
                 Y_target = DATA$Y_target,
                 X_auxiliary = DATA$X_auxiliary[detect_index2],
                 Y_auxiliary = DATA$Y_auxiliary[detect_index2] )
  trans_detect_res =  oracle_trans_sflr(Data = DATA2, grid = grid,
                                         lambda1 = lambda_factor*(source_size*length(detect_index2) )^{-rate},
                                         lambda2 = lambda_factor*(target_size)^{-rate},
                                         alpha = 1.4, ifadap = T)
}
ER[3] = ER_func(X = X_test, beta0 = beta_0, beta_hat = trans_detect_res$beta_hat, grid = grid, alpha_hat = trans_detect_res$alpha_hat)


#### EWATL ####
aggtrans_res = trans_sflr(Data = DATA, grid = grid, lambda_factor = lambda_factor)
ER[4] = ER_func(X = X_test, beta0 = beta_0, beta_hat = aggtrans_res$Qagg_res$beta_hat, grid = grid, alpha_hat = aggtrans_res$Qagg_res$alpha_hat)


#### SATL ####
ER[5] = ER_func(X = X_test, beta0 = beta_0, beta_hat = aggtrans_res$Staragg_res$beta_hat, grid = grid, aggtrans_res$Staragg_res$alpha_hat)


names(ER) = c("non-trans", "trans-oracle", "trans-detect", "trans-Q", "trans-star")
ER


######### plot the beta function ##########
plot(grid$pt, nontrans_res$beta_hat, type = "l", col = "red", 
     ylim = range(nontrans_res$beta_hat, beta_0,
                  aggtrans_res$Qagg_res$beta_hat,
                  aggtrans_res$Staragg_res$beta_hat))
lines(grid$pt,beta_0 ,col = "black")
lines(grid$pt,trans_detect_res$beta_hat, col  = "yellow")
lines(grid$pt,aggtrans_res$Qagg_res$beta_hat, col = "blue")
lines(grid$pt,aggtrans_res$Staragg_res$beta_hat, col = "green")

legend(x="topright",legend=c("True","Non-Trans", "Detect-TL", "EWATL", "SATL"),
       col=c("black", "red", "yellow", "blue", "green"), lty=c(1,2,3,4), lwd=c(rep(2,4)), cex=.5, bty="n",
       seg.len=1, text.font = 1, text.width = .2)



