###################################################################
# This file defines all the functions in simulation.
# To obtain the estimator for the functional linear regression, we follow the representer theorem.
###################################################################


phi_K <- function(x,k){
  return( sqrt(2)*cos(k*pi*x) )
}


RKHS_kernel <- function(x1,x2, r1 = 1)
{
  val <- 0
  for(k in 1:100){
    val <- val + (1/(k)^2*r1)*phi_K(x = x1, k = k) *phi_K(x = x2, k = k) # different from version 1
  }
  return(val)
  
}

Matern_kernel <- function(x1,x2, kernel = "Exp", rho = 1){
  Ds <- abs(x1-x2)
  if(kernel=="Exp"){
    Sig_new <- exp(-Ds/rho)
  }
  if(kernel=="M3/2"){ # kernel of Matern Process nu=3/2
    Sig_new <- (1+sqrt(3)*(-Ds/rho))*exp(-sqrt(3)*Ds/rho)
  }
  if(kernel=="M5/2"){ # kernel of Matern Process nu=5/2
    Sig_new <- (1+(sqrt(5)*Ds/rho)+(5*(Ds^2)/(3*rho^2)))*exp(-sqrt(5)*Ds/rho)
  }
  if(kernel=="M7/2"){ # kernel of Matern Process nu=7/2
    Sig_new <- (1+(sqrt(7)*Ds/rho)+(14*(Ds^2)/(5*rho^2))+(7*sqrt(7)*Ds^3)/(15*rho^3))*exp(-sqrt(7)*Ds/rho)
  }
  if(kernel=="Gau"){
    Sig_new <- exp(-(Ds^2)/(2*rho))
  }
  return( Sig_new )
}


############## self sflr #################
sflr <- function(xmat, ymat, grid, lambda = NA, alpha = 1.4, ifadap = F){
  N.size <- dim(xmat)[1]
  xmat_center <- sweep(xmat, 2 , colMeans(xmat))
  ymat_center <- ymat - mean(ymat)
  
  Kmat <- outer(grid$pt, grid$pt, RKHS_kernel, r1 = r1)
  # Kmat <- outer(grid$pt, grid$pt, RKHS_kernel)
  r <- Kxmat <- Kmat%*%(grid$wt*t(xmat_center))
  q <- t(Kxmat)%*%(grid$wt*t(xmat_center))
  
  
  if(is.na(lambda) == T){
    if( ifadap ){
      lambda <- adap_lambda(xmat = xmat, ymat = ymat, grid = grid)
    }else{
      GCV_res <- GCV_sflr1(xmat = xmat, ymat = ymat, grid = grid, alpha = 1.4)
      lambda <- GCV_res$lambda
    }
  }
  alpha <- solve(q + N.size*lambda*diag(1,N.size))%*%ymat_center
  beta_hat <- t(alpha)%*%t(r)
  
  intercept_hat <- mean(ymat) - sum(grid$wt*colMeans(xmat)*beta_hat)
  
  if(is.na(lambda) == T){
    out <- list(beta_hat = beta_hat,
                alpha_hat = intercept_hat,
                lambda = lambda)
  }else{
    out <- list(beta_hat = beta_hat,
                alpha_hat = intercept_hat,
                lambda = lambda)
  }
  return(out)
}




############### oracle trans sflr #############
oracle_trans_sflr <- function(Data, grid, lambda1 = NA, lambda2 = NA, alpha = 1.4,  ifadap = F){
  ########################## transfer step ##############################
  X_trans <- Data$X_auxiliary[[1]]
  Y_trans <- Data$Y_auxiliary[[1]]
  if(length(Data$Y_auxiliary)>=2){
    for(i in 2:length(Data$X_auxiliary)){
      X_trans <- rbind(X_trans, Data$X_auxiliary[[i]])
      Y_trans <- c(Y_trans, Data$Y_auxiliary[[i]])
    }
  }
  
  
  if(is.na(lambda1) == T){
    if( ifadap ){
      lambda1 <- adap_lambda(xmat = X_trans, ymat = Y_trans, grid = grid)
    }else{
      GCV_res1 <- GCV_sflr1(xmat = X_trans, ymat = Y_trans, grid = grid, alpha = 1.4)
      lambda1 <- GCV_res1$lambda
    }
  }
  
  
  
  trans.res <- sflr(xmat = X_trans, ymat = Y_trans, grid = grid, lambda = lambda1, alpha = alpha, ifadap = ifadap)
  beta_pooled <- trans.res$beta_hat
  alpha_pooled <- trans.res$alpha_hat
  
  ########################## debiased step ##############################
  Y_new <- NULL
  for(i in 1:length(Data$Y_target)){
    Y_new <- c(Y_new, Data$Y_target[i] - alpha_pooled - sum(grid$wt*beta_pooled*Data$X_target[i,]))
  }
  
  if(is.na(lambda2) == T){
    if( ifadap ){
      lambda2 <- adap_lambda(xmat = Data$X_target, ymat = Y_new, grid = grid)
    }else{
      GCV_res2 <- GCV_sflr1(xmat = Data$X_target, ymat = Y_new, grid = grid, alpha = 1.4)
      lambda2 <- GCV_res2$lambda
    }
  }
  debiased.res <- sflr(xmat = Data$X_target, ymat = Y_new, grid = grid, lambda = lambda2, alpha = alpha, ifadap = ifadap)
  
  beta_delta <- debiased.res$beta_hat
  alpha_delta <- debiased.res$alpha_hat
  
  
  
  beta_hat <- beta_pooled + beta_delta
  alpha_hat <- alpha_pooled + alpha_delta
  
  
  out <- list(beta_hat  = beta_hat,
              alpha_hat = alpha_hat,
              beta_pooled = beta_pooled,
              alpha_pooled = alpha_pooled,
              beta_delta = beta_delta,
              alpha_delta = alpha_delta,
              lambda1 = lambda1,
              lambda2 = lambda2)
  
  return(out)
}






############### trans sflr #############
trans_sflr <- function(Data, grid, lambda_factor = 1, ifadap = F){
  X_target <- Data$X_target
  Y_target <- Data$Y_target
  K <- length(Data$Y_auxiliary)
  
  
  
  id1 <- sample.split(Y = Y_target, SplitRatio = .7)
  X_train <- subset(X_target, id1 == T); Y_train <- subset(Y_target, id1 == T)
  X_train_center <- sweep(X_train, 2 , colMeans(X_train));  Y_train_center <- Y_train - mean(Y_train)
  
  X_test <- subset(X_target, id1 == F); Y_test <- subset(Y_target, id1 == F)
  X_test_center <- sweep(X_test, 2 , colMeans(X_test));  Y_test_center <- Y_test - mean(Y_test)
  
  
  N_T <- length(Y_train)
  init_res <- sflr(xmat = X_train, ymat = Y_train, grid = grid, lambda = lambda_factor*N_T^{-6/7}, alpha = 1.4, ifadap = ifadap)
  # calculate the MSE scores
  Rank_norm <- rep(0,K)
  for(k in 1:K){
    X_trans <- Data$X_auxiliary[[k]]; Y_trans <- Data$Y_auxiliary[[k]]
    N_A <- length(Y_trans)
    init_k <- sflr(xmat = X_trans, ymat = Y_trans, grid = grid, lambda = lambda_factor*N_A^{-6/7}, alpha = 1.4, ifadap = ifadap)
    
    Rank_norm[k] <- RKHS_norm(beta = t(init_res$beta_hat - init_k$beta_hat), grid = grid)
  }
  Tset <- list()
  kk_list <- unique(rank(Rank_norm))
  for(k in 1:length(kk_list)){
    Tset[[k]] <- which(rank(Rank_norm) <= kk_list[k])
  }
  Tset <- unique(Tset)
  
  
  alpha_T <- rep(0, 1+length(Tset) )
  beta_T <- list()
  beta_T[[1]] <- as.vector(init_res$beta_hat)
  alpha_T[1] <- init_res$alpha_hat
  for(k in 1:length(Tset)){
    # trans
    DATA1 <- list( X_target = X_train,
                   Y_target = Y_train,
                   X_auxiliary = Data$X_auxiliary[Tset[[k]]],
                   Y_auxiliary = Data$Y_auxiliary[Tset[[k]]] )
    trans_res <- oracle_trans_sflr(Data = DATA1, grid = grid,
                                   lambda1 = lambda_factor*((length(Tset[[k]])*N_A)^{-6/7}),
                                   lambda2 = lambda_factor*(N_T^{-6/7}), alpha = 1.4, ifadap = ifadap)
    beta_T[[k+1]] <- as.vector(trans_res$beta_hat)
    alpha_T[k+1] <- trans_res$alpha_hat
  }
  beta_T<-beta_T[!duplicated((beta_T))]
  beta_T<- as.matrix(as.data.frame(beta_T))
  colnames(beta_T) <- seq(1:ncol(beta_T))
  
  
  
  agg_res <- Qagg_func(B = beta_T, A = alpha_T, X_test = X_test, Y_test = Y_test, total_step = 100, selection = F, grid = grid)
  agg_res$alpha_hat <- mean(Y_train) - sum(grid$wt*colMeans(X_train)*agg_res$beta_hat)
  
  fs <- matrix(data = NA, nrow = nrow(beta_T), ncol = 5)
  fss <- rep(0, 5)
  # half_index <- which(as.numeric(lapply(X = Tset, FUN = length)) <=5)
  for(k in 1:5){
    # Staragg_res <- Staragg_func(B = beta_T[,half_index], A = alpha_T[half_index], X_test = X_test, Y_test = Y_test, grid = grid, conflevel = 0.001, Data = Data)
    Staragg_res <- Staragg_func(B = beta_T, A = alpha_T, X_test = X_test, Y_test = Y_test, grid = grid, conflevel = 0.001, Data = Data)
    fs[,k] <- Staragg_res$beta_hat
    fss[k] <- Staragg_res$alpha_hat
  }
  
  Staragg_res$beta_hat <- rowMeans(fs)
  # Staragg_res$alpha_hat <- mean(fss)
  Staragg_res$alpha_hat <- mean(Y_train) - sum(grid$wt*colMeans(X_train)*Staragg_res$beta_hat)
  
  
  
  return(list( Qagg_res = agg_res,
               Staragg_res = Staragg_res,
               Tset = Tset) )
  
  
}



############# Q_aggregation ####################
Qagg_func <- function(B, A, X_test, Y_test, grid , total_step = 100, selection = F){
  if(sum(B==0)==ncol(B)*nrow(B)){
    return(rep(0,nrow(B)))
  }
  K <- ncol(B)
  theta_hat <- matrix(0,nrow = K, ncol = 1)

  # Three different options for 
  # ## low
  # Temp = 1/5
  # ## mid
  # Temp = 4
  ## high
  Temp = length(Y_test)
  
  if(selection){
    
  }else{
    for(k in 1:K){
      theta_hat[k,] <- exp(-MSE_func(X = X_test, Y = Y_test, beta = B[,k], alpha = A[k], grid = grid)/Temp)
    }
    theta_hat <- theta_hat/(sum(theta_hat))
    theta_old <- theta_hat
    alpha <- A%*%theta_hat
    beta <- B%*%theta_hat
    beta_we <- beta
    for(ss in 1:total_step){
      for(k in 1:K){
        theta_hat[k] <- exp(-MSE_func(X = X_test, Y = Y_test, beta = B[,k], alpha = A[k], grid = grid)/Temp +
                              MSE1_func(X = X_test, Y = Y_test, beta1 = B[,k], beta2 = beta, alpha1 = A[k], alpha2 = alpha, grid = grid)/(4*Temp) )
      }
      theta_hat <- theta_hat/(sum(theta_hat))
      beta <- 1/4*B%*%theta_hat + 3/4*beta
      alpha <- 1/4*A%*%theta_hat + 3/4*alpha
      if( sum(abs(theta_hat - theta_old))<10^{-4}  ){break}
      theta_old <- theta_hat
    }
  }
  # print(ss)
  
  return(list(beta_hat = beta,
              alpha_hat = alpha,
              theta = theta_hat,
              beta_we = beta_we))
}



############## star-shape aggregation ###########
Staragg_func <- function(B, A, X_test, Y_test, grid, conflevel = 0.01, Data){
  
  id2 <- c(1:length(Y_test))
  index <- sample(x = id2, size = .7*length(id2), replace = F)
  X1 <- X_test[index,]; Y1 <- Y_test[index]
  X2 <- X_test[-index,]; Y2 <- Y_test[-index]


  id2 <- sample.split(Y = Y_test, SplitRatio = .7)
  X1 <- subset(X_test, id2 == T); Y1 <- subset(Y_test, id2 == T)
  X2 <- subset(X_test, id2 == F); Y2 <- subset(Y_test, id2 == F)
  

  
  R_func <- function(f, a, X, Y, grid){
    trans_pre <- NULL
    for(j in 1:length(Y)){
      trans_pre <- c(trans_pre, a + sum(grid$wt*X[j,]*f) )
    }
    RMSE_trans <- sum((trans_pre - Y)^2)/length(Y)
    return(RMSE_trans)
  }
  
  fnorm_func <- function(f, a, X, grid){
    trans_pre <- NULL
    for(j in 1:nrow(X)){
      trans_pre <- c(trans_pre, a + sum(grid$wt*X[j,]*f) )
    }
    return( sqrt( sum( trans_pre^2 )/nrow(X) ) )
  }
  
  R1 <- R2 <- rep(0, ncol(B))
  for(k in 1:ncol(B)){
    R1[k] <- R_func(f = B[,k], a = A[k], X = X1, Y = Y1, grid = grid)
    R2[k] <- R_func(f = B[,k], a = A[k], X = X2, Y = Y2, grid = grid)
  }
  
  f1_index <- which.min(R1)
  f1 <- B[,f1_index]
  a1 <- A[f1_index]
  
  f1norm <- f2norm <- rep(0, ncol(B))
  for(k in 1:ncol(B)){
    f1norm[k] <- fnorm_func(f = B[,k] - f1, a = A[k] - a1, X = X1, grid = grid)
    f2norm[k] <- fnorm_func(f = B[,k] - f1, a = A[k] - a1, X = X2, grid = grid)
  }
  
  
  # construct preselection elements
  b_func <- function(X, Y, B, A){
    pred <- matrix(NA, nrow = length(Y), ncol = ncol(B))
    for(k in 1:ncol(B)){
      for(i in 1:length(Y)){
        pred[i,k] <- A[k] + sum(grid$wt*X[i,]*B[,k])
      }
    }
    b <- max( max(Y), max(pred) )
    return( b )
  }
  M <- ncol(B);  n0 <- length(Y_test)
  b <- b_func(X = Data$X_target, Y = Data$Y_target, B = B, A = A)
  Phi <- b*sqrt( (log(M) + conflevel)/n0 )
  c <- 4*(1+9*b)
  
  F1 <- c()
  for(k in 1:M){
    if( R1[k] <= (R1[f1_index] + c*max( Phi*f1norm[k], Phi^2 )) ){
      F1 <- c(F1, k)
    }
  }
  
  convex_mse <- c()
  for(k in F1){
    lambda_k <- min( max(0.5*((R2[k] - R2[f1_index])/(f2norm[k]^2) + 1),0) , 1) 
    convex_mse <- c(convex_mse, R_func(f = (1-lambda_k)*B[,k] + lambda_k*f1, a = (1-lambda_k)*A[k] + lambda_k*a1, X = X2, Y = Y2, grid = grid))
  }
  
  
  j_index <- F1[which.min(convex_mse)]
  lambda_j <- min( max(0.5*((R2[j_index] - R2[f1_index])/(f2norm[j_index]^2) + 1),0) , 1) 
  
  
  
  
  f_titled <- lambda_j*f1 + (1 - lambda_j)*B[,j_index]
  a_titled <- lambda_j*a1 + (1 - lambda_j)*A[j_index]
  
  return(list(beta_hat = f_titled,
              alpha_hat = a_titled,
              index = c(j_index, f1_index),
              lambda = c((1-lambda_j), lambda_j)))
}







################ two MSE functions ###########################
MSE_func <- function(X, Y, beta, alpha, grid){
  trans_pre <- NULL
  for(j in 1:length(Y)){
    trans_pre <- c(trans_pre, alpha + sum(grid$wt*X[j,]*beta) )
  }
  # RMSE_trans <- sum((trans_pre - Y)^2)/length(Y)
  RMSE_trans <- sum((trans_pre - Y)^2)
  return(RMSE_trans)
}

MSE1_func <- function(X, Y, beta1, beta2, alpha1, alpha2, grid){
  trans_pre <- NULL
  for(j in 1:length(Y)){
    trans_pre <- c(trans_pre, sum(grid$wt*X[j,]*(beta1 - beta2)) + alpha1 - alpha2 )
  }
  # RMSE_trans <- sum((trans_pre)^2)/length(Y)
  RMSE_trans <- sum((trans_pre)^2)
  return(RMSE_trans)
}

################# RKHS_norm ######################
RKHS_norm <- function( beta, grid){
  KK <- outer(grid$pt, grid$pt, RKHS_kernel, r1 = r1)
  # KK <- outer(grid$pt, grid$pt, RKHS_kernel)
  Hnorm <- t(beta)%*%solve(KK)%*%beta
  return(Hnorm)
}



###### compute excess risk
ER_func <- function(X, beta0, beta_hat, alpha_hat, grid){
  ER <- 0
  for(j in 1:N_test){
    # ER <- ER + sum(grid$wt*( (beta_hat - beta0)*X$X[j,] )^2  )/N_test
    ER <- ER + (sum(grid$wt*beta0*X$X[j,]) - alpha_hat -  sum(grid$wt*beta_hat*X$X[j,]))^2/N_test
  }
  return(ER)
}



################################### GCV adaptive settings ############################
mysolve = function(R){
  chol2inv(chol( R ))
}

GCV_sflr1 <- function(xmat, ymat, grid, alpha = 1.4){
  N.size <- dim(xmat)[1]
  xmat <- sweep(xmat, 2 , colMeans(xmat))
  ymat <- ymat - mean(ymat)
  
  Kmat <- outer(grid$pt, grid$pt, RKHS_kernel)
  r <- Kxmat <- Kmat%*%(grid$wt*t(xmat))
  q <- t(Kxmat)%*%(grid$wt*t(xmat))
  
  
  # lower_lambda = 0.01*(N.size^{-Rate}); upper_lambada = 10*(N.size^{-Rate})
  lower_lambda = 10^{-6}; upper_lambada = 10^{-1}
  
  lambda = exp(seq(log(lower_lambda),log(upper_lambada),length.out = 20)) 
  
  
  GCV_func <- function(lambda, r, q, ymat, N.size ){
    GCV_score <- 0
    
    s <- cbind(apply(grid$wt*t(xmat), 2, sum), apply( grid$wt*grid$pt*t(xmat), 2, sum))
    W <- q + N.size*lambda*diag(rep(1,N.size))
    d <- mysolve( t(s)%*%mysolve(W)%*%s )%*%t(s)%*%mysolve(W)
    c <- mysolve(W)%*%( diag(rep(1,N.size)) - s%*%mysolve( t(s)%*%mysolve(W)%*%s )%*%t(s)%*%mysolve(W) )
    
    
    A_temp <- s%*%d + q%*%c
    A1_temp <- diag(1,N.size) - A_temp
    # GCV_score <- ( N.size^{-1}* t(A1_temp%*%ymat)%*%(A1_temp%*%ymat) )/(N.size^{-1}*sum(diag( diag(1,N.size)-alpha*A_temp )))^2
    GCV_score <- (  t(A1_temp%*%ymat)%*%(A1_temp%*%ymat)/N.size )/( 1 - sum(diag(A_temp))/N.size  )^2
    # GCV_score <- (  t(A1_temp%*%ymat)%*%(A1_temp%*%ymat)/N.size )/( (1 - 4/N.size)^2 )
    return(GCV_score)
  }
  
  res = sapply(lambda, GCV_func, r = r, q = q, ymat = ymat, N.size = N.size)
  
  
  return(list(lambda = lambda[which.min(res)]))
  
}

#####################################################################
################## transferable source detection ####################
# this function implements the detection algorithm in Ye Tian's paper.
#####################################################################
trans_detect <- function( Data, grid, alpha = 1.4,  ifadap = F, Eps = .01, lambda_factor ){
  
  X_target <- Data$X_target
  Y_target <- Data$Y_target
  
  index <- c()
  index2 <- c()
  total_MSE <- list()
  
  for(i in 1:length(Data$X_auxiliary)){
    X_trans <- Data$X_auxiliary[[i]]
    Y_trans <- Data$Y_auxiliary[[i]]
    
    DATA <- list( X_target = X_target,
                  Y_target = Y_target,
                  X_auxiliary = list(X_trans),
                  Y_auxiliary = list(Y_trans) )
    N_T <- length(Y_target)
    N_A <- length(Y_trans)
    
    MSE_nontrans <- c()
    MSE_trans <- c()
    MSE_trans2 <- c()
    
    for(II in 1:20){
      target_index <- sample( x = c(1:N_T), size = .7*N_T, replace = F)
      
      train_X <- X_target[target_index, ]
      train_Y <- Y_target[target_index]
      test_X <- X_target[-target_index, ]
      test_Y <- Y_target[-target_index]
      
      lambda1 <- lambda_factor*(length(Y_trans) + length(train_Y))^{-6/7}; lambda2 <- lambda_factor*(length(train_Y))^{-6/7}
      
      
      # non-trans
      # nontrans_res <- sflr(xmat = X_target[target_index, ], ymat = Y_target[target_index],
      #                       grid = grid, lambda = NA, alpha = 1.4, ifadap = T)
      nontrans_res <- sflr(xmat = train_X, ymat = train_Y,
                           grid = grid, lambda = lambda2, alpha = 1.4, ifadap = ifadap)
      # non-trans MSE
      nontrans_pre <- NULL
      for(j in 1:length(test_Y)){
        nontrans_pre <- c(nontrans_pre, nontrans_res$alpha_hat + sum(grid$wt*test_X[j,]*nontrans_res$beta_hat))
      }
      RMSE_nontrans <- sqrt(sum((nontrans_pre - test_Y)^2)/length(test_Y))
      MSE_nontrans <- c(MSE_nontrans, RMSE_nontrans)
      
      # trans
      DATA1 <- list( X_target = train_X,
                     Y_target = train_Y,
                     X_auxiliary = list(X_trans),
                     Y_auxiliary = list(Y_trans) )
      # trans_res <- trans_sflr(Data = DATA1,
      #                          grid = grid, lambda1 = NA, lambda2 = NA, alpha = 1.4, ifadap = T)
      trans_res <- oracle_trans_sflr(Data = DATA1,
                                     grid = grid, lambda1 = lambda1, lambda2 = lambda2, alpha = 1.4, ifadap = ifadap)
      
      # debias1
      trans_pre <- NULL
      for(j in 1:length(test_Y)){
        trans_pre <- c(trans_pre, trans_res$alpha_hat + sum(grid$wt*test_X[j,]*trans_res$beta_hat))
      }
      RMSE_trans <- sqrt(sum((trans_pre - test_Y)^2)/length(test_Y))
      MSE_trans <- c(MSE_trans, RMSE_trans)
      
      # debias1
      trans_pre2 <- NULL
      for(j in 1:length(test_Y)){
        trans_pre2 <- c(trans_pre2, trans_res$alpha_pooled + sum(grid$wt*test_X[j,]*trans_res$beta_pooled))
      }
      RMSE_trans2 <- sqrt(sum((trans_pre2 - test_Y)^2)/length(test_Y))
      MSE_trans2 <- c(MSE_trans2, RMSE_trans2)
      
      
    }
    
    source_MSE <- cbind(MSE_nontrans, MSE_trans, MSE_trans2)
    total_MSE[[i]] <- source_MSE
    source_MSE <- Re(source_MSE)
    if( colMeans(source_MSE)[1] > (1 + Eps)*colMeans(source_MSE)[2]){
      index <- c(index,i)
    }
    
    if( colMeans(source_MSE)[1] > (1 + Eps)*colMeans(source_MSE)[3]){
      index2 <- c(index2,i)
    }
  }
  
  source_MSE <- c()
  for(i in 1:length(total_MSE)){
    source_MSE <- rbind(source_MSE,colMeans(total_MSE[[i]]))
  }
  
  # return(index)
  out <- list(index = index,
              index2 = index2,
              source_MSE = source_MSE,
              total_MSE = total_MSE)
  return(out)
  
}











