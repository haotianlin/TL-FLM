################## transferable source detection ####################


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
      target_index <- sample( x = c(1:N_T), size = .5*N_T, replace = F)
      
      train_X <- X_target[target_index, ]
      train_Y <- Y_target[target_index]
      test_X <- X_target[-target_index, ]
      test_Y <- Y_target[-target_index]
      
      lambda1 <- lambda_factor*(length(Y_trans) + length(train_Y))^{-6/7}; lambda2 <- lambda_factor*(length(train_Y))^{-6/7}
      
      
      # non-trans
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


