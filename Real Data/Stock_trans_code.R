require(extraDistr)
require(gss)
require(pracma)
require(RandomFieldsUtils)
require(MASS)
require(truncnorm)
require(caret)
require(caTools)

source("Qagg_functions_new.R")
source("source_detect_new.R")

args <-  as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(2022*args)

lambda_factor <- NA

outsample_MSPE <- function(X,Y, beta, alpha){
  ymat_pre <- NULL
  for(i in 1:dim(X)[1]){
    ymat_pre <- c(ymat_pre, alpha + sum(grid$wt*X[i,]*beta))
  }
  PE_nontrans <- sqrt(sum((Y - ymat_pre)^2)/dim(X)[1])
  return(PE_nontrans)
}


Datasets <- c( "May_to_Jun","Jun_to_Jul", "Jul_to_Aug", "Aug_to_Sept")
Kernels <- c("Exp", "M3/2", "M5/2", "M7/2", "Gau")


eachmonth_res <- list()
for(month in Datasets){
  eachkernel_res <- list()
  for(kernel in Kernels){
    #### load stock datasets
    load(file = paste(month,".RData",sep=""))
    sectors_names <- c("BasicIndustries", "CapitalGoods", "ConsumerDurable", "ConsumerNonDurable",
                       "ConsumerServices", "Energy", "Finance", "HealthCare", "PublicUtility",
                       "Technology", "Transportation")

    grid <- gauss.quad(size = 22, interval = c(0,1))


    ## each sector results for each month and kernel
    eachsector_res <- list()
    
    for(SS in sectors_names){
      print(SS)
      ## target dataset: Energy
      xmat_target <- sectors_data_list[[SS]]$X
      ymat_target <- sectors_data_list[[SS]]$Y

      ## source datasets
      source <- sectors_names[ which(sectors_names != SS) ]
      Xmat_auxiliary <- Ymat_auxiliary <- list()
      for(sector in source){
        Xmat_auxiliary[[ which(sector == source) ]] <- sectors_data_list[[sector]]$X
        Ymat_auxiliary[[ which(sector == source) ]] <- sectors_data_list[[sector]]$Y
      }



      prederror <- rep(NA,6)
      names(prederror) <- c("nontrans", "pooled", "oracle-trans", "Detect Trans" ,  "Q-trans", "star-trans")


      id <- sample.split(Y = ymat_target, SplitRatio = .8)
      xmat_target_train <- subset(xmat_target, id == T); ymat_target_train <- subset(ymat_target, id == T)
      xmat_target_test <- subset(xmat_target, id == F); ymat_target_test <- subset(ymat_target, id == F)
      grid <- gauss.quad(size = dim(xmat_target)[2], interval = c(0,0.01*dim(xmat_target)[2]))

      DATA <- list(X_target = xmat_target_train, Y_target = ymat_target_train,
                   X_auxiliary = Xmat_auxiliary, Y_auxiliary = Ymat_auxiliary)


      nontrans_res <- sflr(xmat = xmat_target_train, ymat = ymat_target_train, grid = grid, lambda = NA, alpha = 1.4, ifadap = F)
      prederror[1] <- outsample_MSPE(X = xmat_target_test, Y = ymat_target_test,
                                          beta = nontrans_res$beta_hat, alpha = nontrans_res$alpha_hat)

      ## pooled
      trans_res <- oracle_trans_sflr(Data = DATA, grid = grid, lambda1 = NA, lambda2 = NA, alpha = 1.4, ifadap = F)
      prederror[2] <- outsample_MSPE(X = xmat_target_test, Y = ymat_target_test,
                                          beta = trans_res$beta_pooled, alpha = trans_res$alpha_pooled)

      ## oracle
      prederror[3] <- outsample_MSPE(X = xmat_target_test, Y = ymat_target_test,
                                          beta = trans_res$beta_hat, alpha = trans_res$alpha_hat)


      ############ detect trans ##########
      detect_res <- trans_detect(Data = DATA, grid = grid, alpha = 1.4, lambda_factor = NA, ifadap = F)
      detect_index <- detect_res$index2
      if(is.null(detect_index)){
        det_res <-  sflr(xmat = DATA$X_target, ymat = DATA$Y_target,
                         grid = grid, lambda = NA, alpha = 1.4, ifadap = F)
      }else{
        DATA2 <- list( X_target = DATA$X_target,
                       Y_target = DATA$Y_target,
                       X_auxiliary = DATA$X_auxiliary[detect_index],
                       Y_auxiliary = DATA$Y_auxiliary[detect_index] )
        det_res <-  oracle_trans_sflr(Data = DATA2, grid = grid,
                                      lambda1 = NA,
                                      lambda2 = NA,
                                      alpha = 1.4, ifadap = F)
      }
      prederror[4] <- outsample_MSPE(X = xmat_target_test, Y = ymat_target_test,
                                          beta = det_res$beta_hat, alpha = det_res$alpha_hat)



      ########### Q-agg ################
      aggtrans_res <- trans_sflr(Data = DATA, grid = grid, lambda_factor = NA, ifadap = F)
      prederror[5] <- outsample_MSPE(X = xmat_target_test, Y = ymat_target_test,
                                          beta = aggtrans_res$Qagg_res$beta_hat, alpha = aggtrans_res$Qagg_res$alpha_hat)

      ######### star-agg ##############
      prederror[6] <- outsample_MSPE(X = xmat_target_test, Y = ymat_target_test,
                                          beta = aggtrans_res$Staragg_res$beta_hat, alpha = aggtrans_res$Staragg_res$alpha_hat)


      eachsector_res[[SS]] <- prederror

    }
    
    eachkernel_res[[kernel]] <- eachsector_res
    
  }
  eachmonth_res[[month]] <- eachkernel_res
}






save(eachmonth_res, file = paste("trans_sflr_res",args, ".RData",sep=""))