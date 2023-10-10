###################################################################
# This file contains the data generation process in simulation.
###################################################################


################## generate true beta #####################
beta_true <- function(x)
{
  
  
  if(beta_ture_type == 1){
    val <- 0
    f_k <- rep(0,K)
    for(k in 1:K){
      f_k[k] <- 4*sqrt(2)*(-1)^{k-1}*k^{-2}
      val <- val + f_k[k]*phi_K(x,k)
    }
    return(val)
  }else if(beta_ture_type == 2){
    return( 4*cos(3*pi*x) )
  }else if(beta_ture_type == 3){
    return( 4*sin(3*pi*x) + 4*cos(3*pi*x) )
  }else if(beta_ture_type == 4){
    return( 10*(2*x - 0.5)^2  )
  }else{
    return( 0.3*(1*dnorm(x, 0.15, 0.04) + 0.5*dnorm(x, 0.5, 0.04) - 1*dnorm(x, 0.8, 0.04)) )
  }
  
}

################## generate transferable beta with h #####################
beta_trans <- function(x,h){
  
  # R <- rsign(K)
  # R <- runif(K, min = -1, 1)
  # R <- rtruncnorm(K, a=-1, b=1, mean = 0, sd = .5)
  val <- 0
  
  f_k <- rep(0,K)
  for(k in 1:K){
    # f_k[k] <- R[k]*((sqrt(6)*h)/(pi * k^2))
    R <- runif(1, min = -1, 1)
    # R <- rsign(1)
    f_k[k] <- R*((sqrt(6)*h)/(pi * k^2))
    val <- val + f_k[k]*phi_K(x,k)
  }
  return(val + beta_true(x))
  
}


beta_trans1 <- function(x,h, kk = c(3,6,9)){
  M <- length(kk)
  
  val <- 0
  for(k in kk){
    # R <- runif(1, min = -1, 1)
    R <- rsign(1)
    f_k <- R*h/(k*sqrt(M))
    val <- val + f_k*phi_K(x,k)
  }
  return(val + beta_true(x))
}




################ generate non-transferable beta ########################
################ using gaussian process ##########################
beta_nontrans <- function(x, mu , K){
  Mu <- mu(x)
  # Sigma <- outer(x, x, K)
  
  Sigma <- sapply(x, function(s1) {
    sapply(x, function(s2) {
      K(s1, s2)
    })
  })
  
  
  beta <- 5*mvrnorm(mu = Mu, Sigma = Sigma)
  return(beta)
}





################### generate predictor ###############
# X_func <- function(x, r2)
# {
#   val <- 0
#   zeta <- rnorm(n = K, mean = 0, sd = (seq(1,K))^{-r2})
#   for(k in 1:K){
#     val <- val + zeta[k]*phi_K(x,k)
#   }
#   val <- val + sin(pi*x)
#   out <- list(X = val,
#               zeta = zeta)
#   return(out)
# }

X_func <- function(x, C )
{
  val <- sin(pi*x)
  Sigma = outer(x,x,C)
  
  val <- val + mvrnorm(mu = rep(0,length(x)),Sigma = Sigma)
  out <- list(X = val)
  return(out)
}




##################### generate X ###########################
gene_target_data <- function(N_target = 50, grid, sigma = 0.5){
  Xmat <- Ymat <- NULL
  for(i in 1:N_target){
    X <- X_func(x = grid$pt, C)
    Xmat <- rbind(Xmat, X$X)
    Y <- sum(grid$wt*beta_0*X$X) + rnorm(n = 1, mean = 0, sd = sigma) 
    Ymat <- c(Ymat,Y)
  }
  out <- list(X = Xmat, Y = Ymat)
  return(out)
}




################ generate all_dataset ######################
gene_transdata <- function(N_target, N_source, total_size, trans_index, h, 
                           nontrans_mu = function(x){ return(5*cos(2*pi*x)) },
                           nontran_kernel = function(s,t){ return(exp(-abs(s-t))) },
                           grid, sigma = 0.5){
  
  ##### for target dataset ########
  Xmat_target <- Ymat_target <- NULL
  for(i in 1:N_target){
    X <- X_func(x = grid$pt, C)
    Xmat_target <- rbind(Xmat_target, X$X)
    Y <- sum(grid$wt*beta_0*X$X) + rnorm(n = 1, mean = 0, sd = sigma)  
    Ymat_target <- c(Ymat_target,Y)
  }
  
  
  ##### for source dataset ########
  Xmat_auxiliary <- Ymat_auxiliary <- list()
  beta_trans <- c()
  beta_nontrans <- c()
  for(I in 1:total_size){
    if(I %in% trans_index){
      beta1 <- beta_trans(x = grid$pt, h = h)
      # beta1 <- beta_trans1(x = grid$pt, h = h)
      Xmat <- Ymat <- NULL
      for(i in 1:N_source){
        X <- X_func(x = grid$pt, C)
        Xmat <- rbind(Xmat, X$X)
        Y <- sum(grid$wt*beta1*X$X) + rnorm(n = 1, mean = 0, sd = sigma) 
        Ymat <- c(Ymat,Y)
      }
      Xmat_auxiliary[[I]] <- Xmat
      Ymat_auxiliary[[I]] <- Ymat
      beta_trans <- rbind(beta_trans,beta1)
      
    }else{
      beta2 <- beta_nontrans(x = grid$pt, mu = nontrans_mu, K = nontran_kernel)
      Xmat <- Ymat <- NULL
      for(i in 1:N_source){
        X <- X_func(x = grid$pt, C)
        Xmat <- rbind(Xmat, X$X)
        Y <- sum(grid$wt*beta2*X$X) + rnorm(n = 1, mean = 0, sd = sigma)
        Ymat <- c(Ymat,Y)
      }
      Xmat_auxiliary[[I]] <- Xmat
      Ymat_auxiliary[[I]] <- Ymat
      beta_nontrans <- rbind(beta_nontrans,beta2)
    }
    
  }
  names(Xmat_auxiliary) <- paste("A", 1:total_size, sep = "")
  names(Ymat_auxiliary) <- paste("A", 1:total_size, sep = "")
  
  
  out <- list(X_target = Xmat_target,
              Y_target = Ymat_target,
              X_auxiliary = Xmat_auxiliary,
              Y_auxiliary = Ymat_auxiliary,
              trans_index = trans_index, 
              beta_trans = beta_trans,
              beta_nontrans = beta_nontrans)
  return(out)
  
}







