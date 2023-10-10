require(extraDistr)
require(gss)
require(pracma)
require(RandomFieldsUtils)
require(MASS)
require(truncnorm)
require(caret)
require(caTools)
# read.csv("C:/Users/haoti/Dropbox/research/transfer learning/code/trans_sflr_code/Q-aggregation/real_data/Stock-Dataset/ConsumerDurable.cvs")


## data preparasion
sectors_names <- c("BasicIndustries", "CapitalGoods", "ConsumerDurable", "ConsumerNonDurable",
             "ConsumerServices", "Energy", "Finance", "HealthCare", "PublicUtility",
             "Technology", "Transportation")
sectors_list <- list()

for(sector in sectors_names){
  path <- paste("C:/Users/haoti/Dropbox/research/transfer learning/code/trans_sflr_code/Q-aggregation/real_data/Stock-Dataset/",sector,".cvs", sep = "" )
  # path <- paste("~/Dropbox/research/transfer learning/code/trans_sflr_code/Q-aggregation/real_data/Stock-Dataset/",sector,".cvs", sep = "" )
  sectors_list[[sector]] <- read.csv(path)
}


## we use April as predictor and May as reponse
sectors_data_list <- list()
Date1 <- "2021-05-03"
Date2 <- "2021-05-28"
Date3 <- "2021-06-30"
for(sector in sectors_names){
  print(sector)
  X <- sectors_list[[sector]]
  companies_names <- unique(X$Name)
  grid <- gauss.quad(size = dim( subset(X, Name == companies_names[1] & as.Date(Date) >= Date1 & as.Date(Date) <= Date2) )[1], interval = c(0,1))
  print( length(companies_names) )
  X_curves <- matrix(NA, nrow = length(companies_names), ncol = length(grid$pt))
  Y <- rep(NA,length(companies_names))
  for(name in companies_names){
    Xcomp <- subset(X, Name == name & as.Date(Date) >= Date1 & as.Date(Date) <= Date2)
    if( length((Xcomp$Close - Xcomp$Close[1])/Xcomp$Close[1]) != length(grid$pt) ){
      next
    }else{
      X_curves[which(name == companies_names),] <- (Xcomp$Close - Xcomp$Close[1])/Xcomp$Close[1]
      
      Close_price1 <- subset(X, Name == name & as.Date(Date) == Date2)$Close
      Close_price2 <- subset(X, Name == name & as.Date(Date) == Date3)$Close
      
      Y[which(name == companies_names)] <- (Close_price2 - Close_price1) /Close_price1
    }
  }
  rownames(X_curves) <- companies_names
  names(Y) <- companies_names
  
  X_curves <- na.omit(X_curves)
  Y <- Y[!is.na(Y)]
  
  sectors_data_list[[sector]] <- list( X = 100*X_curves,
                                       Y = 100*Y)
  
  # plot(x = grid$pt, y = X_curves[1,], type = "l", ylim = range(X_curves))
  # for(i in 1:dim(X_curves)[1]){
  #   lines(x = grid$pt, y = X_curves[i,])
  # }
}


## check if there are NA in dataset
for(sector in sectors_names){
  X <- sectors_data_list[[sector]]$X
  Y <- sectors_data_list[[sector]]$Y
  
  if(sum(is.na(X)) > 0 | sum(is.na(Y)) > 0){
    print(sector)
  }
}


save(sectors_data_list, file = "May_to_Jun.RData")






