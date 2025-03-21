## code to prepare `TBData` dataset
library(dplyr)
#set seed
set.seed(123)

source("/projectnb/cbs/samalate/tbSTATIS/data-raw/functions_sim_data.R")

  N=6
  p=.9
  M=250
  #generate number of clinical measures and states per measure
  D=make_D(N)
  info=D
  bio_info_long<- data.frame(bio=NA, event=NA)
  for(t in 1: dim(info)[1]){
    add <- data.frame(bio=rep(info[t,1]), event = c(1:info[t,2]))
    bio_info_long <- rbind(bio_info_long, add)
  }
  bio_info_long <- na.omit(bio_info_long)
  bio_info_long$var_name <- paste0("V", 1:dim(bio_info_long)[1])

  #make data set
  dat <- make_data(bio_info_long, p, M)
  data<- dat[[1]]

  #true disease sequence
  data.frame(dat[[2]])

  #add covariates associated with stage
  C1 <- rep(NA, M)
  C2 <- rep(NA, M)
  for(i in 1:M){
    C1[i] <- rbinom(1,1, min((data$stage[i]+ runif(1,0,.1))/max(data$stage),1))
    C2[i] <- rnorm(1, 0+(data$stage[i]*.5), 1)
  }
  C2 <- round(C2, digits=3)
  dat2 <- cbind(data, C1, C2)

  colnames(dat2)[2] <- "class"
  colnames(dat2)[3:8] <- paste0("X", 1:6)
  dat2$class <- dat2$class-1
  #set up data frame
  TBData <- dat2

  usethis::use_data(TBData, overwrite = TRUE)

