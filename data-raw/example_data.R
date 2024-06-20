## code to prepare `TBebmData` dataset
library(dplyr)
#set seed
set.seed(123)


  #clinical measures info
  # 6 events, 3 clinical measures
  seq <- data.frame(bio=c("bio1", "bio1", "bio2", "bio2", "bio2", "bio3"), event=c(1,2,1,2,3,1),position=c(1,2,3,4,5,6), order=c(6, 1,3,2,4,5),  sub=c(1,1,1,2,3,4))
  p=.9
  M=250
  #### make data
  N= sum(dim(seq)[1])
  all_seqs  <- expand.grid(rep(list(1:N), N))

  add  <- all_seqs
  add$max <- do.call(pmax, add)

  add$min <-  do.call(pmin, add[,1:N])
  add$sum <- rowSums(add[,1:N])
  add$Var2max <- add$Var2<=2
  add$Var3max <- add$Var3<=3
  add$Var4max <- add$Var4<=4
  add$Var5max <- add$Var5<=5
  add$Var6max <- add$Var6<=6
  add <- add %>% dplyr::filter(min ==1 & sum <= sum(1:N) & Var1==1 & Var2max==T & Var3max==T & Var4max==T & Var5max==T & Var6max==T)

  add$keep =apply(add[,1:N], 1, function(i) all(diff(i) >= 0 & diff( i)  <=1))
  add <- add %>% dplyr::filter(min ==1 & sum <= sum(1:N) & keep==T)

  #save true order and arrange data frame to be bio, event
  S <- seq %>% dplyr::arrange(order)
  true_order <- S$order

  #generate stage for individuals 1:M, randomly sample from all possible stages with replacement
  dat <- data.frame(Index=c(1:M),stage=sample(1:(max(S$sub)+1), M, replace=T))

  #empty matrix of M rows, K-1 columns to store event variables
  bio <- as.data.frame(matrix(nrow=M,ncol=max(S$pos)))

  #combine true stage and empty matrix to fill in
  dat2 <- cbind(dat, bio)

  check <- matrix(nrow=M, ncol=N+1)
  check[,1] <- dat2$stage
  #loop across 1:M
  for(i in 1:M){

    #special case when stage=1, when p=1, all events should = 0
    if(dat2$stage[i]==1){

      for(j in 1:ncol(bio)){
        x<- S$bio[j]
        n <- sum(S$bio==x)
        dat2[i,j+2] <- rbinom(1,1,(1-p))
        check[i,j+1] <- (1-p)
      }

    }


    if(dat2$stage[i]!=1){

      S_temp <- S %>% dplyr::filter(sub <= (dat2$stage[i])-1)
      to_fill <- dim(S_temp)[1]
      for(j in 1:(to_fill)){
        dat2[i,j+2] <- rbinom(1,1,p)
        check[i,j+1] <- p
      }

    }


    for(j in 1:ncol(bio)){
      x<- S$bio[j]
      n <- sum(S$bio==x)
      if(is.na(dat2[i,j+2])){
        dat2[i,j+2] <- rbinom(1,1,1-p)
        check[i,j+1] <- (1-p)
      }

    }

  }
  dat2 <- dat2[,c("Index", "stage", "V2","V4","V3", "V5","V6", "V1")]
  #change col names
  colnames(dat2)[3:8] <- c("X1", "X2", "X3", "X4", "X5", "X6")
  dat2$stage <- dat2$stage-1
  #add covariates associated with stage
  C1 <- rep(NA, M)
  C2 <- rep(NA, M)
  for(i in 1:M){
  C1[i] <- rbinom(1,1, min((dat2$stage[i]+ runif(1,0,.1))/max(dat2$stage),1))
  C2[i] <- rnorm(1, 0+(dat2$stage[i]*.5), 1)
  }
  C2 <- round(C2, digits=3)
  dat3 <- cbind(dat2, C1, C2)

  #set up data frame
  TBebmData <- dat3


  info=measure_info
  nstart=1
  initial_iter=2
  data=TBebmData[,c(1,3:8)]
  p_vec=c(.95,.95,.95,.95,.95,.95)
  ml<- info
  ml$sub <- c(1,2,1,3,4,1)
  ml2<- info
  ml2$sub <- c(1,4,1,2,1,3)
  x=fit_ebm(data, p_vec, info, nstart=4, initial_iter=300 )
  TBebm::get_likelihood(data, ml, p_vec)

  usethis::use_data(TBebmData, overwrite = TRUE)

