## code to prepare `TBebmData` dataset
library(dplyr)
#set seed
set.seed(123)


  #clinical measures info
  # 6 events, 3 clinical measures
  seq <- data.frame(bio=c("bio1", "bio1", "bio2", "bio2", "bio2", "bio3"), event=c(1,2,1,2,3,1),position=c(1,2,3,4,5,6), order=c(6, 1,3,2,4,5),  sub=c(1,1,1,2,3,4))
  p=.95
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

  usethis::use_data(TBebmData, overwrite = TRUE)



  library(ggplot2)
  library(reshape2)

  long <- dat3 %>% dplyr::select(-C1, -C2, -stage) %>% melt(id.var="Index")%>% dplyr::arrange(value)
  dat4 <- dat3 %>% dplyr::arrange(stage)
  #cough, weightloss,smear pos, cavity, sweat,infiltrates, smear 3, hemop

  #heatmaps to visualize data
  long$Index <- factor(long$Index, levels = c(dat4$Index))
  long$variable <- factor(long$variable, levels = c( "X6","X1", "X3","X2", "X4",   "X5"))


ggplot(long, aes(Index,variable)) +
    geom_tile(aes(fill = factor(value)), colour = "white") + scale_fill_manual(name="Levels", values = c("white", "#6099C6", "#C660C6")) +xlab("Individual Participants") + ylab("Clinical Measure") + theme_bw()+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", text=element_text(size=14))
  ggsave("heatmap_ex.png", width=6, height=4)


  #start with first event/biomarker
  get_group <- function( info=data.frame(), all_seqs){
    #D=make_D(N)
    #seq=get_seq(D)
    #info<- seq

    info <- info %>% dplyr::arrange(order)
    N=dim(info)[1]
    #all_seqs=add

    all_seqs2=all_seqs
    all_seqs=data.frame(t(apply(all_seqs[,1:N], 1, function(i) paste(i, info$bio) )))
    all_seqs$unique= apply(all_seqs[,1:N], 1, function(x) length(unique(x)))==N

    final <-  cbind(all_seqs2, data.frame(unique=all_seqs$unique)) %>% dplyr::filter(unique==T)
    select <- sample(1:dim(final)[1],1)

    info$sub<- as.numeric(final[select,1:N])

    #print(info)
    return(info)
  }

  #### Initialize sequence events
  get_seq <- function(D=data.frame()){

    #transform D to long form
    bio_info_long<- data.frame(bio=NA, event=NA)
    for(i in 1: dim(D)[1]){
      add <- data.frame(bio=rep(D[i,1]), event = c(1:D[i,2]))
      bio_info_long <- rbind(bio_info_long, add)

    }

    bio_info_long <- bio_info_long %>% dplyr::filter(is.na(event)==F)

    #bio_info_long <- ml[[2]][[1]]
    N=dim(bio_info_long)[1]
    bio_info_long$pos = c(1:N)


    all_perms <-permn(bio_info_long$bio)
    r <- sample(1:length(all_perms),1)
    temp <- all_perms[[r]]


    bio_info_long_temp <- bio_info_long
    bio_info_long_temp$temp <- temp
    bio_info_long_temp <- bio_info_long_temp %>% dplyr::arrange(event,temp)
    bio_info_long_temp$order <- c(1:dim(bio_info_long_temp)[1])

    bio_info_long_new <- bio_info_long_temp %>% dplyr::arrange(pos) %>% dplyr::select(-temp)
    #bio_info_long_new
    return(bio_info_long_new)

  }

  #start with first event/biomarker


  #S=new_order
  calculate_likelihood <- function(data=data.frame(), S=data.frame(), p=vector()){
    #S=t
    #p<- c(.8,.8,.8, .8, .8)
    #data=dat
    #p=c(.9,.9,.9,.9)
    bio <- data.frame(data[,1:ncol(data)])

    # number of biomarkers
    N = as.numeric(dim(bio)[2])

    #individuals in data set
    M = dim(bio)[1]

    new_bio <- bio

    #multiply each p by each corresponding variable
    for(i in 1:dim(new_bio)[2]){
      x<- S$bio[i]
      n <- sum(S$bio==x)
      new_bio[,i] <- ifelse(new_bio[,i]==1, dbinom(1,1,p[i]), dbinom(0,1,p[i])/max(n,1))
    }


    #empty matrix to store individual stage probabilities
    p_perm_k = matrix(NA,M,N+1)

    #order based on sequence S and length
    order = as.numeric(S$pos)

    group = as.numeric(S$sub)

    #reorder biomarker columns based on S
    bio_order = new_bio

    #add alpha to make more continuous
    #alpha=.1

    #bio_order=bio_order+alpha
    normal <- 1-bio_order

    #special case stage 1
    normal_prob = normal %>% mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
    #normal_prob[which(normal_prob==0)] <- .1

    tot_prob_stage = normal_prob[,1]

    p_perm_k[,1] <- tot_prob_stage
    S2=S

    for (i in unique(group)) {
      #print(i)

      if(i < max(group)){
        #abnormal biomarkers assuming S is true sequence
        abnormal_cols <- as.numeric(row.names(S2[S2$sub<=i,]))

        abnormal <- data.frame(bio_order[,abnormal_cols])

        #normal data assuming S is true sequence
        normal_cols = as.numeric(rownames(S2[!(rownames(S2 )%in% abnormal_cols),]))

        normal <- 1-data.frame(bio_order[,normal_cols])
        normal[normal==-1]<- 1
        abnormal_n <- apply(abnormal, 1, sum)
        normal_n <- apply(normal, 1, function(x) sum(x==0))



        abnormal_prob = abnormal %>% mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
        #abnormal_prob[which(abnormal_prob==0)] <- .1

        normal_prob = normal %>% mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
        #normal_prob[which(normal_prob==0)] <- .1
        #fill in empty matrix
        #abnormal_prob <- ifelse(sum(abnormal_n,normal_n)>0 & abnormal_n==0, .5*normal_prob, abnormal_prob)
        #normal_prob <- ifelse(sum(abnormal_n,normal_n)>0 & normal_n==0, .5*abnormal_prob, normal_prob)


        tot_prob_stage = abnormal_prob[,1]*normal_prob[,1]
        #tot_prob_stage= (abnormal+normal)/dim(bio_order)[1]


        p_perm_k[,i+1] <-   tot_prob_stage

      }

      if(i == max(group)){
        #special case stage 1
        abnormal <- bio_order

        abnormal_prob = abnormal %>% mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
        #abnormal_prob[which(abnormal_prob==0)] <- .1


        p_perm_k[,i+1] <- abnormal_prob[,1]
      }




    }
    #p_perm_k <-  ifelse(p_perm_k==0,y,p_perm_k)

    prob_subj = p_perm_k*(1/(max(group)+1))
    #prob_subj <- prob_subj+0.01
    total_prob_subj = apply(prob_subj,1,sum, na.rm=T)

    loglike = sum(log(total_prob_subj + 1e-250))

    return(list(p_perm_k, prob_subj, total_prob_subj, loglike))
  }



  #Input
  #Data set, measures info, number of start points,
  find_ml <- function(data,p_vec, D, nstart, initial_iter){

    #save likelihoods for prelim sequence search
    prelim_like <- vector(mode='list', length=nstart)

    #save sequences for prelim sequence search
    prelim_seq <- vector(mode='list', length=nstart)

    N=ncol(data)
    all  <- expand.grid(rep(list(N), N))

    add  <- all

    add$max <- do.call(pmax, add)
    #add$keep<-N
    add$min <-  do.call(pmin, add[,1:N])
    add$sum <- rowSums(add[,1:N])
    add$Var2max <- add$Var2<=2
    add$Var3max <- add$Var3<=3
    add$Var4max <- add$Var4<=4
    add$Var5max <- add$Var5<=5
    add$Var6max <- add$Var6<=6

    add <- add %>% dplyr::filter(min ==1 & sum <= sum(1:N) & Var1==1 & Var2max==T & Var3max==T & Var4max==T & Var5max==T& Var6max==T )

    add$keep =apply(add[,1:N], 1, function(i) all(diff(i) >= 0 & diff( i)  <=1))
    add <- add %>% dplyr::filter(min ==1 & sum <= sum(1:N) & keep==T)


    all_likes <- data.frame(start=NA, iter=NA, like=NA)

    max_like <- rep(NA, nstart)
    for(i in 1:nstart){

      all_seqs <- vector(mode='list', length=nstart)

      start_seq <- get_seq(D)

      start_group <- get_group(start_seq, add)

      current_likelihood <- calculate_likelihood(data, start_group,p_vec)[[4]]

      prelim_like_sub <- matrix(nrow=initial_iter, ncol=2)
      prelim_like_sub[1,1] <- current_likelihood
      prelim_like_sub[,2] <- i




      current_seq <- start_group


      for(j in 1:initial_iter){

        print(j)

        bio_info_long <- current_seq  %>% dplyr::arrange(pos)

        selected_pos = sample(x = 1:dim(bio_info_long)[1], size  = 1, replace=F)

        elig_pos <- sample(x = setdiff(c(1:N), c(selected_pos)), size  = N-1, replace=F)


        l=0
        for(k in elig_pos){
          l=l+1
          #print(k)
          bio_info_temp = bio_info_long
          possible_pos <- k
          temp <- bio_info_temp[k,4]
          bio_info_temp[k,4]<- bio_info_temp[selected_pos,4]
          bio_info_temp[selected_pos,4]<- temp
          #print(bio_info_temp)

          #check if eligible sequence based on events within each bio
          check <- bio_info_temp %>% arrange(order) %>% group_by(bio) %>% summarize(Result = all(diff(event) == 1)) %>% ungroup()
          #print(check)
          if(all(check$Result)){
            bio_info_long2 <- bio_info_temp %>% dplyr::arrange(order)
            break
          }
          if(l==length(elig_pos)){
            bio_info_long2 <- bio_info_long
          }

        }

        #start_group
        #bio_info_long2
        #get new group
        new <- get_group(bio_info_long2, add) %>% dplyr::arrange(pos)


        all_seqs[[j]] <- new
        temp_likelihood = calculate_likelihood(data, new,p_vec)[[4]]


        prelim_like_sub[j,1] <- current_likelihood
        all_likes=rbind(all_likes, data.frame(start=i, iter=j, like=current_likelihood))

        #print(new)
        #if current sequence improves likelihood, update current likelihood and sequence
        if(temp_likelihood > current_likelihood){

          prelim_like_sub[j,1] <- temp_likelihood

          current_likelihood <-  temp_likelihood

          current_seq <- new

        }


      }
      #max_pos <- which.max(prelim_like_sub[,1])
      prelim_like[[i]] <- current_likelihood
      prelim_seq[[i]] <- current_seq
      max_like[i] <- current_likelihood



    }

    ml_seq <- prelim_seq[[which.max(max_like)]]

    #find maximum likelihood and corresponding sequence

    return(list(prelim_like, prelim_seq, max_like, ml_seq, all_likes))

  }

  #clinical measures info
  # 6 events, 3 clinical measures
  D <- data.frame(bio=c("bio1", "bio1", "bio2", "bio2", "bio2", "bio3"), events=c(1,2,1,2,3,1),position=c(1,2,3,4,5,6))
  p=.9

  prep <- dat3 %>% dplyr::select(starts_with("X"))
  trust <- find_ml(prep, c(.9, .9, .9, .9, .9, .9), D, 3, 500)
data=prep
nstart=3
niter=500
N=6
