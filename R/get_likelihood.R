
#get_likelihood(data,S,p)
#' Calculate log-likelihood
#'
#' @param data A data frame.
#' @param S Data frame
#' @param p_vec A vector
#' @return A list.

#S=current_seq
get_likelihood <- function(data=data.frame(), S=data.frame(), p_vec=vector()){

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
    new_bio[,i] <- ifelse(new_bio[,i]==1, dbinom(1,1,p_vec[i]), dbinom(0,1,p_vec[i])/max(n,1))
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
  normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
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



      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
      #abnormal_prob[which(abnormal_prob==0)] <- .1

      normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
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

      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
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
