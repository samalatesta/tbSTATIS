
#' Calculate log-likelihood
#' @description
#' Calculate likelihood given observed data, sequence S, vector of clinical measure accuracies p_vec.
#'
#' @param data A data frame.
#' @param S A data frame.
#' @param p_vec A vector
#' @return A list.
#'
get_likelihood <- function(data=data.frame(), S=data.frame(), p_vec=vector()){
  #input data
  bio <- data.frame(data[,1:ncol(data)])

  # number of biomarkers
  N = as.numeric(dim(bio)[2])

  #individuals in data set
  M = dim(bio)[1]

  #new df to become probabilities from 0/1 data
  new_bio <- bio

  #multiply each p by each corresponding variable
  for(i in 1:dim(new_bio)[2]){
    #x<- S$bio[i]
    #n <- sum(S$bio==x)
    new_bio[,i] <- ifelse(new_bio[,i]==1, dbinom(1,1,p_vec[i]), dbinom(0,1,p_vec[i]))
  }


  #empty matrix to store individual stage probabilities
  p_perm_k = matrix(NA,M,N+1)

  #order based on sequence S and length
  group = as.numeric(S$sub)

  #reorder biomarker columns based on S
  bio_order = new_bio

  #bio_order=bio_order+alpha
  normal <- 1-bio_order

  #special case class 1
  normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)

  #save probs for class 1
  tot_prob_stage = normal_prob[,1]

  p_perm_k[,1] <- tot_prob_stage
  S2=S

  for (i in 1:length(unique(group))) {
    #print(i)

    if(i < max(group)){
      #abnormal biomarkers assuming S is true sequence
      abnormal_cols <- S2$var_name[S2$sub<=i]

      abnormal <- data.frame(bio_order[,abnormal_cols])

      #normal data assuming S is true sequence
      normal_cols = S2 %>% filter(!(var_name %in% abnormal_cols)) %>% select(var_name)

      normal <- 1-data.frame(bio_order[,normal_cols[,1]])
      normal[normal==-1]<- 1

      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)

      normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)

      tot_prob_stage = abnormal_prob[,1]*normal_prob[,1]
      #tot_prob_stage= (abnormal+normal)/dim(bio_order)[1]

      #save probs
      p_perm_k[,i+1] <-   tot_prob_stage

    }
    #special case max class, all states should have occurred
    if(i == max(group)){
      #special case stage 1
      abnormal <- bio_order

      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
      #abnormal_prob[which(abnormal_prob==0)] <- .1


      p_perm_k[,i+1] <- abnormal_prob[,1]
    }

  }
  #divide probs by total classess
  prob_subj = p_perm_k*(1/(max(group)+1))

  #multiply across rows to get individual likelihood for class k
  total_prob_subj = apply(prob_subj,1,sum, na.rm=T)

  #full likelhood by summing log likes for each ind.
  loglike = sum(log(total_prob_subj + 1e-250))


  return(list(p_perm_k, prob_subj, total_prob_subj, loglike))
}
