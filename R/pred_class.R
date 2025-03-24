#' Get predicted disease class
#' @description
#' Get predicted disease class for each individual in sample given their observed data and an estimated
#' disease sequence from fitting TB-STATIS.
#'
#' @param data A data frame of observed data.
#' @param S A data frame for the estimated disease sequence.
#' @param p A vector of clinical state accuracies.
#' @return A data frame.

pred_class <- function(data=data.frame(), S=data.frame(), p=vector()){

    #get likelhood
    likelihood=get_likelihood(data, S, p )

    #save class probs
    stage_probs <- likelihood[[2]]
    stage_probs=data.frame(stage_probs)
    #update colnames to class numbers
    colnames(stage_probs) = c(0:max(S$sub))

    #assign disease class with max prob
    pred_stage=as.numeric(colnames(stage_probs)[apply(stage_probs,1,which.max)])

    #cbind to original data
    classdf <- cbind(data,pred_stage)

    return(classdf)
  }

