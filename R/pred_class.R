#' Get disease class
#' @param data A data frame.
#' @param S Data frame
#' @param p A vector
#' @return A data frame.

pred_class <- function(data=data.frame(), S=data.frame(), p=vector()){

#data=dplyr::arrange(data,Index)
likelihood=get_likelihood(data, S, p )

class_probs <- likelihood[[2]]


class_probs=data.frame(class_probs)
colnames(class_probs) = c(0:max(S$sub))

pred_class=as.numeric(colnames(class_probs)[apply(class_probs,1,which.max)])

classdf <- cbind(data,pred_class)

return(classdf)
}


#S=new_order
#measure_info <- data.frame(clinical_measure=c("measure1", "measure1", "measure2", "measure2", "measure2", "measure3"),event_number = c(1,2,1,2,3,1), event_name=c(c("X1", "X2", "X3", "X4", "X5", "X6")))

#data <- TBebmData[,c(1,3:8)]
#S=TBebm::get_seq(measure_info)
#poss=TBebm::possible_seqs(6)
#S <- TBebm::get_group(S, poss)
#p<- c(.8,.8,.8, .8, .8,.8)

#head(stagedf)
