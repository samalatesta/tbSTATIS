#' Fit TB-STATIS
#'
#' @description Estimate maximum likelihood sequence for TB-STATIS. User must input the observed data,
#' set of clinical states and measures to include in the model, and a vector of probabilities indicating the accuracy
#' for each clinical measure.'fit_STATIS' returns the maximum likelihood sequence, and all log-likelihoods calculated during
#' estimation.
#'
#' @param data A data frame of observed data.
#' @param clinical_info A data frame of clinical measures and states.
#' @param p_vec A vector of clinical measure accuracies.
#' @param nstart Number of sequence initializations for maximum likelihood search.
#' @param initial_iter Number of iterations for maximum likelihood search.
#' @return A list of results.

fit_STATIS <- function(data,p_vec,clinical_info, nstart, initial_iter){

  colnames(clinical_info) <- c("clinical_measure", "event_number", "var_name")

  info <- clinical_info %>% dplyr::group_by(clinical_measure) %>% dplyr::summarise(events=dplyr::n())
  info <- data.frame(info)
  colnames(info) <- c("bio", "events")


  #save likelihoods for prelim sequence search
  prelim_like <- vector(mode='list', length=nstart)

  #save sequences for prelim sequence search
  prelim_seq <- vector(mode='list', length=nstart)
  all_likes <- data.frame(start=NA, iter=NA,  like=NA)
  #df of all possible seqs to sample from during likelihood ascent
  add <- possible_seqs(sum(info$events))
  N=dim(info)[1]
  max_like <- rep(NA, nstart)

  #for each nstart, generate an initial sequence then search for mlseq
  for(i in 1:nstart){

    all_seqs <- vector(mode='list', length=nstart)
    #get initial seq to start search
    start_seq <- get_seq(clinical_info)

    start_group <- get_group(start_seq, add)
    #calc first likelihood
    current_likelihood <- get_likelihood(data, start_group,p_vec)[[4]]

    prelim_like_sub <- matrix(nrow=initial_iter, ncol=2)
    prelim_like_sub[1,1] <- current_likelihood
    prelim_like_sub[,2] <- i


    current_seq <- start_group

    #for each iteration generate new seq, calculate likelihood and compare to current likelihood, update if new likelihood is greater
    for(j in 1:initial_iter){

      #print(j)

      bio_info_long <- dplyr::arrange(current_seq, pos)
      #select state to move
      selected_pos = sample(x = 1:dim(bio_info_long)[1], size  = 1, replace=F)
      #elig states to swap with
      elig_pos <- sample(x = setdiff(c(1:dim(bio_info_long)[1]), c(selected_pos)), size  = dim(bio_info_long)[1]-1, replace=F)

      #loop through eligible states, select first one that works
      l=0
      for(k in elig_pos){
        l=l+1
        #k=6
        #print(k)
        bio_info_temp = bio_info_long
        possible_pos <- k
        temp <- bio_info_temp[k,5]
        bio_info_temp[k,5]<- bio_info_temp[selected_pos,5]
        bio_info_temp[selected_pos,5]<- temp
        #print(bio_info_temp)

        #check if eligible sequence based on events within each bio
        check <- bio_info_temp %>% dplyr::arrange(order) %>% dplyr::group_by(bio) %>% dplyr::summarize(Result = all(diff(event) == 1)) %>% dplyr::ungroup()
        #print(check)
        if(all(check$Result)){
          bio_info_long2 <- bio_info_temp %>% dplyr::arrange(order)
          break
        }
        if(l==length(elig_pos)){
          bio_info_long2 <- bio_info_long
        }

      }


      #get new group after swapping states
      new <- get_group(bio_info_long2, add) #%>% dplyr::arrange(pos)

      #get likelihood
      all_seqs[[j]] <- new
      temp_likelihood = get_likelihood(data, new,p_vec)[[4]]


      prelim_like_sub[j,1] <- current_likelihood
      all_likes=rbind(all_likes, data.frame(start=i, iter=j, like=current_likelihood))

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

  #save ml seq
  ml_seq <- prelim_seq[[which.max(max_like)]]

  #save all likes to plot ascent later
  all_likes <- na.omit(all_likes)

  return(list(prelim_like=prelim_like, prelim_seq=prelim_seq, ml=max_like, ml_seq=ml_seq, loglikes=all_likes))

}
