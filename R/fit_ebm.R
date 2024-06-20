#' Find maximum likelihood estimate
#'
#' @param data A data frame.
#' @param info A data frame
#' @param p_vec Unique row identifier.
#' @param nstart Number of initialization points
#' @param initial_iter Number of iterations
#' @return A list
#initial_iter=500
#Input
#info <- D

#fit_ebm(data, p_vec, info, nstart=1, initial_iter=100)
#Data set, measures info, number of start points,
fit_ebm <- function(data,p_vec,clinical_info, nstart, initial_iter){
  #clinical_info=measure_info
  info <- clinical_info %>% dplyr::group_by(clinical_measure) %>% dplyr::summarise(events=dplyr::n())
  info <- data.frame(info)
  colnames(info) <- c("bio", "events")
  #save likelihoods for prelim sequence search
  prelim_like <- vector(mode='list', length=nstart)

  #save sequences for prelim sequence search
  prelim_seq <- vector(mode='list', length=nstart)
  all_likes <- data.frame(start=NA, iter=NA, like=NA)
  add <- possible_seqs(sum(info$events))
  N=dim(info)[1]
  max_like <- rep(NA, nstart)
  for(i in 1:nstart){

    all_seqs <- vector(mode='list', length=nstart)

    start_seq <- get_seq(info)
start_seq
    start_group <- get_group(start_seq, add)

    current_likelihood <- get_likelihood(data, start_group,p_vec)[[4]]

    prelim_like_sub <- matrix(nrow=initial_iter, ncol=2)
    prelim_like_sub[1,1] <- current_likelihood
    prelim_like_sub[,2] <- i




    current_seq <- start_group


    for(j in 1:initial_iter){

      #print(j)

      bio_info_long <- dplyr::arrange(current_seq, pos)

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

      #start_group
      #bio_info_long2
      #get new group
      new <- get_group(bio_info_long2, add) %>% dplyr::arrange(pos)


      all_seqs[[j]] <- new
      temp_likelihood = get_likelihood(data, new,p_vec)[[4]]


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

  ml_seq <- cbind(ml_seq, event_name=data.frame(event_name=clinical_info$event_name))
  #ml_seq$est_seq <- ml_seq$sub
  #ml_seq <- ml_seq %>% dplyr::select(-pos,-order, -sub)

  all_likes <- na.omit(all_likes)
  #find maximum likelihood and corresponding sequence

  return(list(prelim_like=prelim_like, prelim_seq=prelim_seq, ml=max_like, ml_seq=ml_seq, loglikes=all_likes))

}


