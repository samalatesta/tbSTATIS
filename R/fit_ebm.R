#' Find maximum likelihood estimate
#'
#' @param data A data frame.
#' @param info A data frame
#' @param p_vec Unique row identifier.
#' @param nstart Number of initialization points
#' @param initial_iter Number of iterations
#' @return A list


#Input
#Data set, measures info, number of start points,
fit_ebm <- function(data,p_vec, info, nstart, initial_iter){
  #save likelihoods for prelim sequence search
  prelim_like <- vector(mode='list', length=nstart)

  #save sequences for prelim sequence search
  prelim_seq <- vector(mode='list', length=nstart)

 add <- possible_seqs(dim(info)[1])
N=dim(info)[1]
  max_like <- rep(NA, nstart)
  for(i in 1:nstart){

    all_seqs <- vector(mode='list', length=nstart)

    start_seq <- get_seq(info)

    start_group <- get_group(start_seq, add)

    current_likelihood <- get_likelihood(data, start_group,p_vec)[[4]]

    prelim_like_sub <- matrix(nrow=initial_iter, ncol=2)
    prelim_like_sub[1,1] <- current_likelihood
    prelim_like_sub[,2] <- i




    current_seq <- start_group


    for(j in 1:initial_iter){

      print(j)

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
        check <- bio_info_temp %>% arrange(order) %>% group_by(clinical_measure) %>% summarize(Result = all(diff(event_number) == 1)) %>% ungroup()
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

  return(list(prelim_like, prelim_seq, max_like, ml_seq))

}

