#' Propose simultaneous clinical states
#' @description
#' Generate a simulataneous disease sequence given a sequence of clinical states.
#'
#' @param info A data frame.
#' @param all_seqs Data frame of all possible disease sequences.
#' @return A data frame.
#' @keywords internal
#'

get_group <- function( info=data.frame(), all_seqs){
  #arrange sequence of states
  info <- info %>% dplyr::arrange(order)
  N=dim(info)[1]

  #select disease sequence allowing multiple states per disease class
  #code allows for multiple states in one measure to occur at the same class
  #uncomment code below if clinical measures should be unique per class
  all_seqs2=all_seqs
  all_seqs=data.frame(t(apply(all_seqs[,1:N], 1, function(i) paste(i, info$bio) )))
  all_seqs$unique= apply(all_seqs[,1:N], 1, function(x) length(unique(x)))==N

  #all possible seqs after criteria above
  final <-  cbind(all_seqs2, data.frame(unique=all_seqs$unique))#%>% dplyr::filter(unique==T)
  #randomly select one sequence to use
  select <- sample(1:dim(final)[1],1)

  #add to info df
  info$sub<- as.numeric(final[select,1:N])

  # print(info)
  return(info)
}


#' Sequence initialization
#' @description
#' Generate initial sequence as start point for maximum likelihood search.
#'
#' @param info A data frame containing the clinical measures and corresponding clinical
#' states to generate a sequence from.
#' @return A data frame.
#' @keywords internal
#'

get_seq <- function(info=data.frame()){
  #save info df to work with
  bio_info_long <- info
  #update col names
  colnames(bio_info_long) <- c("bio", "event", "var_name")

  N=dim(bio_info_long)[1]
  #save original order of states
  bio_info_long$pos = c(1:N)

  #sample from clinical measures w/o replacement as new order and save as temp df
  temp<-sample(bio_info_long$bio, N)

  #update order in df
  bio_info_long_temp <- bio_info_long
  bio_info_long_temp$temp <- temp
  bio_info_long_temp <- bio_info_long_temp %>% dplyr::arrange(temp)
  bio_info_long_temp$order <- c(1:dim(bio_info_long_temp)[1])
  #bio_info_long_temp
  bio_info_long_new <- bio_info_long #%>% dplyr::arrange(pos) %>% dplyr::select(-temp)
  bio_info_long_new$order <- bio_info_long_temp$pos

  return(bio_info_long_new)
}





#' All possible sequences
#'
#' @description
#' Generate data frame of all possible simultaneous sequences for N clinical states.
#'
#' @param N Number of clinical states.
#' @return A data frame.
#' @keywords internal
#'
#

possible_seqs <- function(N=numeric()){
  #filter such that first disease class is 1 and are monotonically non-decreasing
  filter_rows <- function(g, x) {
    ok  <- function(z) all(diff(z) %in% 0:1)
    out <- g[apply(g, 1, ok), ]
    replace(out, TRUE, lapply(out, \(i) x[i]))
  }
  #generate all seqs
  f2 <- function(x = c(1:N), n=N+1, n1=2) {
    data.frame(as.list(rep(1, n1)),
               gtools::combinations(length(x), n-n1, repeats.allowed = TRUE)) |>
      filter_rows(x)
  }
  #generate all possible sequences of size N and filter using filter_rows and f2 functions
  all=f2()
  all <- all[,-1]
  return(all)
}

