#' Propose simultaneous events
#'
#' @param info A data frame.
#' @param all_seqs List of all possible sequences.
#' @return A ggplot object.
#' @keywords internal
#'


get_group <- function( info=data.frame(), all_seqs=data.frame()){
  #D=make_D(N)
  #seq=get_seq(D)
  #info<- seq

  info <- info %>% dplyr::arrange(order)
  N=dim(info)[1]
  #all_seqs=add

  all_seqs2=all_seqs
  all_seqs=data.frame(t(apply(all_seqs[,1:N], 1, function(i) paste(i, info$clinical_measure) )))
  all_seqs$unique= apply(all_seqs[,1:N], 1, function(x) length(unique(x)))==N

  final <-  cbind(all_seqs2, data.frame(unique=all_seqs$unique)) %>% dplyr::filter(unique==T)
  select <- sample(1:dim(final)[1],1)

  info$sub<- as.numeric(final[select,1:N])

  #print(info)
  return(info)
}

#measure_info <- data.frame(clinical_measure=c("measure1", "measure1", "measure2", "measure2", "measure2", "measure3"),event_number = c(1,2,1,2,3,1), event_name=c(c("X1", "X2", "X3", "X4", "X5", "X6")))


#' Propose sequence as start point for maximum likelihood estimation.
#'
#' @param info A data frame.
#' @return A data frame.
#' @keywords internal
#'
#### Initialize sequence events
get_seq <- function(info=data.frame()){

  bio_info_long <- info

  N=dim(bio_info_long)[1]
  bio_info_long$pos = c(1:N)


  all_perms <-combinat::permn(bio_info_long[,1])
  r <- sample(1:length(all_perms),1)
  temp <- all_perms[[r]]


  bio_info_long_temp <- bio_info_long
  bio_info_long_temp$temp <- temp
  bio_info_long_temp <-  dplyr::arrange(bio_info_long_temp,event_number,temp)
  bio_info_long_temp$order <- c(1:dim(bio_info_long_temp)[1])

  bio_info_long_new <- dplyr::arrange(bio_info_long_temp,pos)
  bio_info_long_new <- dplyr::select(bio_info_long_temp,-temp)

  #bio_info_long_new
  return(bio_info_long_new)

}




#' Generate data frame of all possible simultaneous sequences.
#'
#' @param N Number of events.
#' @return A data frame.
#' @keywords internal
#'
#


possible_seqs <- function(N=numeric()){

  filter_rows <- function(g, x) {
    ok  <- function(z) all(diff(z) %in% 0:1)
    out <- g[apply(g, 1, ok), ]
    replace(out, TRUE, lapply(out, \(i) x[i]))
  }
  f2 <- function(x = c(1:N), n=N+1, n1=2) {
    data.frame(as.list(rep(1, n1)),
               gtools::combinations(length(x), n-n1, repeats.allowed = TRUE)) |>
      filter_rows(x)
  }

  # test runs
  all=f2() # as per question
  all <- all[,-1]
  return(all)
}

