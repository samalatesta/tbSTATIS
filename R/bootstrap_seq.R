#' Quantify sequence uncertainty with bootstrap resampling
#'
#' @param data A data frame of observed data.
#' @param clinical_info A data frame with clinical measures and clinical states.
#' @param p_vec A vector of clinical measure accuracies to fit TB-STATIS.
#' @param nstart Number of initialized sequences for each run of TB-STATIS.
#' @param initial_iter Number of iterations for TB-STATIS.
#' @param z Number of bootstrap resamples to perform.
#' @return A data frame of sequences estimated for bootstrap resamples.


bootstrap_seq <- function(z, data,clinical_info, p_vec, nstart, initial_iter ){
  boot_ml = data.frame()
  for(i in 1:z){
    boot = data[sample(nrow(data),nrow(data), T),]
    model1 <-  fit_STATIS(boot, p_vec, clinical_info, nstart, initial_iter)
    boot_ml = rbind(boot_ml, model1[[4]][,6])

  }
  return(boot_ml)
}
