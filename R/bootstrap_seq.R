#' Find maximum likelihood estimate
#'
#' @param data A data frame.
#' @param clinical_info A data frame
#' @param p_vec Unique row identifier.
#' @param nstart Number of initialization points
#' @param initial_iter Number of iterations
#' @param z number of bootstrap samples
#' @return A list



bootstrap_seq <- function(z, data,clinical_info, p_vec, nstart, initial_iter ){
  boot_ml = data.frame()
  for(i in 1:z){
    boot = data[sample(nrow(data),nrow(data), T),]
    model1 <-  fit_ebm(boot, p_vec, clinical_info, nstart, initial_iter)
    boot_ml = rbind(boot_ml, model1[[4]][,5])

  }
  return(boot_ml)
}
