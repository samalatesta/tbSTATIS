
#' Positional variance diagram
#'
#' @description
#' Plot maximum likelihood sequence with sequences fit using bootstrap resampling to visualize
#' uncertainty about maximum likelihood sequence.
#'
#' @param boot_seqs A data frame of estimated sequences obtained from bootstrap resampling and
#' estimating TB-STATIS on each resample.
#' @param ml Maximum likelihood sequence from fitting TB-STATIS to original data.
#' @return A ggplot object.


pvd <- function(boot_seqs=data.frame(), ml=data.frame()){
#ml=model1$ml_seq
z <- dim(boot_seqs)[1]
ml<- model1$ml_seq

ml <- dplyr::arrange(ml,bio, event)

ns <- data.frame(matrix(ncol=dim(boot_seqs)[2],nrow=dim(boot_seqs)[2]))
#colnames(ns) <- ml$var_name
store=data.frame(pos=factor(1:dim(ml)[1]))

for(i in 1:dim(boot_seqs)[2]){
  event_ns=data.frame(table(boot_seqs[,i]))
  colnames(event_ns) <-c("pos", "n")
  store2 <- dplyr::left_join(store, event_ns, by="pos")
  store2[is.na(store2)] <- 0
  colnames(ns)=ml$var_event
  ns[,i] <- store2$n
}

ns=data.frame(as.matrix(ns)*(1/z))
ns$rowname=rownames(ns)

pvd <- ns %>% tidyr::gather(colname, value,-rowname)
ml <- dplyr::arrange(ml,sub)
ml$var_name <- factor(ml$var_name)
pvd$colname <- factor(pvd$colname, levels=ml$var_name)
#pvd$value <- pvd$value/z
pvd_plot=ggplot2::ggplot(pvd, ggplot2::aes(x = rowname, y = colname, fill = value)) +
  ggplot2::geom_tile() + ggplot2::xlab("Disease Severity Class") + ggplot2::ylab("Clinical State") + ggplot2::scale_fill_gradient(name="",low = "white", high = "#3A68AB", breaks=c(0,.25,.5,.75,1)) + ggplot2::theme_bw() + ggplot2::theme(text=ggplot2::element_text(size=14, color="black"))

d <- data.frame(Y=ml$var_name, X=ml$sub)
dat <- merge(pvd, d)

pvd_plot2=pvd_plot + ggplot2::geom_tile(data=dat, ggplot2::aes(X,Y), fill="transparent", colour="black", linewidth=1.5)

return(pvd_plot2)

}


