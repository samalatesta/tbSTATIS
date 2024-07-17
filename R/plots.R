#' Plot clinical states
#'
#' @param data A data frame.
#' @param id.var Unique row identifier.
#' @param state.vars Variable names for clinical states.
#' @return A ggplot object.


plot_states <- function(data=data.frame(), id.var=character(), state.vars=vector()){

  states <- data[,c(id.var, state.vars)]

  states <- dplyr::arrange(states,states[,state.vars])

  props <- data.frame(states=c(state.vars), prop=rep(NA,length(state.vars)))

  for(i in state.vars){
    props[props$states==i,2] <- prop.table(table(states[,i]))[2]

  }

  props <- dplyr::arrange(props, -prop)
  props$states <- factor(props$states)

  long <- reshape2::melt(states, id.var=id.var)

  long <- dplyr::arrange(long,value)

  long[,id.var] <- factor( long[,id.var], levels=states[,id.var])

  long$variable <- factor(long$variable, levels = props$states)

  long$value <- factor(long$value)

  event_plot <- ggplot2::ggplot(long, ggplot2::aes(long[,id.var],variable)) +
                ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white") +
                ggplot2::scale_fill_manual(name="Levels", values = c("white", "#6099C6")) +
                ggplot2::xlab("Individual Participants") +
                ggplot2::ylab("Clinical State") +
                ggplot2::theme_bw()+
                ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", text=ggplot2::element_text(size=14))

  return(event_plot)

}


#' Plot disease class
#'
#' @param pred_class A vector.
#' @return A ggplot object.


plot_class <- function(pred_class=vector()){

  classdf <- data.frame(class=c(0:max(pred_class)), prop = as.vector(prop.table(table(pred_class))))
  class_plot <- ggplot2::ggplot(data=classdf) +
                ggplot2::geom_bar(ggplot2::aes(x=class, y=prop), stat="identity", fill="#3A68AB") +
                ggplot2::theme_bw() + ggplot2::xlab("Disease Severity Class")+
                ggplot2::ylab("Proportion") + ggplot2::scale_x_continuous(breaks=c(0:max(pred_class)))+
                ggplot2::theme(text=ggplot2::element_text(size=14))
  return(class_plot)

}



#' Plot likelihood ascent
#' @description
#' The function `plot_likes` function generates a plot of the likelihood ascent from fitting the EBM. The x-axis is the iteration number, the y-axis is the log-likelihood, and each line represents the log-likelihood at each iteration for a single initialized sequence.
#'
#' @param likes A data frame with three columns (start, iter, like).
#' @return A ggplot object.


plot_likes <- function(likes=data.frame()){

  likes_plot <- ggplot2::ggplot(data=likes) +
                ggplot2::geom_line(ggplot2::aes(x=iter, y=like, color=factor(start))) +
                ggplot2::theme_bw() + ggplot2::ylab("Log-likelihood") +
                ggplot2::xlab("Iteration") +
                ggplot2::scale_color_discrete(name="Initialized Sequence") +
                ggplot2::xlim(0,max(likes$iter)) +
                ggplot2::theme(legend.position = "top", text = ggplot2::element_text(size=14))

  return(likes_plot)

}

