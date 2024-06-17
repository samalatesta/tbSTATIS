#' Plot clinical events
#'
#' @param data A data frame.
#' @param id.var Unique row identifier.
#' @param event.vars Variable names for clinical events.
#' @return A ggplot object.


plot_events <- function(data=data.frame(), id.var=character(), event.vars=vector()){

  events <- data[,c(id.var, event.vars)]

  events <- dplyr::arrange(events,events[,event.vars])

  props <- data.frame(events=c(event.vars), prop=rep(NA,length(event.vars)))

  for(i in event.vars){
    props[props$events==i,2] <- prop.table(table(events[,i]))[2]

  }

  props <- dplyr::arrange(props, -prop)
  props$events <- factor(props$events)

  long <- reshape2::melt(events, id.var=id.var)

  long <- dplyr::arrange(long,value)

  long[,id.var] <- factor( long[,id.var], levels=events[,id.var])

  long$variable <- factor(long$variable, levels = props$events)

  long$value <- factor(long$value)

  event_plot <- ggplot2::ggplot(long, ggplot2::aes(long[,id.var],variable)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white") + ggplot2::scale_fill_manual(name="Levels", values = c("white", "#6099C6")) + ggplot2::xlab("Individual Participants") + ggplot2::ylab("Clinical Event") + ggplot2::theme_bw()+ ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", text=ggplot2::element_text(size=14))

  return(event_plot)

}


#' Plot disease stage
#'
#' @param pred_stage A vector.
#' @return A ggplot object.


plot_stage <- function(pred_stage=vector()){

  stagedf <- data.frame(stage=c(0:max(pred_stage)), prop = as.vector(prop.table(table(pred_stage))))
  stage_plot <- ggplot2::ggplot(data=stagedf) + ggplot2::geom_bar(ggplot2::aes(x=stage, y=prop), stat="identity", fill="#3A68AB") + ggplot2::theme_bw() + ggplot2::xlab("Disease Stage")+ ggplot2::ylab("Proportion") + ggplot2::scale_x_continuous(breaks=c(0:max(pred_stage)))+ ggplot2::theme(text=ggplot2::element_text(size=14))
  return(stage_plot)

}



#' Plot likelihood ascent
#'
#' @param likes A data frame.
#' @return A ggplot object.


plot_likes <- function(likes=data.frame()){

  likes_plot <- ggplot2::ggplot(data=likes) + ggplot2::geom_line(ggplot2::aes(x=iter, y=like, color=factor(start))) + ggplot2::theme_bw() + ggplot2::ylab("Log-likelihood") + ggplot2::xlab("Iteration") + ggplot2::scale_color_discrete(name="Initialized Sequence") + ggplot2::xlim(0,max(likes$iter)) + ggplot2::theme(legend.position = "top", text = ggplot2::element_text(size=14))

  return(likes_plot)

}

