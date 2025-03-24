#' Clinical states visualization
#'
#' @description
#' Generate plot of individual-level data arranged by clinical state frequency where clinical states are ordered from most frequent to least
#' frequent.
#'
#' @param data A data frame.
#' @param id.var Unique row identifier.
#' @param state.vars Variable names for clinical states.
#' @return A ggplot object.


plot_states <- function(data=data.frame(), id.var=character(), state.vars=vector()){

  #subset df to id.var and clinical states
  states <- data[,c(id.var, state.vars)]

  #arrange clinical states
  states <- dplyr::arrange(states,states[,state.vars])

  #initialize df to save props for each clinical state
  props <- data.frame(states=c(state.vars), prop=rep(NA,length(state.vars)))

  #calculate the proportion for each clinical state
  for(i in state.vars){
    props[props$states==i,2] <- prop.table(table(states[,i]))[2]
  }

  #arrange proportions by highest to lowest
  props <- dplyr::arrange(props, -prop)

  #factor to preserve arranged order
  props$states <- factor(props$states)

  #melt data by id.var to get long format for plotting
  long <- reshape2::melt(states, id.var=id.var)

  #arrange in long form
  long <- dplyr::arrange(long,value)

  #factor again to preserve order at id.var level
  long[,id.var] <- factor( long[,id.var], levels=states[,id.var])

  #factor states var to preserve state order
  long$variable <- factor(long$variable, levels = props$states)

  long$value <- factor(long$value)

  #generate plot
  state_plot <- ggplot2::ggplot(long, ggplot2::aes(long[,id.var],variable)) +
                ggplot2::geom_tile(ggplot2::aes(fill = value)) +
                ggplot2::scale_fill_manual(name="Levels", values = c("white", "#6099C6")) +
                ggplot2::xlab("Individual Participants") +
                ggplot2::ylab("Clinical State") +
                ggplot2::theme_bw()+
                ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", text=ggplot2::element_text(size=14))

  return(state_plot)

}


#' Disease class distribution
#' @description
#' Generate a bar plot of the disease class distribution using the vector containing the predicted
#' disease class for each individual after fitting TB-STATIS.
#'
#' @param pred_class A vector of predicted disease classes.
#' @return A ggplot object.


plot_class <- function(pred_class=vector()){

  #make df that has nrows = number of disease classes and 2 columns
  #col1 is disease class, col2 is proportion of sample assigned to disease class
  classdf <- data.frame(class=c(0:max(pred_class)), prop = as.vector(prop.table(table(pred_class))))

  #generate bar plot for disease class distribution
  class_plot <- ggplot2::ggplot(data=classdf) +
                ggplot2::geom_bar(ggplot2::aes(x=class, y=prop), stat="identity", fill="#3A68AB") +
                ggplot2::theme_bw() + ggplot2::xlab("Disease Severity Class")+
                ggplot2::ylab("Proportion") + ggplot2::scale_x_continuous(breaks=c(0:max(pred_class)))+
                ggplot2::theme(text=ggplot2::element_text(size=14))

  return(class_plot)

}



#' Likelihood ascent plot
#'
#' @description
#' Generate a line plot of the likelihood ascent from fitting TB-STATIS. The x-axis is the iteration number, the y-axis is the log-likelihood,
#' and each line represents the log-likelihood at each iteration for a single initialized sequence.
#'
#' @param likes A data frame with three columns containing the log-likelihood for each iteration and each initialization.
#' @return A ggplot object.


plot_likes <- function(likes=data.frame()){

  #generate line plot for likelihood ascent
  likes_plot <- ggplot2::ggplot(data=likes) +
                ggplot2::geom_line(ggplot2::aes(x=iter, y=like, color=factor(start))) +
                ggplot2::theme_bw() + ggplot2::ylab("Log-likelihood") +
                ggplot2::xlab("Iteration") +
                ggplot2::scale_color_discrete(name="Initialized Sequence") +
                ggplot2::xlim(0,max(likes$iter)) +
                ggplot2::theme(legend.position = "top", text = ggplot2::element_text(size=14))

  return(likes_plot)

}

