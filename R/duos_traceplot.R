#' Traceplots from DUOS.
#'
#' Traceplots
#'
#' @param duos_output The list returned by \code{duos}
#' @param type The type of desired graph (pdf, cdf)
#' @param burnin The desired burnin to discard from including in the estimate
#' @param cr_i An option to include credible intervals
#' @param data An option to overlay histogram of original data on plot
#'
#' @export
#' @importFrom ggplot2 ggplot aes labs theme_bw theme geom_histogram
#' geom_line expand_limits ylab xlab element_text
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#' @importFrom ggforce facet_wrap_paginate

duos_traceplot <- function(duos_output,parameters="c", plots="all", burnin=NA, start=NA){



  C_output <- duos_output[[1]]
  P_output <- duos_output[[2]]

  if(is.na(burnin)&is.na(start)){
    start <- nrow(C)/2
  }

  if(parameters=="c"){
    if(plots=="all"){
    C <- data.frame(C_output[start:nrow(C_output),])
    C$Iteration <- start:nrow(C_output)
    C_plot <- C %>% gather(Parameter, Simulation,-Iteration)
    C_plot$Parameter <- gsub("X", "", C_plot$Parameter)

    ggplot(C_plot)+
      geom_line(aes(Iteration, Simulation, group=Parameter, color=Parameter))+
      theme_bw()+theme(axis.title = element_text(size = 12))
    }else if (plots=="indiv"){
      C <- data.frame(C_output[start:nrow(C_output),])
      C$Iteration <- start:nrow(C_output)
      C_plot <- C %>% gather(Parameter, Simulation,-Iteration)
      C_plot$Parameter <- as.factor(gsub("X", "", C_plot$Parameter))

      graph_index <- ceiling(ncol(C_output)/4)
      for (i in 1:graph_index) {

        print(ggplot(C_plot)+
                geom_line(aes(Iteration, Simulation))+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 2, nrow = 2, page = i, scales="free"))


      }
    }
  }else if(parameters=="p"){
    if(plots=="all"){
      P <- data.frame(P_output[start:nrow(P_output),])
      P$Iteration <- start:nrow(P_output)
      P_plot <- P %>% gather(Parameter, Simulation,-Iteration)
      P_plot$Parameter <- gsub("X", "", P_plot$Parameter)

      ggplot(P_plot)+
        geom_line(aes(Iteration, Simulation, group=Parameter, color=Parameter))+
        theme_bw()+theme(axis.title = element_text(size = 12))
    }else if (plots=="indiv"){
      P <- data.frame(P_output[start:nrow(P_output),])
      P$Iteration <- start:nrow(P_output)
      P_plot <- P %>% gather(Parameter, Simulation,-Iteration)
      P_plot$Parameter <- as.factor(gsub("X", "", P_plot$Parameter))

      graph_index <- ceiling(ncol(P_output)/4)
      for (i in 1:graph_index) {

        print(ggplot(P_plot)+
                geom_line(aes(Iteration, Simulation))+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 2, nrow = 2, page = i, scales="free"))


      }
    }
  }
}
