#' Traceplots from DUOS.
#'
#' Traceplots
#'
#' @param gold_output The list returned by \code{duos}
#' @param type The type of desired graph (pdf, cdf)
#' @param burnin The desired burnin to discard from including in the estimate
#' @param cr_i An option to include credible intervals
#' @param data An option to overlay histogram of original data on plot
#'
#' @export


gold_traceplot <- function(gold_output,npar=9, burnin=NA, start=NA){



  G_output <- gold_output[[1]]

  npar_intervals <- floor(ncol(G_output)/npar)

  parameters <- NA
  parameters[1] <- 1
  for(i in 2:npar){
    if((parameters[{i-1}] +npar_intervals)<ncol(G_output)){
      parameters[i] <- parameters[{i-1}] +npar_intervals
    }
  }

  if(is.na(burnin)&is.na(start)){
    start <- nrow(G)/2
  }

  G <- data.frame(G_output[start:nrow(G_output),])
  G <- G[,parameters]
  G$Iteration <- start:nrow(G_output)



  G_plot <- G %>% gather(Parameter, Simulation,-Iteration)
  G_plot$Parameter <- as.factor(gsub("X", "", G_plot$Parameter))
  G_plot$Parameter <- as.numeric(as.character(G_plot$Parameter))

  graph_index <- ceiling((ncol(G)-1)/9)
  for (i in 1:graph_index) {

        print(ggplot(G_plot)+
                geom_line(aes(Iteration, Simulation))+
                #geom_point(aes(Iteration, Simulation))+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i,scales="free"))

  }

}
