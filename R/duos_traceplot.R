#' Plot Traceplots
#'
#' Plots the traceplots of the paramters from \code{duos}.
#'
#'
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param parameters The group of paramters to plot (see details).
#' @param plots An option to plot paramters on the same plot or as individuals (see details).
#' @param burnin The desired burnin to discard from the results. By default it is 1 so that all iterations are plotted.
#'
#' @export
#' @importFrom ggforce facet_wrap_paginate
#' @importFrom dplyr group_by summarise %>%
#' @importFrom tidyr gather
#' @details
#'
#'
#' \strong{Options for} \code{parameters}
#'
#' There are two sets of paramters that can plotted in the traceplots: the cut-points and the proportion parameters.
#' \itemize{
#'     \item \code{"c"}: Plots the traceplots of the cut-points that are ordered and between 0 and 1 (DEFAULT).
#'     \item \code{"p"}: Plots the traceplots of the proportion paramters that sum to one.
#'   }
#'
#' \strong{Options for} \code{plots}
#'
#' There are several options on how to display the traceplots.
#' \itemize{
#'     \item \code{"all"}: Overlays all traceplots in one plot (works well for the cut-point parameters) (DEFAULT).
#'     \item \code{"indiv"}: Create a grid of traceplots where each plot contains a single paramter's traceplot. Four plots are included in a traceplot so multiple panels are created there are more than 4 paramters.
#'   }
#'
#'
#' @return A plot of traceplots.
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Uniform Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Uniform(0,1) distribution with 50 data points.
#' duos_unif <- duos(runif(50), k=4, MH_N=20000)
#'
#' #Plot the traceplots of the cut-points all in a single graph
#' duos_traceplot(duos_unif)
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(0.5,0.5) distribution with 300 data points.
#' duos_arcsin <- duos(rbeta(300, 0.5, 0.5), k=10, MH_N=20000)
#'
#' #Plot the traceplots for the cut-points as individual graphs
#' #Note: The plots are printed four at a time so multiple panels are printed
#' duos_traceplot(duos_arcsin, plots="indiv")
#'
#' #Plot the CDF with the credible intervals
#' duos_plot(duos_arcsin, type="cdf", cri=TRUE)
#'
#' ## --------------------------------------------------------------------------------
#' ## Bimodal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' #Sample 150 random uniforms
#' U =runif(150)
#' y = rep(NA,150)
#' #Sampling from the mixture
#' for(i in 1:150){
#'   if(U[i]<.3){
#'    y[i] = rnorm(1,0,1)
#'   }else {
#'    y[i] = rnorm(1,4,1)
#'   }
#' }
#' # First run 'duos' on data sampled from a bimodal distribution with 150 data points.
#' duos_bimodal <- duos(y, k=6, MH_N=20000)
#'
#' #Plot the individual traceplots for the proportion paramters
#' duos_traceplot(duos_bimodal, parameters="p", plots="indiv")
#'
#' #Plot the individual tracepltos for the cut-points with a burnin of 10,000
#' duos_traceplot(duos_bimodal, parameters="c", burnin=10000)
#'


duos_traceplot <- function(duos_output, parameters="c", plots="all", burnin=1){



  C_output <- duos_output$C
  P_output <- duos_output$P


  if(parameters=="c"){
    if(plots=="all"){
    C <- data.frame(C_output[burnin:nrow(C_output),])
    C$Iteration <- burnin:nrow(C_output)
    C_plot <- C %>% gather(Parameter, Simulation,-Iteration)
    C_plot$Parameter <- gsub("X", "", C_plot$Parameter)

    C_plot$Parameter <- as.factor(C_plot$Parameter)
    C_plot$Parameter <- factor(C_plot$Parameter, levels=c(1:ncol(C)))

    print(ggplot(C_plot)+
      geom_line(aes(Iteration, Simulation, group=Parameter, color=Parameter))+
      theme_bw()+theme(axis.title = element_text(size = 12)))
    }else if (plots=="indiv"){
      C <- data.frame(C_output[burnin:nrow(C_output),])
      C$Iteration <- burnin:nrow(C_output)
      C_plot <- C %>% gather(Parameter, Simulation,-Iteration)
      C_plot$Parameter <- as.factor(gsub("X", "", C_plot$Parameter))

      C_plot$Parameter <- as.factor(C_plot$Parameter)
      C_plot$Parameter <- factor(C_plot$Parameter, levels=c(1:ncol(C)))

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
      P <- data.frame(P_output[burnin:nrow(P_output),])
      P$Iteration <- burnin:nrow(P_output)
      P_plot <- P %>% gather(Parameter, Simulation,-Iteration)
      P_plot$Parameter <- gsub("X", "", P_plot$Parameter)

      P_plot$Parameter <- as.factor(P_plot$Parameter)
      P_plot$Parameter <- factor(P_plot$Parameter, levels=c(1:ncol(P)))

      print(ggplot(P_plot)+
        geom_line(aes(Iteration, Simulation, group=Parameter, color=Parameter))+
        theme_bw()+theme(axis.title = element_text(size = 12)))
    }else if (plots=="indiv"){
      P <- data.frame(P_output[burnin:nrow(P_output),])
      P$Iteration <- burnin:nrow(P_output)
      P_plot <- P %>% gather(Parameter, Simulation,-Iteration)
      P_plot$Parameter <- as.factor(gsub("X", "", P_plot$Parameter))

      P_plot$Parameter <- as.factor(P_plot$Parameter)
      P_plot$Parameter <- factor(P_plot$Parameter, levels=c(1:ncol(P)))

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
