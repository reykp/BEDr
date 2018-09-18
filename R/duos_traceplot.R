#' Plot Trace Plots
#'
#' Plots the trace plots of the parameters from \code{duos}.
#'
#' @usage 
#' duos_traceplot(duos_output, parameters = "c", plots = "all", burnin = 1)
#'
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param parameters The group of parameters to plot (see details).
#' @param plots An option to plot parameters on the same plot or as individuals in a grid (see details).
#' @param burnin The desired burnin to discard from the results. By default, it is 1 so that all iterations are plotted.
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
#' There are two sets of parameters that can plotted in the trace plots: the cut-points and the bin proportion parameters.
#' \itemize{
#'     \item \code{"c"}: Plots the trace plots of the cut-points that are ordered and between 0 and 1 (DEFAULT).
#'     \item \code{"p"}: Plots the trace plots of the proportion parameters that sum to one.
#'   }
#'
#' \strong{Options for} \code{plots}
#'
#' There are several options on how to display the trace plots.
#' \itemize{
#'     \item \code{"all"}: Overlays all trace plots in one plot (works well for the cut-point parameters) (DEFAULT).
#'     \item \code{"indiv"}: Create a grid of trace plots where each plot contains a single parameter's trace plot. Nine plots are allowed in a grid so multiple graphs are created if there are more than 9 parameters.
#'   }
#'
#'
#' @return A plot of trace plots.
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Uniform Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Uniform(0,1) distribution with 50 data points.
#' y <- runif(50)
#' duos_unif <- duos(y = y)
#'
#' # Plot the trace plots of the cut-points on a single graph
#' duos_traceplot(duos_unif)
#' # Plot the trace plots of the cut-points on separate graphs
#' duos_traceplot(duos_unif, plots = "indiv")
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(0.5,0.5) distribution with 300 data points.
#' y <- rbeta(300, 0.5, 0.5)
#' duos_arcsin <- duos(y = y, k = 10, MH_N = 20000)
#'
#' #Plot the trace plots for the bin proportions as individual graphs
#' #Note: The plots are printed nine at a time so multiple panels are printed
#' duos_traceplot(duos_arcsin, parameters = "p", plots = "indiv")
#'
#' ## --------------------------------------------------------------------------------
#' ## Bimodal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # Sample 150 random uniforms
#' u <- runif(150)
#' y <- rep(NA, 150)
#' # Sampling from the mixture
#' for(i in 1:150){
#'   if(u[i]<.3){
#'    y[i] <- rnorm(1, 0, 1)
#'   }else {
#'    y[i] <- rnorm(1, 4, 1)
#'   }
#' }
#' # First run 'duos' on data sampled from a bimodal distribution with 150 data points.
#' duos_bimodal <- duos(y = y, k = 8, MH_N = 20000)
#'
#' # Plot the individual trace plots for the bin proportion parameters
#' duos_traceplot(duos_bimodal, parameters = "p", plots = "indiv")
#'
#' # Plot the individual trace plots for the cut-points parameters with a burnin of 10,000
#' duos_traceplot(duos_bimodal, parameters = "c", burnin = 10000)


duos_traceplot <- function(duos_output, parameters="c", plots="all", burnin=1){
  
  # Get cut-point parameters
  C_output <- duos_output$C
  # Get bin proportion parameters
  P_output <- duos_output$P

  
  if(burnin>nrow(C_output)){
    stop("The specified burnin is greater than the number of iterations.")
  }
  
  # Check if burnin is an integer
  if(burnin/ceiling(burnin) < 1){
    stop("burnin should be an integer greater than zero.")
  }
  
  # Check if burnin is an positive number
  if(burnin <= 0){
    stop("burnin should be an integer greater than zero.")
  }
  
  # Plot for cut-points
  if(parameters=="c"){
    # Single plot
    if(plots=="all"){
    
    # Remove burnin
    if(burnin>1){
    C <- data.frame(C_output[{burnin+1}:nrow(C_output),])
    # Create iteration variable for x-axs
    C$Iteration <- {burnin+1}:nrow(C_output)
    }else{
      C <- data.frame(C_output[{burnin}:nrow(C_output),])
      # Create iteration variable for x-axs
      C$Iteration <- burnin:nrow(C_output)
      
    }
    # Stack all variables into one column
    C_plot <- C %>% gather(Parameter, Simulation,-Iteration)
    # Remove X from the name
    C_plot$Parameter <- gsub("X", "", C_plot$Parameter)

    # Order parameters correctly
    C_plot$Parameter <- as.factor(C_plot$Parameter)
    C_plot$Parameter <- factor(C_plot$Parameter, levels=c(1:ncol(C)))

    # Print plot
    print(ggplot(C_plot)+
      geom_line(aes(Iteration, Simulation, group=Parameter, color=Parameter))+
      theme_bw()+theme(axis.title = element_text(size = 12)))
    }else if (plots=="indiv"){
      # Same as above with a facet_wrap
      C <- data.frame(C_output[burnin:nrow(C_output),])
      C$Iteration <- burnin:nrow(C_output)
      C_plot <- C %>% gather(Parameter, Simulation,-Iteration)
      C_plot$Parameter <- as.factor(gsub("X", "", C_plot$Parameter))

      C_plot$Parameter <- as.factor(C_plot$Parameter)
      C_plot$Parameter <- factor(C_plot$Parameter, levels=c(1:ncol(C)))

      graph_index <- ceiling(ncol(C_output)/9)
      for (i in 1:graph_index) {

        print(ggplot(C_plot)+
                geom_line(aes(Iteration, Simulation))+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i, scales="free"))


      }
    }
  }else if(parameters=="p"){
    # Repeate the  same process for P
    if(plots=="all"){
      if(burnin>1){
      P <- data.frame(P_output[{burnin+1}:nrow(P_output),])
      P$Iteration <- {burnin+1}:nrow(P_output)
      }else{
        P <- data.frame(P_output[burnin:nrow(P_output),])
        P$Iteration <- burnin:nrow(P_output)
      }
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

      graph_index <- ceiling(ncol(P_output)/9)
      for (i in 1:graph_index) {

        print(ggplot(P_plot)+
                geom_line(aes(Iteration, Simulation))+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i, scales="free"))


      }
    }
  }
}
