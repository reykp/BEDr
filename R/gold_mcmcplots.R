#' Plots a Variety of Convergence Diagnostic Plots
#'
#' Plots the convergence plots on the parameters from \code{gold}.
#'
#' @usage 
#' gold_mcmcplots(duos_output, type = "traceplot", plots = "all", burnin = 1)
#'
#' @param gold_output The list returned by \code{duos} containing the density estimate results.
#' @param type The type of convergent plot to create (see details).
#' @param npar The number of paramters to plot. 3X3 grids of tracepltos are created (DEFAULT is 9).
#' @param burnin The desired burnin to discard from the results. By default, it is 1 so that all iterations are plotted.
#'
#' @details
#'
#' \strong{Options for} \code{type}
#'
#' There are several options on which convergence plots to create.
#' \itemize{
#'     \item \code{"traceplot"}: Creates trace plots on the parameters (DEFAULT).
#'     \item \code{"acf"}: Creates autocorrelation plots using lag=1000. The burnin should be set to the desired value to discard. 
#'     \item \code{"rm"}: Creates running mean plots on the parameters.
#'   }
#'   
#' If \code{npar} is specified to be 9, 9 paramters (as equally spaced as possible) are chosen to plot.
#'
#' @return A plot of trace plots, acf plots, or running mean pltos.
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
#' # Plot the trace plots of the cut-points on a single graph (trace plot is the default)
#' duos_mcmcplots(duos_unif)
#' 
#' # Plot the acf plots of the cut-points on separate graphs
#' duos_mcmcplots(duos_unif, type = "acf", plots = "indiv", burnin = 10000)
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(0.5,0.5) distribution with 300 data points.
#' y <- rbeta(300, 0.5, 0.5)
#' duos_arcsin <- duos(y = y, k = 10, MH_N = 20000)
#'
#' #Plot the trace plots for the cut-points as individual graphs
#' #Note: The plots are printed six at a time so multiple panels are printed
#' duos_mcmcplots(duos_arcsin, plots = "indiv")
#' 
#' #Plot the trace plots for the bin proportions as individual graphs
#' #Note: The plots are printed six at a time so multiple panels are printed
#' duos_mcmcplots(duos_arcsin, parameters = "p", plots = "indiv")
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
#' # Plot the running mean plots for the cut-point parameters 
#' duos_mcmcplots(duos_bimodal, type = "rm", parameters = "c")
#'
#' # Plot the autocorrelation plots for the bin proportions with a burnin of 10,000
#' duos_mcmcplots(duos_bimodal, type = "acf", parameters = "p", burnin = 10000)


gold_mcmcplots <- function(gold_output, type = "traceplot",  npar = 9, burnin=NA){
  
  G_output <- gold_output$G
  
  
  burnin_na <- FALSE
  if(is.na(burnin)){
    burnin_na <- TRUE
    burnin <- 1
  }
  
  if(burnin>nrow(G_output)){
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
  
  
  npar_intervals <- floor(ncol(G_output)/npar)
  
  parameters <- NA
  parameters[1] <- 1
  for(i in 2:npar){
    if((parameters[{i-1}] +npar_intervals)<ncol(G_output)){
      parameters[i] <- parameters[{i-1}] +npar_intervals
    }
  }
  
  
  G <- data.frame(G_output[(burnin+1):nrow(G_output),])
  G <- G[,parameters]
  G$Iteration <- (burnin+1):nrow(G_output)
  
  
  
  G_plot <- G %>% gather(Parameter, Simulation,-Iteration)
  G_plot$Parameter <- as.factor(gsub("X", "", G_plot$Parameter))
  G_plot$Parameter <- as.numeric(as.character(G_plot$Parameter))
  

  
  if(type == "traceplot"){
    graph_index <- ceiling((ncol(G)-1)/9)
    for (i in 1:graph_index) {
      
      
      print(ggplot(G_plot)+
              geom_line(aes(Iteration, Simulation))+
              #geom_point(aes(Iteration, Simulation))+
              theme_bw()+theme(axis.title = element_text(size = 12))+
              facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i,scales="free"))
      
    }
  }else if (type == "acf"){
    
    if(burnin_na){
      burnin <- nrow(G_output)/2
    }
    # Create data set to plot
    acf_plot <- data.frame(0, 0, 0)
    names(acf_plot) <- c("Parameter", "Lag", "ACF")
    
    for(i in parameters){
      get_acf <- acf(G_output[burnin:nrow(G_output),i], plot = FALSE, lag = 1000)
      acf_plot_add <- data.frame(rep(i, 1001), with(get_acf, data.frame(lag, acf)))
      names(acf_plot_add) <- c("Parameter", "Lag", "ACF")
      
      acf_plot <- rbind(acf_plot, acf_plot_add)
    }
    
    # Remove the row with zeros
    acf_plot <- acf_plot[-1,]
    
    # Convert parameter to factor
    acf_plot$Parameter <- as.factor(acf_plot$Parameter)
    acf_plot$Parameter <- factor(acf_plot$Parameter, levels = c(parameters))
    
    graph_index <- ceiling((ncol(G)-1)/9)
    
    for (i in 1:graph_index) {
      
      print(ggplot(acf_plot,aes(x = Lag, y = ACF))+
              geom_hline(aes(yintercept = 0)) +
              geom_segment(mapping = aes(xend = Lag, yend = 0, color = Parameter))+
              theme_bw()+theme(axis.title = element_text(size = 12))+
              facet_wrap_paginate(~Parameter, ncol = 2, nrow = 3, page = i, scales="free"))
      
      
    }
      
    
 }else if (type %in% c("rm", "RM")){
    
      
      if(parameters %in% c("c", "C")){
      # Remove burnin
      if(burnin>1){
        C <- data.frame(C_output[{burnin+1}:nrow(C_output),])
      }else{
        C <- data.frame(C_output[{burnin}:nrow(C_output),])
        
      }
      
      # Calculate running means
      C_RM <- apply(C, 2, function(x) cumsum(x)/seq_along(x))
      
      # Convert to data frame
      C_RM <- data.frame(C_RM)
      # Add iteration
      if(burnin>1){
        C_RM$Iteration <- {burnin+1}:nrow(C_output)
      }else{
        C_RM$Iteration <- burnin:nrow(C_output)
      }
      
      
      # Stack all variables into one column
      C_plot <- C_RM %>% gather(Parameter, RunningMean,-Iteration)
      # Remove X from the name
      C_plot$Parameter <- gsub("X", "", C_plot$Parameter)
      
      # Order parameters correctly
      C_plot$Parameter <- as.factor(C_plot$Parameter)
      C_plot$Parameter <- factor(C_plot$Parameter, levels=c(1:ncol(C)))
      
      if(plots=="all"){  
      # Print plot
      print(ggplot(C_plot)+
              geom_line(aes(Iteration, RunningMean, group=Parameter, color=Parameter))+
              theme_bw()+theme(axis.title = element_text(size = 12))+ylab("Running mean"))
    }else if (plots=="indiv"){
      
      graph_index <- ceiling(ncol(C)/6)
      for (i in 1:graph_index) {
        
        print((ggplot(C_plot)+
                   geom_line(aes(Iteration, RunningMean, group=Parameter, color=Parameter))+
                   theme_bw()+theme(axis.title = element_text(size = 12))+ylab("Running mean")+
                facet_wrap_paginate(~Parameter, ncol = 2, nrow = 3, page = i, scales="free")))
        
        
      }
    }
  }else if(parameters %in% c("p", "P")){
    # Remove burnin
    if(burnin>1){
      P <- data.frame(P_output[{burnin+1}:nrow(P_output),])
    }else{
      P <- data.frame(P_output[{burnin}:nrow(P_output),])
      
    }
    
    # Calculate running means
    P_RM <- apply(P, 2, function(x) cumsum(x)/seq_along(x))
    
    # Convert to data frame
    P_RM <- data.frame(P_RM)
    # Add iteration
    if(burnin>1){
      P_RM$Iteration <- {burnin+1}:nrow(P_output)
    }else{
      P_RM$Iteration <- burnin:nrow(P_output)
    }
    
    
    # Stack all variables into one column
    P_plot <- P_RM %>% gather(Parameter, RunningMean,-Iteration)
    # Remove X from the name
    P_plot$Parameter <- gsub("X", "", P_plot$Parameter)
    
    # Order parameters correctly
    P_plot$Parameter <- as.factor(P_plot$Parameter)
    P_plot$Parameter <- factor(P_plot$Parameter, levels=c(1:ncol(P)))
    
    if(plots=="all"){  
      # Print plot
      print(ggplot(P_plot)+
              geom_line(aes(Iteration, RunningMean, group=Parameter, color=Parameter))+
              theme_bw()+theme(axis.title = element_text(size = 12))+ylab("Running mean"))
    }else if (plots=="indiv"){
      
      graph_index <- ceiling(ncol(P)/6)
      for (i in 1:graph_index) {
        
        print((ggplot(P_plot)+
                 geom_line(aes(Iteration, RunningMean, group=Parameter, color=Parameter))+
                 theme_bw()+theme(axis.title = element_text(size = 12))+ylab("Running mean")+
                 facet_wrap_paginate(~Parameter, ncol = 2, nrow = 3, page = i, scales="free")))
        
        
      }
    
      }
  }
  }
}
  
