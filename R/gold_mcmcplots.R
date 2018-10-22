#' Plots a Variety of Convergence Diagnostic Plots
#'
#' Plots the convergence plots on the parameters from \code{gold}.
#'
#' @usage 
#' gold_mcmcplots(gold_output, type = "traceplot", npar = 6, burnin = 1)
#'
#' @param gold_output The list returned by \code{duos} containing the density estimate results.
#' @param type The type of convergent plot to create (see details).
#' @param npar The number of paramters to plot. 3X3 grids of tracepltos are created (DEFAULT is 9).
#' @param burnin The desired burnin to discard from the results. By default, it is 1 so that all iterations are plotted.
#'
#' @export
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


gold_mcmcplots <- function(gold_output, type = "traceplot",  npar = 6, burnin=NA){
  
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
    graph_index <- ceiling((ncol(G)-1)/6)
    for (i in 1:graph_index) {
      
      
      print(ggplot(G_plot)+
              geom_line(aes(Iteration, Simulation))+
              #geom_point(aes(Iteration, Simulation))+
              theme_bw()+theme(axis.title = element_text(size = 12))+
              facet_wrap_paginate(~Parameter, ncol = 2, nrow = 3, page = i,scales="free")+
              theme(text = element_text(size=25)))
      
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
    
    graph_index <- ceiling((ncol(G)-1)/6)
    
    for (i in 1:graph_index) {
      
      print(ggplot(acf_plot,aes(x = Lag, y = ACF))+
              geom_hline(aes(yintercept = 0))+
              geom_segment(mapping = aes(xend = Lag, yend = 0))+
              theme_bw()+theme(axis.title = element_text(size = 12))+
              facet_wrap_paginate(~Parameter, ncol = 2, nrow = 3, page = i, scales="free")+
              theme(text = element_text(size=25)))
      
      
    }
      
    
 }else if (type %in% c("rm", "RM")){
    
      
      # Remove burnin
      if(burnin>1){
        G <- data.frame(G_output[{burnin+1}:nrow(G_output),])[,parameters]
      }else{
        G <- data.frame(G_output[{burnin}:nrow(G_output),])[,parameters]
        
      }
      
      # Calculate running means
      G_RM <- apply(G, 2, function(x) cumsum(x)/seq_along(x))
      
      # Convert to data frame
      G_RM <- data.frame(G_RM)
      # Add iteration
      if(burnin>1){
        G_RM$Iteration <- {burnin+1}:nrow(G_output)
      }else{
        G_RM$Iteration <- burnin:nrow(G_output)
      }
      
    
      # Stack all variables into one column
      G_plot <- G_RM %>% gather(Parameter, RunningMean,-Iteration)
      # Remove X from the name
      G_plot$Parameter <- gsub("X", "", G_plot$Parameter)
      
      # Order parameters correctly
      G_plot$Parameter <- as.factor(G_plot$Parameter)
      G_plot$Parameter <- factor(G_plot$Parameter, levels=c(parameters))
      
      graph_index <- ceiling(ncol(G)/6)
      for (i in 1:graph_index) {
        
        print((ggplot(G_plot)+
                   geom_line(aes(Iteration, RunningMean, group=Parameter))+
                   theme_bw()+theme(axis.title = element_text(size = 12))+ylab("Running mean")+
                facet_wrap_paginate(~Parameter, ncol = 2, nrow = 3, page = i, scales="free"))+
                theme(text = element_text(size=25)))
        
        
      }
    }
  }

