#' Plot the Prior vs. the Posterior
#'
#' Plots the histograms of simulations from the prior verses the posterior from \code{duos}.
#' 
#' @usage
#' duos_pp(duos_output, parameters = "c", burnin = NA)
#'
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param parameters The group of parameters to plot (see details).
#' @param burnin The desired burnin to discard from the results. By default, it is half the number of iterations.
#'
#' @export
#' @details
#' 
#' The results are designed to plot of a 3X3 grid. If there are more than 9 parameters, separate plots are created and printed and can be viewed by clicking the arrow through the results in the 'Plots' window.
#'
#' \strong{Options for} \code{parameters}
#'
#' There are two sets of parameters that can plotted in the histograms: the cut-points and the bin proportion parameters.
#' \itemize{
#'     \item \code{"c"}: Plots the histograms of the cut-points. (DEFAULT).
#'     \item \code{"p"}: Plots the histograms of the bin proportion parameters.
#'   }
#'
#' @return A plot of overlaid histograms.
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a beat(2,5) distribution with 150 data points.
#' y <- rbeta(150, 2, 5)
#' duos_beta <- duos(y)
#'
#' # Plot histograms of the priors vs. the posteriors of the cut-point parameters 
#' duos_pp(duos_beta)
#' 
#' # Plot histograms of the priors vs. the posteriors of the proportion parameters 
#' duos_pp(duos_beta, parameters = "p")

duos_pp <- function(duos_output, parameters="c", burnin=NA){
  
  # Get the cut-point parameters
  C_output <- duos_output$C
  # Get the bin proportion parameters
  P_output <- duos_output$P

  # Set burnin to default is no value was entered
  if(is.na(burnin)){
    burnin <- ceiling(nrow(C_output)/2)
  }
  
  # Plot for cut-point parameters
  if(parameters=="c"){
      #Create data set with burnin removed
      C <- data.frame(C_output[burnin:nrow(C_output),])
      # Stack the columns
      C_plot <- C %>% gather(Parameter, Simulation)
      # Create a variable called posterior
      C_plot$Source <- "Posterior"
      
      # Get the number of columns
      k <<- ncol(C)
      
      # Create data set to contain prior simulations
      C_prior <- data.frame(Source = rep("Prior", nrow(C)))
      
      # Function to simulate from cut-point prior
      duos_prior <- function(x){
        prior_sim <- runif(k, 0, 1)
        return(sort(prior_sim))
      }
      
      # Run function to get the same number of simulations as in the posterior data set
      C_prior <- cbind(C_prior,t(apply(C_prior, 1, duos_prior)))
      # Add the names to match the posterior data set
      names(C_prior)[2:ncol(C_prior)] <- names(C)
      
      # Stack the prior columns
      C_plot_prior <- C_prior %>% gather(Parameter, Simulation, -Source)
      
      # Combine the two data sets
      C_plot <- rbind(C_plot, C_plot_prior)
      
      # Remove the 'X' from the names
      C_plot$Parameter <- as.factor(gsub("X", "", C_plot$Parameter))
      
      # Convert to a factor and fix order
      C_plot$Parameter <- as.factor(C_plot$Parameter)
      C_plot$Parameter <- factor(C_plot$Parameter, levels=c(1:ncol(C)))
      
      # Create an index for the loop
      graph_index <- ceiling(ncol(C_output)/9)
      for (i in 1:graph_index) {
        
        print(ggplot(C_plot, aes(Simulation,fill=Source, color=Source))+
                geom_histogram(position = "identity", alpha = .5)+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i, scales="free"))
        
        
            }
  }else if(parameters=="p"){ # Repeat same code above except for the bin proportions
    P <- data.frame(P_output[burnin:nrow(P_output),])
    P_plot <- P %>% gather(Parameter, Simulation)
    P_plot$Source <- "Posterior"
    
    k_p <<- ncol(P)
    
    P_prior <- data.frame(Source = rep("Prior", nrow(P)))
    
    # Simulate from the Dirichleet distribution using a gamma distribution
    duos_prior <- function(x){
      prior_sim <- rgamma(k_p, 1, 1)
      prior_sum <- sum(prior_sim)
      
      return(prior_sim/prior_sum)
    }
    
    P_prior <- cbind(P_prior,t(apply(P_prior, 1, duos_prior)))
    names(P_prior)[2:ncol(P_prior)] <- names(P)
    
    P_plot_prior <- P_prior %>% gather(Parameter, Simulation, -Source)
    
    P_plot <- rbind(P_plot, P_plot_prior)
    
    P_plot$Parameter <- as.factor(gsub("X", "", P_plot$Parameter))
    
    P_plot$Parameter <- as.factor(P_plot$Parameter)
    P_plot$Parameter <- factor(P_plot$Parameter, levels=c(1:ncol(P)))
    
    graph_index <- ceiling(ncol(P_output)/9)
    for (i in 1:graph_index) {
      

      print(ggplot(P_plot, aes(Simulation,fill=Source, color=Source))+
              geom_histogram(position = "identity", alpha = .5)+
              theme_bw()+theme(axis.title = element_text(size = 12))+
              facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i, scales="free"))
      
      
    }
  }
}
