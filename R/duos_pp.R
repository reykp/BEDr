#' Plot the Prior vs. the Posterior
#'
#' Plots the histograms and densities of simulations from the prior verses the posterior from \code{duos}.
#'
#'
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param parameters The group of paramters to plot (see details).
#' @param burnin The desired burnin to discard from the results. By default it is half the number of iterations.
#'
#' @export
#' @details
#'
#'
#' \strong{Options for} \code{parameters}
#'
#' There are two sets of paramters that can plotted in the histograms: the cut-points and the proportion parameters.
#' \itemize{
#'     \item \code{"c"}: Plots the histograms of the cut-points that are ordered and between 0 and 1 (DEFAULT).
#'     \item \code{"p"}: Plots the histograms of the proportion paramters that sum to one.
#'   }
#'
#'
#' @return A plot of histograms.
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a beat(2,5) distribution with 50 data points.
#' y <- rbeta(50, 2, 5)
#' duos_beta <- duos(y, k=4, MH_N=20000)
#'
#' #Plot histograms of the priors vs. the posteriors of the cut-point parameters 
#' duos_pp(duos_beta)
#' 
#' #Plot histograms of the priors vs. the posteriors of the proportion parameters 
#' duos_pp(duos_beta, parameters = "p")



duos_pp <- function(duos_output, parameters="c", burnin=NA){
  
  C_output <- duos_output$C
  P_output <- duos_output$P

  if(is.na(burnin)){
    burnin <- ceiling(nrow(C_output)/2)
  }
  
  
  if(parameters=="c"){
      C <- data.frame(C_output[burnin:nrow(C_output),])
      C_plot <- C %>% gather(Parameter, Simulation)
      C_plot$Source <- "Posterior"
      
      k <<- ncol(C)
      
      C_prior <- data.frame(Source = rep("Prior", nrow(C)))
      
      duos_prior <- function(x){
        prior_sim <- runif(k, 0, 1)
        return(sort(prior_sim))
      }
      
      C_prior <- cbind(C_prior,t(apply(C_prior, 1, duos_prior)))
      names(C_prior)[2:ncol(C_prior)] <- names(C)
      
      C_plot_prior <- C_prior %>% gather(Parameter, Simulation, -Source)
      
      C_plot <- rbind(C_plot, C_plot_prior)
      
      C_plot$Parameter <- as.factor(gsub("X", "", C_plot$Parameter))
      
      C_plot$Parameter <- as.factor(C_plot$Parameter)
      C_plot$Parameter <- factor(C_plot$Parameter, levels=c(1:ncol(C)))
      
      graph_index <- ceiling(ncol(C_output)/4)
      for (i in 1:graph_index) {
        
        print(ggplot(C_plot, aes(Simulation,fill=Source, color=Source))+
                geom_histogram(position = "identity", alpha = .5)+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i, scales="free"))
        
        
            }
  }else if(parameters=="p"){
    P <- data.frame(P_output[burnin:nrow(P_output),])
    P_plot <- P %>% gather(Parameter, Simulation)
    P_plot$Source <- "Posterior"
    
    k_p <<- ncol(P)
    
    P_prior <- data.frame(Source = rep("Prior", nrow(P)))
    
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
    
    graph_index <- ceiling(ncol(P_output)/4)
    for (i in 1:graph_index) {
      

      print(ggplot(P_plot, aes(Simulation,fill=Source, color=Source))+
              geom_histogram(position = "identity", alpha = .5)+
              theme_bw()+theme(axis.title = element_text(size = 12))+
              facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i, scales="free"))
      
      
    }
  }
}
