#' Plot the Prior vs. the Posterior
#'
#' Plots the histograms and densities of simulations from the prior verses the posterior from \code{gold}.
#'
#'
#' @param gold_output The list returned by \code{gold} containing the density estimate results.
#' @param npar The number of paramters to plot. 3X3 grids of tracepltos are created (DEFAULT is 9).
#' @param burnin The desired burnin to discard from the results. By default it is half the number of iterations.
#'
#' @export
#' @details
#'
#'
#' If \code{npar} is specified to be 9, 9 paramters (as equally spaced as possible) are chosen to plot.
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



gold_pp <- function(gold_output, npar=9, burnin=NA){
  
  G_output <- gold_output$G
  
  if(is.na(burnin)){
    burnin <- ceiling(nrow(G_output)/2)
  }
  
  npar_intervals <- floor(ncol(G_output)/npar)
  
  parameters <- NA
  parameters[1] <- 1
  for(i in 2:npar){
    if((parameters[{i-1}] +npar_intervals)<ncol(G_output)){
      parameters[i] <- parameters[{i-1}] +npar_intervals
    }
  }
  parameters <<- parameters
  
  G <- data.frame(G_output[(burnin+1):nrow(G_output),parameters])

  names(G) <- parameters
  G_plot <- G %>% gather(Parameter, Simulation)
  G_plot$Parameter <- as.numeric(as.character(G_plot$Parameter))
  G_plot$Source <- "Posterior"
  
  s1 <- gold_output$prior[which(names(gold_output$prior)=="s1")]
  c1 <- gold_output$prior[which(names(gold_output$prior)=="c1")]
  s2 <- gold_output$prior[which(names(gold_output$prior)=="s2")]
  c2 <- gold_output$prior[which(names(gold_output$prior)=="c2")]
  x <- gold_output$x
  
  cov_prior <- matrix(0, nrow=length(x), ncol=length(x))
  for (i in 1:length(x)){
    for(j in 1:length(x)){
      d <- abs(x[i]-x[j])
      if(c1*d<.5){
        cov_prior[i,j] <- s1^2*(1-6*(c1*d)^2+6*(c1*d)^3)
      }else if (c1*d <1){
        cov_prior[i,j] <- s1^2*(2*(1-c1*d)^3)
        
      } else {
        cov_prior[i,j] <- 0
        
      }
    }
  }
  
  cov_prior <- cov_prior[parameters, parameters]
  
      n_x <<- length(gold_output$x[parameters])

      n <- 1
      mu <- rep(0, n_x)
      tol = 1e-06
      p <- length(mu)
      eS <- svd(cov_prior)
      ev <- eS$d
      if (!all(ev >= -tol * abs(ev[1L])))
        stop("'Sigma' is not positive definite")
      
      mvrnorm_calc <<- eS$u %*% diag(sqrt(pmax(ev, 0)), p)

            
      G_prior <- data.frame(Source = rep("Prior", nrow(G)))
      
      gold_prior <- function(x){
        rw <- matrix(rnorm(n_x * 1), 1)
        prior_sim <- mvrnorm_calc %*% t(rw)
        
        return(prior_sim)
      }
      
      G_prior <- cbind(G_prior,t(apply(G_prior, 1, gold_prior)))
      names(G_prior)[2:ncol(G_prior)] <- names(G)
      
      G_plot_prior <- G_prior %>% gather(Parameter, Simulation, -Source)
      
      G_plot <- rbind(G_plot, G_plot_prior)
      
      #G_plot$Parameter <- as.factor(gsub("X", "", G_plot$Parameter))
      
      G_plot$Parameter <- as.factor(parameters)
      G_plot$Parameter <- factor(G_plot$Parameter, levels=c(parameters))
      
      graph_index <- ceiling(npar/9)
      for (i in 1:graph_index) {
        
        print(ggplot(G_plot, aes(Simulation,fill=Source, color=Source))+
                geom_histogram(position = "identity", alpha = .5, bins = 60)+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i))
        
        
            }
}
