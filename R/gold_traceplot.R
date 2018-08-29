#' Plot Traceplots
#'
#' Plots the traceplots of the paramters from \code{gold}.
#'
#' @param gold_output The list returned by \code{gold} containing the density estimate results.
#' @param npar The number of paramters to plot. 3X3 grids of tracepltos are created (DEFAULT is 9).
#' @param burnin The desired burnin to discard from the results. By default it is 1 so that all iterations are plotted.
#'
#' @export
#' @details
#'
#' If \code{npar} is specified to be 9, 9 paramters (as equally spaced as possible) are chosen to plot.
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
#' ## Gamma Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a gamma(1, 3) distribution with 300 data points.
#' y <- rgamma(300, 1, 3)
#' gold_gamma <- gold(y, s1 = 2, c1 = 4, s2 = 0.6, c2 = 1, MH_N = 20000)
#'
#' #By default, 9 parameters are plotted
#' gold_traceplot(gold_gamma)
#'
#' #Plot the PDF with the credible intervals
#' gold_plot(gold_gamma, type="pdf", cri=TRUE)
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
#' gold_bimodal <- gold(y, s1 = 1, c1 = 5, s2 = 0.4, c2 = 4, MH_N = 20000)
#'
#' #Plot 15 of the paramters after discarding the first 5000
#' gold_traceplot(gold_bimodal, npar = 15, burnin = 5000)
#'
#' #Plot the CDF with the emprical CDF
#' gold_plot(gold_bimodal, type = "cdf", data = TRUE)



gold_traceplot <- function(gold_output,npar=9, burnin=1){



  G_output <- gold_output$G

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

  graph_index <- ceiling((ncol(G)-1)/9)
  for (i in 1:graph_index) {


        print(ggplot(G_plot)+
                geom_line(aes(Iteration, Simulation))+
                #geom_point(aes(Iteration, Simulation))+
                theme_bw()+theme(axis.title = element_text(size = 12))+
                facet_wrap_paginate(~Parameter, ncol = 3, nrow = 3, page = i,scales="free"))

  }

}
