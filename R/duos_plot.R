#' Plot PDF or CDF
#'
#' Plots the posterior mean PDF or CDF based on the output from \code{duos}.
#'
#' @usage
#' duos_plot(duos_output, type = "pdf", burnin = NA, cri = FALSE, data = FALSE, interact = FALSE)
#' 
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param type The desired type of graph (see details).
#' @param burnin The desired burnin to discard from the results. If no value is entered, the default is half the number of iterations.
#' @param cri An option to include credible intervals. The default is FALSE.
#' @param data An option to include data in the graph (see details).
#' @param interact An option to make the plots interactive. The default is FALSE.
#'
#' @export
#' @importFrom ggplot2 ggplot aes labs theme_bw theme geom_histogram
#' geom_line expand_limits ylab xlab element_text geom_point
#' @importFrom plotly ggplotly
#' @importFrom gridExtra grid.arrange
#'
#' @details
#'
#' \strong{Options for} \code{type}
#'
#' The input from \code{duos} can be used to plot a Bayesian estimate of the PDF or the CDF.
#' \itemize{
#'     \item \code{"pdf"}: The density, f(x), is calculated at a grid of 999 points between 0 and 1 at each iteration. The density at each grid point is then averaged across the iterations to produce the results in the graph (DEFAULT).
#'     \item \code{"cdf"}: The CDF, F(x), is calculated at a grid of 999 points between 0 and 1 at each iteration. The CDF at each grid point is then averaged across the iterations to produce the results in the graph.
#'   }
#'
#' \strong{Options for} \code{cri}
#'
#' Credible intervals can also be added to the plot of the PDF or CDF. 
#' \itemize{
#'     \item \code{"FALSE"}: No credible intervals lines are plotted (DEFAULT).
#'     \item \code{"TRUE"}: Credible interval lines are plotted in red on the PDF or CDF. These are calculated by taking the 0.025th and 0.975th quantiles of the iterations from \code{duos} after burnin.
#'   }
#'
#' \strong{Options for} \code{data}
#'
#' Incorporates the data into the PDF and CDF.
#' \itemize{
#'     \item \code{"FALSE"}: The data is not included (DEFAULT).
#'     \item \code{"TRUE"}: If the PDF is plotted, a histogram is overlaid with the density estimate. If the CDF is plotted, the empirical CDF is overlaid with the \code{duos} CDF estimate.
#'   }
#'   
#' \strong{Options for} \code{interact}
#'
#' Provides the PDF or CDF plot in an interactive setting using \code{plotly}.
#' \itemize{
#'     \item \code{"FALSE"}: No interactivity (DEFAULT).
#'     \item \code{"TRUE"}: Allows for zooming in on locations and running mouse over plot to see PDF or CDF values.
#'   }
#'
#' @return A plot of the PDF or CDF estimate.
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Uniform Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Uniform(0,1) distribution with 50 data points.
#' y <- runif(50)
#' duos_unif <- duos(y = y, k = 4, MH_N = 20000)
#'
#' # Plot the PDF with all default values
#' duos_plot(duos_unif)
#' 
#' # Plot the PDF with the data overlayed
#' duos_plot(duos_unif, data = TRUE)
#'
#' # Plot the CDF with the credible intervals and the empirical CDF
#' duos_plot(duos_unif, type="cdf", cri=TRUE)
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(0.5,0.5) distribution with 300 data points.
#' y <- rbeta(300, 0.5, 0.5)
#' duos_arcsin <- duos(y, k = 10, MH_N = 20000)
#' 
#' # Plot the PDF with the data in an interactive setting
#' duos_plot(duos_arcsin, type = "pdf", data = TRUE, interact = TRUE)
#'
#' # Plot the CDF with the empirical cdf
#' duos_plot(duos_arcsin, type = "cdf", data = TRUE)
#'
#' ## --------------------------------------------------------------------------------
#' ## Bimodal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # Sample 150 random uniforms
#' u  <- runif(150)
#' y <- rep(NA,150)
#' # Sampling from the mixture
#' for(i in 1:150){
#'   if(u[i] < 0.3){
#'    y[i]  <-  rnorm(1, 0, 1)
#'   }else {
#'    y[i] <- rnorm(1, 4, 1)
#'   }
#' }
#' 
#' # First run 'duos' on data sampled from a bimodal distribution with 150 data points.
#' duos_bimodal <- duos(y, k = 8, MH_N = 20000, scale_l = 0.5*sd(y), scale_u = 0.5*sd(y))
#'
#' # Plot the PDF with the data and credible intervals
#' duos_plot(duos_bimodal, data = TRUE, cri = TRUE)
#'
#' #Plot the PDF with credible intervals and a histogram of the data
#' duos_plot(duos_bimodal, cri=TRUE, data=TRUE)
#'
#' # Plot the CDF interactively
#' duos_plot(duos_bimodal, type = "cdf", interact = TRUE)
#' 


duos_plot <- function(duos_output,type="pdf",burnin=NA, cri=FALSE, data=FALSE, interact=FALSE){
  
  if(!(type%in%c("pdf","Pdf", "PDF", "cdf", "Cdf", "CDF"))){
    stop("Please choose from the available graph types: 'pdf' or 'cdf'.")
  }
  
  # Plotting implemented in helper functions
  if(type=="pdf"){
    helper_duos_pdf_plot(duos_output, burnin, cri, data, interact)
  }else if (type=="cdf"){
    helper_duos_cdf_plot(duos_output, burnin, cri, data, interact)
  }
}
