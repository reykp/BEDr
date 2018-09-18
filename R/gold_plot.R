#' Plot PDF or CDF
#'
#' Plots the posterior mean PDF or CDF based on the output from \code{gold}.
#'
#' @param gold_output The list returned by \code{gold} containing the density estimate results.
#' @param type The type of desired graph (see details).
#' @param burnin The desired burnin to discard from the results. If no values is entered, the default is half the number of iterations.
#' @param cri An option to include credible intervals.
#' @param data An option to include data in the graph (see details).
#' @param interact An option to make the plots interactive.
#'
#' @export
#'
#' @details
#'
#' The form of the density whose paramters are estimated in \code{gold} is below:
#'
#' \deqn{f(x) =}
#' \deqn{exp(g(x)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' where g(x) is an unknown log density.
#'
#' Given that g(x) is unknown, the normalizing constant is estimated using a weighted average and the set of unknown paramters that recieve a prior is g(x) at a finite set of points.
#'
#' \deqn{F(x) =}
#' \deqn{\int_{0}^{x} exp(g(y)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' where g(x) is an unknown log density.
#'
#' Given that g(x) is unknown, the normalizing constant is estimated using a weighted average and the set of unknown paramters that recieve a prior is g(x) at a finite set of points.
#' These weights are also used in the estimate of the integral in the numerator.
#' \strong{Options for} \code{type}
#'
#' The input from \code{gold} can be used to plot a Bayesian estimate of the PDF or the CDF.
#' \itemize{
#'     \item \code{"pdf"}: The density, f(x), is calculated at a grid overlayed with the data at each iteration to produce the results in the graph (DEFAULT).
#'     \item \code{"cdf"}: The CDF, F(x), is calculated at a grid overlayed with the data at each iteration. The CDF at each grid point and data point is then averaged across the iterations to produce the results in the graph.
#'   }
#'
#' \strong{Options for} \code{cri}
#'
#' Credibal intervals can also be added to the plot of the PDF or CDF.
#' \itemize{
#'     \item \code{"FALSE"}: No credible intervals lines are plotted (DEFAULT).
#'     \item \code{"TRUE"}: Credible interval lines are plotted in red on the PDF or CDF. These are calculated by taking the 0.025th and 0.975th quantiles of the iterations from \code{gold} after burnin.
#'   }
#'
#' \strong{Options for} \code{data}
#'
#' Incoproates the data into the PDF and CDF.
#' \itemize{
#'     \item \code{"FALSE"}: The data is not included (DEFAULT).
#'     \item \code{"TRUE"}: If the PDF is plotted, a histogram is overlayed with the density estimate. If the CDF is plotted, the empirical CDF is overlayed with the \code{gold} CDF estimate.
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
#' # First run 'gold' on data sampled from a Uniform(0,1) distribution with 50 data points.
#' y <- runif(50)
#' gold_unif <- gold(y, s1 = 1, c1 = 1, s2 = 0.9, c2 = 0.8, MH_N = 20000)
#'
#' #Plot the PDF with the data and credible intervals
#' gold_plot(gold_unif, type="pdf", data=TRUE)
#'
#' #Plot the CDF with the credible intervals and the empirical CDF
#' gold_plot(gold_unif, type="cdf", cri=TRUE, data=TRUE)
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'gold' on data sampled from a Beta(0.5,0.5) distribution with 300 data points.
#' y <- rbeta(300, 0.5, 0.5)
#' gold_arcsin <- gold(y, s1 = 2, c1 = 6, s2 = .2, c2 = 1.5, MH_N = 20000)
#' #Plot the PDF with the data
#' gold_plot(gold_arcsin, type="pdf", data=TRUE)
#'
#' #Plot the CDF with the credible intervals
#' gold_plot(gold_arcsin, type="cdf", cri=TRUE)
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
#' # First run 'gold' on data sampled from a bimodal distribution with 150 data points.
#' gold_bimodal <- gold(y, s1 = 1, c1 = 5, s2 = .4, c2 = 3, MH_N = 20000)
#'
#' #Plot the PDF
#' gold_plot(gold_bimodal)
#'
#' #Plot the PDF with credible intervals and a histogram of the data
#' gold_plot(gold_bimodal, cri=TRUE, data=TRUE)
#'
#' #Plot the CDF with the empirical CDF
#' gold_plot(gold_bimodal, type="cdf", data=TRUE)
#'


gold_plot <- function(gold_output,type="pdf",burnin=NA, cri=FALSE, data=FALSE, interact = FALSE, scale = FALSE){

  #Check if call is to plot pdf or cdf
  if(type=="pdf"){
    helper_gold_pdf_plot(gold_output, burnin, cri, data, interact, scale)
  }else if (type=="cdf"){
    helper_gold_cdf_plot(gold_output, burnin, cri, data, interact, scale)
  }
}
