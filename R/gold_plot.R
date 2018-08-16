#' Plot Probability Density from GOLD.
#'
#' Plots the density estimate based on the output from gold.
#'
#' @param gold_output The list returned by \code{gold}
#' @param type The type of desired graph.
#' @param burnin The desired burnin to discard from including in the estimate
#' @param cri An option to include credible intervals
#' @param data An option to overlay histogram of original data on plot
#'
#' @export
#' @importFrom ggplot2 ggplot aes labs theme_bw theme geom_histogram
#' geom_line expand_limits ylab xlab element_text geom_point

gold_plot <- function(gold_output,type="pdf",burnin=NA, cri=FALSE, data=FALSE){

  if(type=="pdf"){
    helper_gold_pdf_plot(gold_output, burnin, cri, data)
  }else if (type=="cdf"){
    helper_gold_cdf_plot(gold_output, burnin, cri, data)
  }
}
