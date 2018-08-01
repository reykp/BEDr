#' Plot Probability Density from DUOS.
#'
#' Plots the density estimate based on the output from DUOS.
#'
#' @param duos_output The list returned by \code{duos}
#' @param type The type of desired graph.
#' @param burnin The desired burnin to discard from including in the estimate
#' @param cr_i An option to include credible intervals
#' @param data An option to overlay histogram of original data on plot
#'
#' @export
#' @importFrom ggplot2 ggplot aes labs theme_bw theme geom_histogram
#' geom_line expand_limits ylab xlab element_text

duos_plot <- function(duos_output,type="pdf",burnin=NA, cr_i=FALSE, data=FALSE){

  if(type=="pdf"){
    helper_duos_pdf_plot(duos_output, burnin, cr_i, data)
  }else if (type=="cdf"){
    helper_duos_cdf_plot(duos_output, burnin, cr_i, data)
  }
}
