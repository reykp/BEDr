#' Bayesina Density Estimate of PDF from DUOS.
#'
#' Calculates an estimate of the density based on the output from \code{duos} for an individual
#' or vector of values
#'
#' @param x A single value or vector of values to calculate the density at
#' @param duos_output The list returned by \code{duos}
#' @param burnin The desired burnin to discard from including in the estimate
