#' Bayesina Density Estimate of CDF from GOLD.
#'
#' Calculates an estimate of the CDF based on the output from \code{duos} for an individual
#' or vector of values
#'
#' @param x A single value or vector of values to calculate the density at
#' @param gold_output The list returned by \code{duos}
#' @param burnin The desired burnin to discard from including in the estimate
#' @param scale If the data is not already between 0 or 1, this determines if the pdf is returned on the (0,1) transformed data or the original data scale
#'
#' @export

gold_cdf <- function(x, gold_output, burnin=NA,scale=FALSE){


  G <- gold_output[[1]]
  widths <- gold_output[[2]]
  x <- gold_output[[3]]
  y_orig <- gold_output[[4]]

  input_scaled <- x

  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if((min_y<0 | max_y>1)){
    input <- x*(max_y+.00001-(min_y-.00001))+(min_y-.00001)
  }


  if(is.na(burnin)){
    burnin <- nrow(G)/2
  }

  if(burnin>nrow(G)){
    stop("The specified burnin is greater than the number of iterations.")
  }

  G_burnin <- G[burnin:nrow(G),]

  #Exponentiate for numerator
  G_burnin_exp <- exp(G_burnin)
  #Calculate normalizing constat
  nc <- (widths%*%t(G_burnin_exp))

  #Calcualte pdf at each iteration
  G_CDF_start <- G_burnin_exp
  for (i in 1:nrow(G_burnin_exp)){
    G_CDF_start[i,] <- widths*G_CDF_start[i,]/nc[i]
  }


  G_CDF <- t(apply(G_CDF_start, 1,cumsum))


  #Calculate posterior mean pdf
  cdf_y <- NA
  for (j in 1:ncol(G_CDF)){
    cdf_y[j] <- mean(G_CDF[,j])
  }

  #Calculate percentiles

  cdf_y_perc <- apply(G_CDF, 2, quantile, probs=c(.025, .975))


  if(scale==FALSE & (min_y<0 | max_y>1)){
    return(list(cdf_matrix=G_CDF, cdf_y=cdf_y, cdf_percentiles=t(cdf_y_perc), x=input))


  }else{
    return(list(cdf_matrix=cdf_matrix, cdf_y=cdf_y, cdf_percentiles=t(cdf_y_perc), x=input_scaled))
}

}
