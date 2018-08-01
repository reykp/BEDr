#' Bayesina Density Estimate of PDF from GOLD.
#'
#' Calculates an estimate of the density based on the output from \code{duos} for an individual
#' or vector of values
#'
#' @param x A single value or vector of values to calculate the density at
#' @param gold_output The list returned by \code{duos}
#' @param burnin The desired burnin to discard from including in the estimate
#' @param scale If the data is not already between 0 or 1, this determines if the pdf is returned on the (0,1) transformed data or the original data scale
#'
#' @export

gold_pdf <- function(x, gold_output, burnin=NA,scale=FALSE){


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
  G_PDF <- G_burnin_exp
  for (i in 1:nrow(G_burnin_exp)){
    G_PDF[i,] <- G_PDF[i,]/nc[i]
  }

  #Calculate posterior mean pdf
  pdf_y <- NA
  for (j in 1:ncol(G_PDF)){
    pdf_y[j] <- mean(G_PDF[,j])
  }

  #Calculate percentiles

  pdf_y_perc <- apply(G_PDF, 2, quantile, probs=c(.025, .975))


  if(scale==FALSE & (min_y<0 | max_y>1)){
    pdf_y <- pdf_y/(max_y+.00001-(min_y-.00001))
    pdf_y_perc[1,] <- pdf_y_perc[1,]/(max_y+.00001-(min_y-.00001))
    pdf_y_perc[2,] <- pdf_y_perc[2,]/(max_y+.00001-(min_y-.00001))

    return(list(pdf_matrix=G_PDF, pdf_y=pdf_y, pdf_percentiles=t(pdf_y_perc), x=input))

  }else{
    return(list(pdf_matrix=G_PDF, pdf_y=pdf_y, pdf_percentiles=t(pdf_y_perc), x=input_scaled))
  }


}
