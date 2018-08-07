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
  x_gold <- gold_output[[3]]
  y_orig <- gold_output[[4]]

  input <- x

  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if((min_y<0 | max_y>1)){
    input_scaled <- (input-(min_y-.00001))/(max_y+.00001-(min_y-.00001))
  }else{
    input_scaled <- input
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

  G_CDF_return <- matrix(0, nrow=nrow(G_CDF), ncol=length(input))
  cdf_y_return <- NA
  cdf_y_perc_return <- matrix(0, nrow=length(input), ncol=2)

  x_cdfs <- NA
  for(i in 1:length(input)){
    if((input[i]<=max_y)&(input[i]>=min_y)){
    diff1 <- input_scaled[i]-x_gold[max(which(x_gold<=input_scaled[i]))]
    if(!is.na(x_gold[max(which(x_gold<=input_scaled[i]))+1])){
    diff2 <- x_gold[max(which(x_gold<=input_scaled[i]))+1]-input_scaled[i]
    }else{
      diff2 <- 0
    }

    if(diff1<=diff2){
      x_cdfs[i] <- max(which(x_gold<=input_scaled[i]))
      G_CDF_return[,i] <- G_CDF[,x_cdfs[i]]
      cdf_y_return[i] <- cdf_y[x_cdfs[i]]
      cdf_y_perc_return[i,] <- t(cdf_y_perc)[x_cdfs[i],]
    }else{
      x_cdfs[i] <- max(which(x_gold<=input_scaled[i]))+1
      G_CDF_return[,i] <- G_CDF[,x_cdfs[i]]
      cdf_y_return[i] <- cdf_y[x_cdfs[i]]
      cdf_y_perc_return[i,] <- t(cdf_y_perc)[x_cdfs[i],]
    }
    }else{
      if(input[i]<min_y){
        G_CDF_return[,i] <- rep(0,nrow(G_CDF_return))
        cdf_y_return[i] <- 0
        cdf_y_perc_return[i,] <- c(0,0)
      }else{
        G_CDF_return[,i] <- rep(1,nrow(G_CDF_return))
        cdf_y_return[i] <- 1
        cdf_y_perc_return[i,] <- c(1,1)
      }
    }
  }


  if(scale==FALSE & (min_y<0 | max_y>1)){
    return(list(cdf_matrix=G_CDF_return, cdf_y=cdf_y_return, cdf_percentiles=cdf_y_perc_return, x=input))


  }else{
    return(list(cdf_matrix=G_CDF_return, cdf_y=cdf_y_return, cdf_percentiles=cdf_y_perc_return, x=input_scaled))
  }

}
