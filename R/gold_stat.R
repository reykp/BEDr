#' Bayesina Density Estimate of PDF from DUOS.
#'
#' Calculates an estimate of the density based on the output from \code{duos} for an individual
#' or vector of values
#'
#' @param x A single value or vector of values to calculate the density at
#' @param duos_output The list returned by \code{duos}
#' @param burnin The desired burnin to discard from including in the estimate
#' @param scale If the data is not already between 0 or 1, this determines if the pdf is returned on the (0,1) transformed data or the original data scale
#'
#' @export
#' @useDynLib BEDr


gold_stat <- function(stat,gold_output, p=NA,burnin=NA,scale=FALSE){



  G <- gold_output[[1]]
  widths <- gold_output[[2]]

  if(is.na(burnin)){
    burnin <- nrow(G)/2
  }

  if(burnin>nrow(G)){
    stop("The specified burnin is greater than the number of iterations.")
  }


  if(!is.na(p[1])&!(stat%in%c("q", "quantile", "quant"))){
  stop("p specified, but quantile not requested as the statistic.")
  }

  y_orig <- gold_output[[4]]
  min_y <<- min(y_orig)
  max_y <<- max(y_orig)

  x_gold <- gold_output[[3]]

  G_burnin <- G[burnin:nrow(G),]

  #Exponentiate for numerator
  G_burnin_exp <-  exp(G_burnin)
  #Calculate normalizing constat
  nc <- (widths%*%t(G_burnin_exp))

  #Calcualte pdf at each iteration
  G_PDF <- G_burnin_exp
  for (i in 1:nrow(G_burnin_exp)){
    #Need W to get weighted average
    G_PDF[i,] <- (G_PDF[i,]*widths)/nc[i]
  }

  if((min_y<0) | (max_y>1)){
    mean_matrix <- (G_PDF %*% x_gold)*(max_y+.00001-(min_y-.00001))+(min_y-.00001)
    mean_matrix_scaled <- (G_PDF %*% x_gold)
    var_matrix <- ((G_PDF %*% (x_gold^2))-mean_matrix_scaled^2)*(max_y+.00001-(min_y-.00001))^2
  }else{
    mean_matrix <- (G_PDF %*% x_gold)
    var_matrix <- ((G_PDF %*% (x_gold^2))-mean_matrix^2)
  }


  if((min_y<0) | (max_y>1)){
    input <- (x_gold)*(max_y+.00001-(min_y-.00001))+(min_y-.00001);
    duos_CDF <- gold_cdf(input,gold_output,burnin)
  }else{
    input <- x_gold
    duos_CDF <- gold_cdf(x_gold,gold_output,burnin)
  }

  x_cdf <- duos_CDF[[4]]
  CDF <- duos_CDF[[2]]

  #Create empty list to contain quantiles
  quantiles <- list()

  for(i in 1:length(p)){

    which_close <- min(which(CDF>=p[i]))

    diff1 <- abs(p[i]-CDF[which_close])
    if(which_close>1){
    diff2 <- abs(p[i]-CDF[{which_close-1}])

    if(diff1<=diff2){
      quantiles[[i]] <- x_cdf[which_close]
    }else{
      quantiles[[i]] <- x_cdf[{which_close-1}]
    }

    }else{
      quantiles[[i]] <- x_cdf[which_close]
    }

  }

  if(stat=="mean"){
    return(mean(mean_matrix))
  }else if (stat=="var"){
    return(mean(var_matrix))
  }else{
    return(quantiles)
  }

}
