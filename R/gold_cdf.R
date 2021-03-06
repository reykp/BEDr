#' Estimate of CDF from \code{gold}
#'
#' Calculates a posterior estimate of the CDF based on the output from \code{gold} for an individual
#' or vector of values.
#'
#' @param x A single value or vector of values at which to calculate the CDF. These values are to be entered on the scale of the data (i.e. values can fall outside of 0 and 1).
#' @param gold_output The list returned by \code{gold} containing the density estimate results.
#' @param burnin The desired burnin to discard from the results. If no values is entered, the default is half the number of iterations.
#'
#' @export
#'
#' @details
#'
#' The function \code{gold_cdf} returns the posterior mean CDF. The CDF is calculated based on the following equation at each iteration:
#'
#' \deqn{F(x) =}
#' \deqn{\int_{0}^{x} exp(g(y)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' where g(x) is an unknown log density.
#'
#' Given that g(x) is unknown, the normalizing constant is estimated using a weighted average and the set of unknown paramters that recieve a prior is g(x) at a finite set of points.
#' These weights are also used in the estimate of the integral in the numerator using Riemann sums.
#'
#' @return
#'
#' \code{gold_cdf} returns a list of the CDF results from \code{gold}.
#'
#' \item{\code{cdf}}{A vector of the posterior mean CDF at each value in \code{x}.}
#' \item{\code{cri}}{A matrix with 2 columns and rows equaling the length of \code{x} containing the 95\% credible interval for the CDF at each of the points in \code{x}.}
#' \item{\code{mat}}{A matrix containing the CDF values for each \code{x} at EACH iteratation after the burnin is discarded. The number of columns is the length of \code{x}.}
#' \item{\code{x}}{A vector containing the values at which to estimate the CDF. }
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(2,5) distribution wiht 100 data points.
#' y <- rbeta(100, 2, 5)
#' gold_beta <- gold(y, s1 = 1, c1 = 2, s2 = 1, c2= 0.6, N = 20000, graves = FALSE)
#' # Calculate CDF at a variety of values
#' cdf_beta <- gold_cdf(x = c(.01, .25, .6, .9), gold_beta)
#'
#' # Examine the CDF at 'x'
#' cdf_beta$cdf
#'
#' # Examine the credibal intervals of the CDF at 'x'
#' cdf_beta$cri
#'
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'gold' on data sampled from a Normal(0,1) distribution with 50 data points.
#' y <- rnorm(50, 0, 1)
#' # Use all defaults: enter '0.8' as a point of interest in order to examine the histogram
#' # of its CDF estimate
#' gold_norm <- gold(y, poi = 0.8)
#' # Estimate the CDF at (-2, -1, 0, .8, 1.8)
#' cdf_norm <- gold_cdf(x=c(-2, -1, 0, .8, 1.8), gold_norm)
#'
#' # Examine the CDF at 'x'
#' cdf_norm$cdf
#'
#' # Examine the credibal intervals of the CDF at 'x'
#' cdf_norm$cri
#'
#' # Histogram of distribution of the CDF density estimate at 0.8
#' hist(cdf_norm$mat[, 4])


gold_cdf <- function(x, gold_output, burnin=NA){

  # error check 'x'
  if(length(which(is.na(x)))>0){
    stop("'x' contains at least one missing value.")
  }
  
  #Get paramters
  G <- gold_output$G
  
  #If burnin not specified, it is half the iterations
  if(is.na(burnin)){
    burnin <- nrow(G)/2
  }
  
  #Error check to see if burnin greater than the number of iterations
  if(burnin>nrow(G)){
    stop("The specified 'burnin' is greater than the number of iterations.")
  }
  
  if(burnin <= 0){
    stop("Please specify a positive integer for 'burnin'.")
  }
  
  if(burnin/ceiling(burnin)<1){
    stop("Please specify a positive integer for 'burnin'.")
  }
  
  #Get widths for estimate of normalizing constant
  widths <- gold_output$widths
  #Get points at which density estimated
  x_gold <- gold_output$x
  #Get original data set
  y_orig <- gold_output$y

  #These are the points at wich reqested to estimate cdf
  input <- x

  #Calculate minimum and maximum
  #Find min and max of data
  if(!is.null(gold_output[["poi"]])){
    max_y <- max(y_orig, gold_output$poi)
    min_y <- min(y_orig, gold_output$poi)
  }else{
    min_y <- min(y_orig)
    max_y <- max(y_orig)
  }

  scale_l <- gold_output$scale_l
  scale_u <- gold_output$scale_u
  
  
  # Check to make sure can estimate p at x
  if(max(gold_output$y)>1 | min(gold_output$y)< 0){
    
    outside_range_u <- which(x > (max(gold_output$y)+gold_output$scale_u))
    outside_range_l <- which(x < (min(gold_output$y)-gold_output$scale_l))
    if((length(outside_range_u)>0) &(length(outside_range_l)==0)){
      print(x[outside_range_u])
      stop("The requested 'x' vector contains the above value(s) outside the range of max(y)+scale_u.")
    }
    
    if((length(outside_range_l)>0) &(length(outside_range_u)==0)){
      print(x[outside_range_l])
      stop("The requested 'x' vector contains the above value(s) outside the range of min(y)-scale_l.")
    }
    
    if((length(outside_range_l)>0) &(length(outside_range_u)>0)){
      print(x[c(outside_range_l, outside_range_u)])
      stop("The requested 'x' vector contains the above value(s) outside the range of min(y)-scale_l and max(y)+scale_u.")
    }
    
    
  }else{
    outside_range_u <- which(x > 1)
    outside_range_l <- which(x < 0)
    if((length(outside_range_u)>0)|(length(outside_range_l)>0)){
      print(x[c(outside_range_u, outside_range_l)])
      stop("The requested 'x' vector contains the above value(s) outside the range of (0, 1).")
    }
    
  }
  
  
  #If min and max of y outside of range of 1, data needs to be standardized
  if((min_y<0 | max_y>1)){
    input_scaled <- (input-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l))
  }else{
    input_scaled <- input
  }

  
  #Collect the iterations after burnin
  G_burnin <- G[{burnin+1}:nrow(G),]

  #Exponentiate for numerator
  G_burnin_exp <- exp(G_burnin)
  #Calculate normalizing constant
  nc <- (widths%*%t(G_burnin_exp))

  #Calcualte first part of CDF
  G_CDF_start <- G_burnin_exp
  for (i in 1:nrow(G_burnin_exp)){
    G_CDF_start[i,] <- widths*G_CDF_start[i,]/nc[i]
  }

  #Need to cumulative sum to get values of CDF
  G_CDF <- t(apply(G_CDF_start, 1,cumsum))

  #Calculate posterior mean cdf
  cdf_y <- NA
  for (j in 1:ncol(G_CDF)){
    cdf_y[j] <- mean(G_CDF[,j])
  }


  #Calculate percentiles
  cdf_y_perc <- apply(G_CDF, 2, quantile, probs=c(.025, .975))

  #The following steps are how to handle values in 'x' outside the range of the
  #original data
  G_CDF_return <- matrix(0, nrow=nrow(G_CDF), ncol=length(input))
  cdf_y_return <- NA
  cdf_y_perc_return <- matrix(0, nrow=length(input), ncol=2)

  
  for(i in 1:length(input_scaled)){
    
    find_x <- which(x_gold==input_scaled[i]) 
    
    if((min_y<0|max_y>1) & (input[i]<(min_y-scale_l)|input[i]>(max_y+scale_u))){
      print(input[i])
      print("is out of the range of the data in 'y'.")
    }else if(length(find_x)>0){
      cdf_y_return[i] <- cdf_y[find_x]
      G_CDF_return[,i] <- G_CDF[,find_x]
      cdf_y_perc_return[i,] <- t(cdf_y_perc)[find_x,]
    }else{
      ss_indiv <- smooth.spline(x_gold, cdf_y)
      cdf_y_return[i] <- predict(ss_indiv, input_scaled[i])$y
      G_CDF_return[,i] <- rep(NA, nrow(G_CDF_return))
      ss_indiv_q1 <- smooth.spline(x_gold, t(cdf_y_perc)[,1])
      cdf_y_perc_return[i,1] <- predict(ss_indiv_q1, input_scaled[i])$y
      ss_indiv_q2 <- smooth.spline(x_gold, t(cdf_y_perc)[,2])
      cdf_y_perc_return[i,2] <- predict(ss_indiv_q2, input_scaled[i])$y
    }
  }
  
  
  
  # #Since the density is estimated at the grid and data points, only estimate
  # #for values in x is the data point or grid point nearest to them
  # x_cdfs <- NA
  # for(i in 1:length(input)){
  #   #Check if input is within appropriate range (if outside of 0 and 1, is outside range
  #   #of original data)
  #   if((input_scaled[i]<1)&(input_scaled[i]>0)){
  #     #Find data point or grid point nearest
  #   diff1 <- input_scaled[i]-x_gold[max(which(x_gold<=input_scaled[i]))]
  #   if(!is.na(x_gold[max(which(x_gold<=input_scaled[i]))+1])){
  #   diff2 <- x_gold[max(which(x_gold<=input_scaled[i]))+1]-input_scaled[i]
  #   }else{
  #     diff2 <- 0
  #   }
  # 
  #   if(diff1<=diff2){
  #     x_cdfs[i] <- max(which(x_gold<=input_scaled[i]))
  #     G_CDF_return[,i] <- G_CDF[,x_cdfs[i]]
  #     cdf_y_return[i] <- cdf_y[x_cdfs[i]]
  #     cdf_y_perc_return[i,] <- t(cdf_y_perc)[x_cdfs[i],]
  #   }else{
  #     x_cdfs[i] <- max(which(x_gold<=input_scaled[i]))+1
  #     G_CDF_return[,i] <- G_CDF[,x_cdfs[i]]
  #     cdf_y_return[i] <- cdf_y[x_cdfs[i]]
  #     cdf_y_perc_return[i,] <- t(cdf_y_perc)[x_cdfs[i],]
  #   }
  #   }else{
  #     #if 'x' outside of range of original data, CDF value is 0 or 1
  #     if(input[i]<min_y){
  #       G_CDF_return[,i] <- rep(0,nrow(G_CDF_return))
  #       cdf_y_return[i] <- 0
  #       cdf_y_perc_return[i,] <- c(0,0)
  #     }else{
  #       G_CDF_return[,i] <- rep(1,nrow(G_CDF_return))
  #       cdf_y_return[i] <- 1
  #       cdf_y_perc_return[i,] <- c(1,1)
  #     }
  #   }
  # }

  #return results with original 'x' or 'x' after scaling
  if((min_y<0 | max_y>1)){
    return(list(cdf=cdf_y_return, cri=cdf_y_perc_return, mat=G_CDF_return, x=input))

  }else{
    return(list(cdf=cdf_y_return, cri=cdf_y_perc_return, mat=G_CDF_return, x=input_scaled))
  }

}
