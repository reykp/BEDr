#' Estimate of the CDF from \code{duos}
#'
#' Calculates a posterior estimate of the CDF based on the output from \code{duos} for an individual
#' or vector of values.
#'
#' @usage
#' duos_cdf(x, duos_output, burnin = NA, scale = FALSE)
#' 
#' @param x A single value or vector of values at which to calculate the CDF. These values should be entered on the scale of the data (i.e. values can fall outside of 0 and 1 if data does).
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param burnin The desired burnin to discard from the results. If no value is entered, the default is half the number of iterations.
#' @param scale This value TRUE/FALSE indicates whether to return scaled or unscaled results IF the original data does not fall between 0 and 1. The default is FALSE (i.e. returns results on the original data scale).
#'
#' @export
#'
#' @details
#'
#' The function \code{duos_cdf} returns the posterior mean CDF. The CDF is calculated based on the following equation at each iteration:
#'
#' \eqn{F(x) =}
#'
#' \deqn{(\pi_1) / (\gamma_1) * x , 0 \le x < \gamma_1}
#' \deqn{\pi_1 + (\pi_2) / (\gamma_2-\gamma_1) * (x-\gamma_1) , \gamma_1 \le x  < \gamma_2}
#' \deqn{\pi_1 + \pi_2 + (\pi_3) / (\gamma_3-\gamma_2) * (x-\gamma_2) , \gamma_2 \le x  < \gamma_3}
#' \deqn{...  , ... \le x  < ...}
#' \deqn{\sum_{i=1}^{k-1} \pi_i + (\pi_k) / (\gamma_k-\gamma_{k-1}) * (x-\gamma_{k-1}) , \gamma_{k-1} \le x  < \gamma_k}
#' \deqn{\sum_{i=1}^{k}\pi_i + (\pi_{k+1}) / (1-\gamma_k) * (x-\gamma_k) , \gamma_k \le x  < 1}
#'
#' where \eqn{\gamma_1 < \gamma_2 < ... < \gamma_k  is in (0,1) and \pi_1 + \pi_2 + ... + \pi_{k+1} = 1}
#' @return
#'
#' \code{duos_cdf} returns a list of the CDF results from \code{duos}.
#'
#' \item{\code{cdf}}{A vector of the posterior mean CDF at each value in \code{x}.}
#' \item{\code{cri}}{A matrix with 2 columns and a row for each value in \code{x}. It contains the 95\% credible interval for the CDF at each point in \code{x}.}
#' \item{\code{mat}}{A matrix containing the CDF values for \code{x} at EACH iteration after the burnin is discarded. There is a column for each value in \code{x}.}
#' \item{\code{x}}{A vector containing the values at which to estimate the CDF. }

#' @examples

#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(2,5) distribution with 100 data points.
#' y <- rbeta(100, 2, 5)
#' duos_beta <- duos(y = y, k = 5, MH_N = 20000)
#' # Estimate the CDF at a vector of values
#' cdf_beta <- duos_cdf(x = c(.01, .25, .6, .9), duos_beta)
#'
#' # Examine the CDF at each value in 'x'
#' cdf_beta$cdf
#'
#' # Examine the credible intervals of the CDF at 'x'
#' cdf_beta$cri
#'
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Normal(0,1) distribution with 50 data points.
#' y <- rnorm(50, 0, 1)
#' duos_norm <- duos(y = y, k = 4, MH_N = 20000, scale_l = sd(y), scale_u = sd(y))
#' # Estimate the CDF at a vector of values
#' cdf_norm <- duos_cdf(x = c(-2, -1, 0, .8, 1.8), duos_norm)
#'
#' # Examine the CDF at 'x'
#' cdf_norm$cdf
#'
#' # Examine the credible intervals of the CDF at 'x'
#' cdf_norm$cri
#'
#' # Plot a histogram of distribution of the posterior draws for the CDF estimate at 0.8
#' hist(cdf_norm$mat[, 4])
#' 
#' # Find probability of being greater than 2
#' 1 - duos_cdf(2, duos_norm)$cdf
#' 


duos_cdf <- function(x, duos_output, burnin=NA, estimate = "mean"){

  # Initial set up ##############
  
  # If burnin is NA, use half as default
  if(is.na(burnin)){
    burnin <- ceiling(nrow(duos_output$C)/2)
  }

  #Get the Cutpoints
  C <- duos_output$C

  #Get the bin proportions
  P <- duos_output$P
  
  # Error checks ##############
  
  # Do an error check to make sure burnin not greater than number of iterations
  if(burnin>nrow(C)){
    stop("The specified burnin is greater than the number of iterations.")
  }
  # Check if burnin is an integer
  if(burnin/ceiling(burnin) < 1){
    stop("burnin should be an integer greater than zero.")
  }
  
  # Check if burnin is an positive number
  if(burnin <= 0){
    stop("burnin should be an integer greater than zero.")
  }
  
  # Check to make sure can estimate cdf at x
  if(max(duos_output$y)>1 | min(duos_output$y)< 0){
    
    outside_range_u <- which(x > (max(duos_output$y)+duos_output$scale_u))
    outside_range_l <- which(x < (min(duos_output$y)-duos_output$scale_l))
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
  
  # Number of cutpoints
  k <<- ncol(C)

  # Calculate cdf at:
  input <<- x

  # Get the data
  y_orig <- duos_output$y
  
  # Calculate the maximum and minimum
  min_y <- min(y_orig)
  max_y <- max(y_orig)

  # Get the scale paramters
  scale_l <- duos_output$scale_l
  scale_u <- duos_output$scale_u
  
  # If the data is not between 0 and 1, scale 'x' to be on the same scale the density was estimated on
  if((min_y<0 | max_y>1)){
    input <<- (x-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l))
  }
  # Function that calculates CDF
  # input and k are global paramters
  cdf_forapply <- function(x){
    # Vector for CDF results
    pr <- rep(0,length(input))
    # Vector for cutpoints from rows with 0 and 1 added in
    c_full <- c(0, x[1:k],1)
    #Bin proportions
    p <- x[{k+1}:length(x)]
    # Loop to get CDF estimates
    for (j in 1:(k+1)){
      pr[which(input<c_full[(j+1)] & input>=c_full[j])]<- sum(p[1:{j-1}])*ifelse(j>1,1,0)+(p[j]/(c_full[{j+1}]-c_full[{j}]))*(input[which(input<c_full[{j+1}] & input>=c_full[j])]-c_full[{j}])
    }
    return(pr)
  }
  
  # Only use rows past burnin
  C_sub <- C[{burnin+1}:nrow(C),]
  P_sub <- P[{burnin+1}:nrow(C),]

  # Get matrix of cdf values
  if(length(x)==1){
  cdf_matrix <- matrix(apply(cbind(C_sub,P_sub), 1, cdf_forapply), nrow=1)
  }else{
    cdf_matrix <- apply(cbind(C_sub,P_sub), 1, cdf_forapply)
  }
  
  # Calculate mean at each 'x'
  if(estimate == "mean"){
    cdf_y <- apply(cdf_matrix, 1, mean)
  }else{
    cdf_y <- apply(cdf_matrix, 1, median)
  }
  
  # Calculate percentiles at each 'x'
  cdf_y_perc <- apply(cdf_matrix, 1, quantile, probs=c(.025, .975))

  #If data is outside (0,1) and scaling was not requested, return all results on data scale, else
  #return scaled values

  if((min_y<0 | max_y>1)){
    return(list(cdf=cdf_y, cri=t(cdf_y_perc), mat=t(cdf_matrix),x=x))
  }else{
    return(list(cdf=cdf_y, cri=t(cdf_y_perc), mat=t(cdf_matrix), x=input))
  }
}
