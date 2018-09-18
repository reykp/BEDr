#' Estimate of PDF from \code{duos}
#'
#' Calculates a posterior estimate of the PDF based on the output from \code{duos} for an individual
#' or vector of values.
#' 
#' @usage
#' duos_pdf(x, duos_output, burnin = NA, scale = FALSE)
#' 
#' @param x A single value or vector of values at which to calculate the PDF. These values should be entered on the scale of the data (i.e. values can fall outside of 0 and 1 if data does).
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param burnin The desired burnin to discard from the results. If no value is entered, the default is half the number of iterations.
#' @param scale This value TRUE/FALSE indicates whether to return scaled or unscaled results IF the original data does not fall between 0 and 1. The default is FALSE (i.e. returns results on the original data scale).
#'
#' @export
#'
#' @details
#'
#' The function \code{duos_pdf} returns the posterior mean PDF. The PDF is calculated based on the following equation at each iteration:
#'
#' \eqn{f(x) =}
#'
#' \deqn{(\pi_1) / (\gamma_1) , 0 \le x < \gamma_1}
#' \deqn{(\pi_2) / (\gamma_2-\gamma_1) , \gamma_1 \le x  < \gamma_2}
#' \deqn{(\pi_3) / (\gamma_3-\gamma_2) , \gamma_2 \le x  < \gamma_3}
#' \deqn{...  , ... \le x  < ...}
#' \deqn{(\pi_k) / (\gamma_k-\gamma_{k-1}) , \gamma_{k-1} \le x  < \gamma_k}
#' \deqn{(\pi_{k+1}) / (1-\gamma_k) , \gamma_k \le x  < 1}
#'
#' where \eqn{\gamma_1 < \gamma_2 < ... < \gamma_k  is in (0,1) and \pi_1 + \pi_2 + ... + \pi_{k+1} = 1}
#'
#' @return \code{duos_pdf} returns a list of the PDF results from \code{duos}.
#'
#' \item{\code{pdf}}{A vector of the posterior mean PDF at each value in \code{x}.}
#' \item{\code{cri}}{A matrix with 2 columns and a row for each value in \code{x}. It contains the 95\% credible interval for the CDF at each point in \code{x}.}
#' \item{\code{mat}}{A matrix containing the PDF values for \code{x} at EACH iteration after the burnin is discarded. There is a column for each value in \code{x}.}
#' \item{\code{x}}{A vector containing the values at which to estimate the PDF. If the data is not between 0 and 1 and scale=TRUE, the scaled version of \code{x} is returned.}

#' @examples

#' ## --------------------------------------------------------------------------------
#' ## Uniform Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run \code{duos} on data sampled from a Unif(0,1) distribution with 70 data points.
#' y <- runif(70)
#' duos_unif <- duos(y = y, k = 4, MH_N = 20000)
#' 
#' # Estimate the PDF at a vector of values
#' pdf_unif <- duos_pdf(x = c(.1, .5, .65, .98), duos_unif)
#'
#' # Examine the PDF at \code{x}
#' pdf_unif$pdf
#'
#' # Examine the credible intervals of the PDF at \code{x}
#' pdf_unif$cri
#'
#' ## --------------------------------------------------------------------------------
#' ## Gamma Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run \code{duos} on data sampled from a Gamma(2,2) distribution with 90 data points.
#' y <- rgamma(90, 2, 2)
#' duos_gamma <- duos(y = y, MH_N = 20000)
#' 
#' # Estimate the PDF at a vector of values
#' pdf_gamma <- duos_pdf(x = c(0.4, 1, 2, 3), duos_gamma)
#'
#' # Examine the PDF at \code{x}
#' pdf_gamma$pdf
#'
#' # Examine the credible intervals of the PDF at \code{x}
#' pdf_gamma$cri
#'
#' #Plot a histogram of distribution of the posterior draws for the CDF estimate at 1
#' hist(pdf_gamma$mat[,2])
#' 
#' # Data is scaled between 0 and 1 for the density estimation.
#' # If \code{scale} = TRUE is requested, the \code{x} is returned on the scaled range.
#' duos_pdf(x = c(0.4, 1, 2, 3), duos_gamma, scale = TRUE)$x1

duos_pdf <- function(x, duos_output, burnin=NA,scale=FALSE){

  # Get the cut-point matrix
  C <- duos_output$C
  
  
  # Set burnin to half the iterations if it is set to NA
  if(is.na(burnin)){
    burnin <- ceiling(nrow(C)/2)
  }
  
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
  
  # Check to make sure can estimate p at x
  if(max(duos_output$y)>1 | min(duos_output$y)< 0){
  
  outside_range_u <- which(x > (max(duos_output$y)+duos_output$scale_u))
  if(length(outside_range_u)>0){
    print(x[outside_range_u])
    stop("The requested 'x' vector contains a value outside the range of max(y)+scale_u.")
  }
  
  outside_range_l <- which(x < (min(duos_output$y)-duos_output$scale_l))
  if(length(outside_range_l)>0){
    print(x[outside_range_l])
    stop("The requested 'x' vector contains a value outside the range of min(y)-scale_l.")
  }
  }else{
    outside_range_u <- which(x > 1)
    if(length(outside_range_u)>0){
      print(x[outside_range_u])
      stop("The requested 'x' vector contains a value outside the range of (0, 1).")
    }
    
    outside_range_l <- which(x < 0)
    if(length(outside_range_l)>0){
      print(x[outside_range_l])
      stop("The requested 'x' vector contains a value outside the range of (0, 1).")
    }
  }
  
  
  # Bin proportions
  P <- duos_output$P

  # Get the number of cut-points
  k <<- ncol(C)
  
  # Calculate pdf at:
  input <<- x

  # Get data
  y_orig <- duos_output$y
  # Find max and min of data
  min_y <- min(y_orig)
  max_y <- max(y_orig)
  
  # Get scale parameters
  scale_l <- duos_output$scale_l
  scale_u <- duos_output$scale_u
  
  # Scale 'x' to be between 0 and 1 if necessary
  if((min_y<0 | max_y>1)){
    input <<- (x-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l))
  }else{
    if(scale == TRUE){
      stop("Scaling requested when no scaling was used.")
    }
  }

  # Function to calculate pdf
  pdf_forapply <- function(x){
    #Vector for probabilities
    pr <- rep(0,length(input))
    #Vector for cutpoints from rows with 0 and 1 added in
    c_full <- c(0, x[1:k],1)
    #Bin probabilities
    p <- x[{k+1}:length(x)]
    #Loop through density estimates
    for (j in 1:(k+1)){
      pr[which(input<c_full[(j+1)] & input>=c_full[j])]<- p[j]/(c_full[(j+1)]-c_full[j])
    }
    return(pr)
  }

  # Take burnin into account
  C_sub <- C[{burnin+1}:nrow(C),]
  P_sub <- P[{burnin+1}:nrow(C),]

  # Get matrix of pdf values
  if(length(x)==1){
    pdf_matrix <- matrix(apply(cbind(C_sub,P_sub), 1, pdf_forapply), nrow=1)
  }else{
    pdf_matrix <- apply(cbind(C_sub,P_sub), 1, pdf_forapply)
  }
  # Average to get mean estimates
  pdf_y <- apply(pdf_matrix, 1, mean)
  # Find credibal intervals
  pdf_y_perc <- apply(pdf_matrix, 1, quantile, probs=c(.025, .975))

  # PDF needs to be scaled correctly if data not between 0 and 1
  if(scale==FALSE & (min_y<0 | max_y>1)){
    pdf_y <- pdf_y/(max_y+scale_u-(min_y-scale_l))
    pdf_y_perc[1,] <- pdf_y_perc[1,]/(max_y+scale_u-(min_y-scale_l))
    pdf_y_perc[2,] <- pdf_y_perc[2,]/(max_y+scale_u-(min_y-scale_l))

    pdf_matrix <- pdf_matrix/(max_y+scale_u-(min_y-scale_l))

    return(list(pdf=pdf_y, cri=t(pdf_y_perc), mat=t(pdf_matrix), x=x))

  }else{
    return(list(pdf=pdf_y, cri=t(pdf_y_perc), mat=t(pdf_matrix), x=input))
  }


}

