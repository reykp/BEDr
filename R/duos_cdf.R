#' Estimate of CDF from DUOS.
#'
#' Calculates a posterior estimate of the CDF based on the output from \code{duos} for an individual
#' or vector of values.
#'
#' @param x A single value or vector of values at which to calculate the CDF. These values are to be entered on the scale of the data (i.e. values can fall outside of 0 and 1).
#' @param duos_output The list returned by \code{duos} contain the density estimate results.
#' @param burnin The desired burnin to discard from the results. If left as NA, the default is half the number of iterations.
#' @param scale This value TRUE/FALSE indicates whether to return scaled or unscalled results IF the original data does not fall between 0 and 1. The default is FALSE (i.e. returns results on the original data scale).
#'
#' @export
#'
#' @details
#'
#' The function \code{duos_cdf} returns the posterior mean CDF. The CDF is calculated based on the following equation at each iteration:
#'
#' \deqn{F(x) =}
#'
#' \deqn{(\pi_1) / (\gamma_1) * x , 0 \le x < \gamma_1}
#' \deqn{\pi_1 + (\pi_2) / (\gamma_2-\gamma_1) * (x-\gamma_1) , \gamma_1 \le x  < \gamma_2}
#' \deqn{\pi_1 + \pi_2 + (\pi_3) / (\gamma_3-\gamma_2) * (x-\gamma_2) , \gamma_2 \le x  < \gamma_3}
#' \deqn{...  , ... \le x  < ...}
#' \deqn{\sum_{i=1}^{k-1} \pi_i + (\pi_k) / (\gamma_k-\gamma_{k-1}) * (x-\gamma_{k-1}) , \gamma_{k-1} \le x  < \gamma_k}
#' \deqn{\sum_{i=1}^{k}\pi_i + (\pi_{k+1}) / (1-\gamma_k) * (x-\gamma_k) , \gamma_k \le x  < 1}
#'
#' @return \code{duos_cdf} returns a list of the CDF results from DUOS.
#'
#' \item{\code{cdf}}{A vector of the posterior mean CDF values at each value in \code{x}.}
#' \item{\code{cri}}{A matrix with 2 columns and rows equaling the length of \code{x} containing the 95\% credible interval for the CDF at each of the points in \code{x}.}
#' \item{\code{mat}}{A matrix containing the CDF values for each \code{x} at EACH itertation after the burnin is discarded. The number of columns is the length of \code{x}.}
#' \item{\code{x}}{A vector containing the values at which to estimate the CDF. If the data is not between 0 and 1 and scale=TRUE, the scaled version of \code{x} is returned.}

#' @examples

#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(2,5) distribution wiht 100 data points.
#' duos_beta <- duos(rbeta(100, 2, 5), k=5, MH_N=20000)
#' cdf_beta <- duos_cdf(x = c(.01, .25, .6, .9), duos_beta)
#'
#' #Examine the CDF at 'x'
#' cdf_beta$cdf
#'
#' #Examine the credibal intervals of the CDF at 'x'
#' cdf_beta$cri
#'
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Normal(0,1) distribution with 50 data points.
#' duos_norm <- duos(rnorm(50, 0, 1), k=4, MH_N=20000)
#' cdf_norm <- duos_cdf(x=c(-2, -1, 0, .8, 1.8), duos_norm)
#'
#' #Examine the CDF at 'x'
#' cdf_norm$cdf
#'
#' #Examine the credibal intervals of the CDF at 'x'
#' cdf_norm$cri
#'
#Histogram of distribution of the PDF density estimate at 0.8
#' hist(pdf_norm$mat[,4])

duos_cdf <- function(x, duos_output, burnin=NA, scale=FALSE){

  #if burnin is NA, use half as default
  if(is.na(burnin)){
    burnin <- nrow(duos_output$C)/2
  }

  #Cutpoints
  C <- duos_output$C
  if(burnin>nrow(C)){
    stop("The specified burnin is greater than the number of iterations.")
  }

  #Bin probabilities
  P <- duos_output$P
  #Number of cutpoints
  k <<- ncol(C)
  #Calculate cdf at:
  input <<- x

  y_orig <- duos_output$y
  min_y <- min(y_orig)
  max_y <- max(y_orig)
  if((min_y<0 | max_y>1)){
    input <<- (x-(min_y-.00001))/(max_y+.00001-(min_y-.00001))
  }

  cdf_forapply <- function(x){
    #Vector for probabilities
    pr <- rep(0,length(input))
    #Vector for cutpoints from rows with 0 and 1 added in
    c_full <- c(0, x[1:k],1)
    #Bin probabiliites
    p <- x[{k+1}:length(x)]
    #Loop to get CDF estimates
    for (j in 1:(k+1)){
      pr[which(input<c_full[(j+1)] & input>=c_full[j])]<- sum(p[1:{j-1}])*ifelse(j>1,1,0)+(p[j]/(c_full[{j+1}]-c_full[{j}]))*(input[which(input<c_full[{j+1}] & input>=c_full[j])]-c_full[{j}])
    }
    return(pr)
  }
  #Take burnin into account
  C_sub <- C[burnin:nrow(C),]
  P_sub <- P[burnin:nrow(C),]

  #Get matrix of cdf values
  if(length(x)==1){
  cdf_matrix <- matrix(apply(cbind(C_sub,P_sub), 1, cdf_forapply), nrow=1)
  }else{
    cdf_matrix <- apply(cbind(C_sub,P_sub), 1, cdf_forapply)
  }
  cdf_y <- apply(cdf_matrix, 1, mean)
  cdf_y_perc <- apply(cdf_matrix, 1, quantile, probs=c(.025, .975))


  if(scale==FALSE & (min_y<0 | max_y>1)){
    return(list(cdf=cdf_y, cri=t(cdf_y_perc), mat=cdf_matrix,x=x))


  }else{
    return(list(cdf=cdf_y, cri=t(cdf_y_perc), mat=cdf_matrix, x=input))
  }
}
