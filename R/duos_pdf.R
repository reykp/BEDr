#' Estimate of PDF from DUOS.
#'
#' Calculates a posterior estimate of the PDF based on the output from \code{duos} for an individual
#' or vector of values.
#'
#' @param x A single value or vector of values at which to calculate the PDF. These values are to be entered on the scale of the data (i.e. values can fall outside of 0 and 1).
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param burnin The desired burnin to discard from the results. If no values is entered, the default is half the number of iterations.
#' @param scale This value TRUE/FALSE indicates whether to return scaled or unscalled results IF the original data does not fall between 0 and 1. The default is FALSE (i.e. returns results on the original data scale).
#'
#' @export
#'
#' @details
#'
#' The function \code{duos_pdf} returns the posterior mean PDF. The PDF is calculated based on the following equation at each iteration:
#'
#' \deqn{f(x) =}
#'
#' \deqn{(\pi_1) / (\gamma_1) , 0 \le x < \gamma_1}
#' \deqn{(\pi_2) / (\gamma_2-\gamma_1) , \gamma_1 \le x  < \gamma_2}
#' \deqn{(\pi_3) / (\gamma_3-\gamma_2) , \gamma_2 \le x  < \gamma_3}
#' \deqn{...  , ... \le x  < ...}
#' \deqn{(\pi_k) / (\gamma_k-\gamma_{k-1}) , \gamma_{k-1} \le x  < \gamma_k}
#' \deqn{(\pi_{k+1}) / (1-\gamma_k) , \gamma_k \le x  < 1}
#'
#' @return \code{duos_pdf} returns a list of the PDF results from DUOS.
#'
#' \item{\code{pdf}}{A vector of the posterior mean PDF values at each value in \code{x}.}
#' \item{\code{cri}}{A matrix with 2 columns and rows equaling the length of \code{x} containing the 95\% credible interval for the PDF at each of the points in \code{x}.}
#' \item{\code{mat}}{A matrix containing the PDF values for each \code{x} at EACH itertation after the burnin is discarded. The number of columns is the length of \code{x}.}
#' \item{\code{x}}{A vector containing the values at which to estimate the PDF. If the data is not between 0 and 1 and scale=TRUE, the scaled version of \code{x} is returned.}

#' @examples

#' ## --------------------------------------------------------------------------------
#' ## Uniform Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Unif(0,1) distribution wiht 70 data points.
#' y <- runif(70)
#' duos_unif <- duos(y, k=4, MH_N=20000)
#' pdf_unif <- duos_pdf(x = c(.1, .5, .65, .98), duos_unif)
#'
#' #Examine the PDF at 'x'
#' pdf_unif$pdf
#'
#' #Examine the credibal intervals of the PDF at 'x'
#' pdf_unif$pri
#'
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Normal(2,2) distribution with 90 data points.
#' y <- rnorm(90, 2, 2)
#' duos_norm <- duos(y, k=5, MH_N=20000)
#' pdf_norm <- duos_pdf(x=c(-1.5, -.3, 0, .3, 2.1), duos_norm)
#'
#' #Examine the PDF at 'x'
#' pdf_norm$pdf
#'
#' #Examine the credibal intervals of the PDF at 'x'
#' pdf_norm$cri
#'
#' #Histogram of distribution of the PDF density estimate at -0.3
#' hist(pdf_norm$mat[,2])



duos_pdf <- function(x, duos_output, burnin=NA,scale=FALSE){



  #Cutpoints
  C <- duos_output$C

  if(is.na(burnin)){
    burnin <- nrow(C)/2
  }
  if(burnin>nrow(C)){
    stop("The specified burnin is greater than the number of iterations.")
  }
  #Bin probabilities
  P <- duos_output$P


  #Number of cutpoints
  k <<- ncol(C)
  #Calculate pdf at:
  input <<- x

  y_orig <- duos_output$y
  min_y <- min(y_orig)
  max_y <- max(y_orig)
  
  scale_l <- duos_output$scale_l
  scale_u <- duos_output$scale_u
  
  
  
  if(min_y<0|max_y>1){
    for (i in 1:length(x)){
      if(x[i]<(min_y-scale_l)|x[i]>(max_y+scale_u)){
        print(x[i])
        print("is out of the range of the data in 'y'.")
      }
    }    
  }
  
  if((min_y<0 | max_y>1)){
    input <<- (x-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l))
  }

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

  #Take burnin into account
  C_sub <- C[burnin:nrow(C),]
  P_sub <- P[burnin:nrow(C),]

  #Get matrix of pdf values
  #Get matrix of cdf values
  if(length(x)==1){
    pdf_matrix <- matrix(apply(cbind(C_sub,P_sub), 1, pdf_forapply), nrow=1)
  }else{
    pdf_matrix <- apply(cbind(C_sub,P_sub), 1, pdf_forapply)
  }
  pdf_y <- apply(pdf_matrix, 1, mean)
  pdf_y_perc <- apply(pdf_matrix, 1, quantile, probs=c(.025, .975))


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

