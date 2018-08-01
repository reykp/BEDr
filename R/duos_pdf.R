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
#' @importFrom Rcpp sourceCpp



duos_pdf <- function(x, duos_output, burnin,scale=FALSE){



  #Cutpoints
  C <- duos_output[[1]]
  if(burnin>nrow(C)){
    stop("The specified burnin is greater than the number of iterations.")
  }
  #Bin probabilities
  P <- duos_output[[2]]
  #Number of cutpoints
  k <<- ncol(C)
  #Calculate pdf at:
  input <<- x

  y_orig <- duos_output[[3]]
  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if((min_y<0 | max_y>1)){
    input <<- (x-(min_y-.00001))/(max_y+.00001-(min_y-.00001))
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
  pdf_matrix <- apply(cbind(C_sub,P_sub), 1, pdf_forapply)
  pdf_y <- apply(pdf_matrix, 1, mean)
  pdf_y_perc <- apply(pdf_matrix, 1, quantile, probs=c(.025, .975))


  if(scale==FALSE & (min_y<0 | max_y>1)){
    pdf_y <- pdf_y/(max_y+.00001-(min_y-.00001))
    pdf_y_perc[1,] <- pdf_y_perc[1,]/(max_y+.00001-(min_y-.00001))
    pdf_y_perc[2,] <- pdf_y_perc[2,]/(max_y+.00001-(min_y-.00001))

    return(list(pdf_matrix=pdf_matrix, pdf_y=pdf_y, pdf_percentiles=t(pdf_y_perc), x=x))

  }else{
    return(list(pdf_matrix=pdf_matrix, pdf_y=pdf_y, pdf_percentiles=t(pdf_y_perc), x=input))
  }


}
