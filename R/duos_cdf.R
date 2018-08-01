#' Bayesina Density Estimate of CDF from DUOS.
#'
#' Calculates an estimate of the CDF based on the output from \code{duos} for an individual
#' or vector of values
#'
#' @param x A single value or vector of values to calculate the CDF at
#' @param duos_output The list returned by \code{duos}
#' @param burnin The desired burnin to discard from including in the estimate
#'
#' @export

duos_cdf <- function(x, duos_output, burnin=NA){

  #if burnin is NA, use half as default
  if(is.na(burnin)){
    burnin <- nrow(duos_output[[1]])/2
  }

  #Cutpoints
  C <- duos_output[[1]]
  if(burnin>nrow(C)){
    stop("The specified burnin is greater than the number of iterations.")
  }
  #Bin probabilities
  P <- duos_output[[2]]
  #Number of cutpoints
  k <<- ncol(C)
  #Calculate cdf at:
  input <<- x

  y_orig <- duos_output[[3]]
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
  cdf_matrix <- apply(cbind(C_sub,P_sub), 1, cdf_forapply)
  cdf_y <- apply(cdf_matrix, 1, mean)
  cdf_y_perc <- apply(cdf_matrix, 1, quantile, probs=c(.025, .975))


  if(scale==FALSE & (min_y<0 | max_y>1)){
    return(list(cdf_matrix=cdf_matrix, cdf_y=cdf_y, cdf_percentiles=t(cdf_y_perc), x=x))


  }else{
    return(list(pdf_matrix=pdf_matrix, pdf_y=pdf_y, pdf_percentiles=t(pdf_y_perc), x=input))
  }
}
