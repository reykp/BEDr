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
#' @importFrom Hmisc Lag



duos_stat <- function(stat,duos_output, p=NA,burnin=NA,scale=FALSE){



  #Cutpoints
  C <- duos_output[[1]]

  if(is.na(burnin)){
    burnin <- nrow(C)/2
  }

  if(burnin>nrow(C)){
    stop("The specified burnin is greater than the number of iterations.")
  }
  #Bin probabilities
  P <- duos_output[[2]]

  #Number of cutpoints
  k <<- ncol(C)

  if(!is.na(p)&!(stat%in%c("q", "quantile", "quant"))){
  stop("p specified, but quantile not requested as the statistic.")
  }


  y_orig <- duos_output[[3]]
  min_y <<- min(y_orig)
  max_y <<- max(y_orig)

  #Calculate mean and variance
  stats_apply <- function(x){
    c_full <- c(0, x[1:k],1)
    p <- x[{k+1}:length(x)]
    c_full_lag <- Lag(c_full)[-1]
    c_full <- c_full[-1]

    if((min_y<0) | (max_y>1)){
      #mean
      m_scaled <- (sum((c_full+c_full_lag)*p)/2)
      m <- (sum((c_full+c_full_lag)*p)/2)*(max_y+.00001-(min_y-.00001))+(min_y-.00001)

      c_full_lag3 <- c_full_lag^3
      c_full3 <- c_full^3

      v <- ((sum((c_full3-c_full_lag3)/3*(p/(c_full-c_full_lag))))-m_scaled^2)*(max_y+.00001-(min_y-.00001))^2
    }else{
      #mean
      m <- sum((c_full+c_full_lag)*p)/2

      c_full_lag3 <- c_full_lag^3
      c_full3 <- c_full^3

      v <- sum((c_full3-c_full_lag3)/3*(p/(c_full-c_full_lag)))-m^2

    }

    return(c(m,v))
  }

  #Calculates quantiels
  duos_quant <- function(x){

    c_full <- c(0, x[1:k],1)
    p <- x[{k+1}:length(x)]
    p_bounds <- c(0, cumsum(p))

    for(i in 1:{k+1}){
      if((q_loop>=p_bounds[i])&(q_loop<p_bounds[{i+1}])){
        if((min_y<0) | (max_y>1)){
        return((((q_loop-p_bounds[i])*(c_full[{i+1}]-c_full[i])/p[i])+c_full[i])*(max_y+.00001-(min_y-.00001))+(min_y-.00001))
        }else{
          return((q_loop-p_bounds[i])*(c_full[{i+1}]-c_full[i])/p[i]+c_full[i])
        }
      }
    }
  }

    #Take burnin into account
  C_sub <- C[burnin:nrow(C),]
  P_sub <- P[burnin:nrow(C),]

  #Get matrix of pdf values
  stats_matrix <- apply(cbind(C_sub,P_sub), 1, stats_apply)
  mean(stats_matrix[1,])
  mean(stats_matrix[2,])

  #Create empty list to contain quantiles
  quantiles <- list()


  for(i in 1:length(q)){
    q_loop <<- q[i]
    quant_matrix <- apply(cbind(C_sub,P_sub), 1, duos_quant)
    quantiles[[i]] <- mean(quant_matrix)
  }

  if(stat=="mean"){
    return(mean(stats_matrix[1,]))
  }else if (stat=="var"){
    return(mean(stats_matrix[2,]))
  }else{
    return(quantiles)
  }

}
