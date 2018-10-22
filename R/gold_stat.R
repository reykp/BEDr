#' Statistics from GOLD.
#'
#' Calculates a variety of statistics from the \code{gold} density estimate.
#'
#' @param stat A value indicating choice of statistic (see details).
#' @param gold_output The list returned by \code{gold} containing the density estimate results.
#' @param p A list of quantiles if quantiles are requested (see details).
#' @param burnin The desired burnin to discard from the results. If no values is entered, the default is half the number of iterations.
#'
#' @export
#'
#' @details
#'
#' The form of the density is below:
#'
#' \deqn{f(x) =}
#' \deqn{exp(g(x)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' where g(x) is an unknown log density.
#'
#' Given that g(x) is unknown, the normalizing constant is estimated using a weighted average and the set of unknown paramters that recieve a prior is g(x) at a finite set of points.
#'
#' \eqn{E[X] = \int_{0}^{1} x * exp(g(x)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' \eqn{Var[X] = \int_{0}^{1} (x-x^2) * exp(g(x)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' The inverse CDF is specified below:
#'
#' \strong{Options for} \code{stat}
#'
#' Several of the standard statistics are available through the function \code{gold_stat}.
#' \itemize{
#'     \item \code{"mean" or "m"}: The E[x] is calculated as in the deatils section for each iteration, and the average of this result is returned as the mean.
#'     \item \code{"var" or "v"}: The Var[x] is calculated as in the deatils section for each iteration, and the average of this result is returned as the variance.
#'     \item \code{"quant" or "q"}: Returns the quantiles specified in 'p' using the inverse CDF described in the details section.
#'     }
#'
#' \strong{Options for} \code{p}
#'
#' Default is NA is quantiles are not the desired statistics. If 'quant' is specified for \code{stat}, p is a single value or vector of values between 0 and 1.
#'
#' @return
#'
#' \code{gold_stat} the statistic of choice and the credible intervals associated with it.
#'
#'If 'mean' or 'm' is specified for \code{stat}:
#' \item{\code{mean}}{The E[X] from the details is calulated at each iteration and the mean of these is returned.}
#' \item{\code{cri}}{The credible intervals for E[X] calculated using the 0.025th and 0.975th quantiles.}
#'
#'If 'var' or 'v' is specified for \code{stat}:
#' \item{\code{variance}}{The Var[X] from the details is calulated at each iteration and the mean of these is returned.}
#' \item{\code{cri}}{The credible intervals for Var[X] calculated using the 0.025th and 0.975th quantiles.}
#'
#'If 'quant' or 'q' is specified for \code{stat}:
#' \item{\code{quantiles}}{The quantils specified in 'p' is calulated using the inverse CDF in the details at each iteration and the mean of these is returned.}
#' \item{\code{cri}}{The credible intervals for the quantiles calculated using the 0.025th and 0.975th quantiles.}
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'gold' on data sampled from a Beta(2, 6) distribution with 100 data points.
#' y <- rbeta(100, 2, 6)
#' gold_beta <- gold(y, s1 = 1, c1 = 1, s2 = 1, c2 = 0.8, MH_N = 20000)
#'
#' #Get an estimate of the mean and its credible intervals
#' gold_stat(gold_beta, stat="mean")
#'
#' #Get an estimate of the variance and its credible intervals
#' gold_stat(gold_beta, stat="var")
#'
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'gold' on data sampled from a Normal(0, 1) distribution with 200 data points.
#' y <- rnorm(200, 0, 1)
#' gold_norm <- gold(y, s1 = 1, c1 = 1, s2 = 0.5, c2 = 0.8, MH_N = 20000)
#'
#' #Get an estimate of the quantils and their credible intervals
#' gold_stat(gold_norm, stat="q", p=c(0.1, 0.5, 0.9))



gold_stat <- function(gold_output, stat = "m", p=NA,burnin=NA, print = TRUE){



  G <- gold_output$G
  widths <- gold_output$widths

  if(is.na(burnin)){
    burnin <- nrow(G)/2
  }

  if(burnin>nrow(G)){
    stop("The specified burnin is greater than the number of iterations.")
  }


  if(!is.na(p[1])&!(stat%in%c("q", "quantile", "quant"))){
    stop("p specified, but quantile not requested as the statistic.")
  }

  # Get the original data
  y_orig <- gold_output$y
  # Find max and minimums
  #Find min and max of data
  if(!is.null(gold_output[["poi"]])){
    max_y <- max(y_orig, gold_output$poi)
    min_y <- min(y_orig, gold_output$poi)
  }else{
    min_y <- min(y_orig)
    max_y <- max(y_orig)
  }

  # Get points at which to 
  x_gold <- gold_output$x

  G_burnin <- G[{burnin+1}:nrow(G),]

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

  scale_l <- gold_output$scale_l
  scale_u <- gold_output$scale_u
  
  if((min_y<0) | (max_y>1)){
    mean_matrix <- (G_PDF %*% x_gold)*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l)
    mean_matrix_scaled <- (G_PDF %*% x_gold)
    var_matrix <- ((G_PDF %*% (x_gold^2))-mean_matrix_scaled^2)*(max_y+scale_u-(min_y-scale_l))^2
  }else{
    mean_matrix <- (G_PDF %*% x_gold)
    var_matrix <- ((G_PDF %*% (x_gold^2))-mean_matrix^2)
  }

# 
#   if((min_y<0) | (max_y>1)){
#     input <- (x_gold)*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l);
#     duos_CDF <- gold_cdf(input,gold_output,burnin)
#   }else{
#     input <- x_gold
#     duos_CDF <- gold_cdf(x_gold,gold_output,burnin)
#   }
# 
#   x_gold <- duos_CDF$x
#   cdf_y <- duos_CDF$cdf
# 
#   CDF_mat <- duos_CDF$mat
#   cdf_y_perc <- apply(CDF_mat, 2, quantile, probs=c(.025, .975))


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
  CDF_mat <- t(apply(G_CDF_start, 1,cumsum))
  
  #Calculate posterior mean cdf
  cdf_y <- NA
  for (j in 1:ncol(CDF_mat)){
    cdf_y[j] <- mean(CDF_mat[,j])
  }
  
  
  #Calculate percentiles
  cdf_y_perc <- apply(CDF_mat, 2, quantile, probs=c(.025, .975))
  
  
  #Create empty list to contain quantiles

  if(!is.na(p[1])){
    quantiles <- NA
    quantiles_cri <- matrix(NA, nrow=length(p), ncol=2)

    for(i in 1:length(p)){
      if((max_y>1)|(min_y<0)){
        ss_indiv <- smooth.spline(x_gold, cdf_y)
        finv1 <- splinefun(ss_indiv$y, ss_indiv$x)
        quantiles[i] <- finv1(p[i])*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l)
        
        
        ss_indiv_q1 <- smooth.spline(x_gold, t(cdf_y_perc)[,1])
        finv2 <- splinefun(ss_indiv_q1$y, ss_indiv_q1$x)
        quantiles_cri[i,2] <- finv2(p[i])*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l)
        
        ss_indiv_q2 <- smooth.spline(x_gold, t(cdf_y_perc)[,2])
        finv3 <- splinefun(ss_indiv_q2$y, ss_indiv_q2$x)
        quantiles_cri[i,1] <- finv3(p[i])*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l)
       
      }else{
        ss_indiv <- smooth.spline(x_gold, cdf_y)
        finv1 <- splinefun(ss_indiv$y, ss_indiv$x)
        quantiles[i] <- finv1(p[i])
    
      
        ss_indiv_q1 <- smooth.spline(x_gold, t(cdf_y_perc)[,1])
        finv2 <- splinefun(ss_indiv_q1$y, ss_indiv_q1$x)
        quantiles_cri[i,1] <- finv2(p[i])
     
        ss_indiv_q2 <- smooth.spline(x_gold, t(cdf_y_perc)[,2])
        finv3 <- splinefun(ss_indiv_q2$y, ss_indiv_q2$x)
        quantiles_cri[i,2] <- finv3(p[i])
      } 
    }
  }


  if(stat%in%c("mean", "m")){
    mean_list <- list(mean=mean(mean_matrix), cri=quantile(mean_matrix,c(0.025, 0.975)), mat = mean_matrix)
    if(print == TRUE){
      print(mean_list[1])
      print(mean_list[2])
    }
    return(mean_list)
  }else if (stat%in%c("var", "v")){
    var_list <- list(variance=mean(var_matrix), cri=quantile(var_matrix,c(0.025, 0.975)), mat = var_matrix)
    if(print == TRUE){
      print(var_list[1])
      print(var_list[2])
    }
    return(var_list)
  }else{
    quant_list <- list(quantiles=quantiles, cri=quantiles_cri)
    if(print == TRUE){
      print(quant_list[1])
      print(quant_list[2])
    }
    return(quant_list)
  }

}

