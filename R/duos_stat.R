#' Statistics for \code{duos}
#'
#' Calculates a variety of statistics from the \code{duos} density estimate.
#'
#' @usage 
#' duos_stat(duos_output, stat = "m", p = NA, burnin = NA, scale = FALSE)
#' 
#' @param duos_output The list returned by \code{duos} containing the density estimate results.
#' @param stat A value indicating choice of statistic (see details).
#' @param p A list of quantiles if quantiles are requested (see details).
#' @param burnin The desired burnin to discard from the results. If no value is entered, the default is half the number of iterations.
#' @param print By default, the posterior mean values and their credible intervals are printed. Change to FALSE if no printing is desired.
#'
#' @export
#' @useDynLib biRd
#' @importFrom Hmisc Lag
#'
#' @details
#'
#' The form of the density is below:
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
#' \eqn{E[X] = \sum_{i=1}^{k+1} ((\gamma_i+\gamma_{i-1}) * \pi_i) / 2}
#'
#' \eqn{Var[X] = \sum_{i=1}^{k+1} ((\gamma_i^3-\gamma_{i-1}^3) / 3) * \pi_i / (\gamma_i-\gamma_{i-1}) - (E[x])^2}
#'
#' The inverse CDF is specified below:
#'
#'  \eqn{F^{-1}(x) =}
#'  \deqn{x * (\gamma_1) / (\pi_1) , 0 \le x < \pi_1}
#'  \deqn{(x - \pi_1) * (\gamma_2 - \gamma_1) / (\pi_2) + \gamma_1 , \pi_1 \le x  < \pi_1 + \pi_2}
#'  \deqn{(x - \pi_1 + \pi_2) * (\gamma_3-\gamma_2) / (\pi_3) + \gamma_2 , \pi_1 + \pi_2 \le x  < \pi_1 + \pi_2 + \pi_3}
#'  \deqn{... ,  ... \le x  < ...}
#'  \deqn{(x - \sum_{i=1}^{k-1} \pi_i) * (\gamma_k - \gamma_{k-1}) / (\pi_k) + \gamma_{k-1} , \sum_{i=1}^{k-1} \pi_i \le x  < \sum_{i=1}^{k} \pi_i}
#'  \deqn{(x - \sum_{i=1}^{k} \pi_i) * (1 - \gamma_{k}) / (\pi_{k+1}) + \gamma_{k} ,  \sum_{i=1}^{k} \pi_i \le x  < 1}
#'
#' \strong{Options for} \code{stat}
#'
#' Several of the standard statistics are available through the function \code{duos_stat}.
#' \itemize{
#'     \item \code{"mean" or "m"}: The E[x] is calculated as in the details section for each iteration, and the average of this result is returned as the mean.
#'     \item \code{"var" or "v"}: The Var[x] is calculated as in the details section for each iteration, and the average of this result is returned as the variance.
#'     \item \code{"quant" or "q"}: Returns the quantiles specified in 'p' using the inverse CDF described in the details section.
#'     }
#'
#' \strong{Options for} \code{p}
#'
#' Default is NA if quantiles are not the desired statistics. If 'quant' is specified for \code{stat}, p is a single value or vector of values between 0 and 1.
#'
#' @return
#'
#' \code{duos_stat} the statistic of choice and the credible intervals associated with it.
#'
#'If 'mean' or 'm' is specified for \code{stat}:
#' \item{\code{mean}}{The E[X] from the details is calculated at each iteration and the mean of these is returned.}
#' \item{\code{cri}}{The credible intervals for E[X] calculated using the 0.025th and 0.975th quantiles.}
#'
#'If 'var' or 'v' is specified for \code{stat}:
#' \item{\code{variance}}{The Var[X] from the details is calculated at each iteration and the mean of these is returned.}
#' \item{\code{cri}}{The credible intervals for Var[X] calculated using the 0.025th and 0.975th quantiles.}
#'
#'If 'quant' or 'q' is specified for \code{stat}:
#' \item{\code{quantiles}}{The quantiles specified in 'p' is calculated using the inverse CDF in the details at each iteration and the mean of these is returned.}
#' \item{\code{cri}}{The credible intervals for the quantiles calculated using the 0.025th and 0.975th quantiles.}
#'
#' @return \code{duos_pdf} returns a list of the PDF results from \code{duos}.
#'
#' \item{\code{"mean", "var", or "quantiles"}}{Contains the mean, the variance, or a vector of the quantile values from \code{duos_stat}.}
#' \item{\code{cri}}{Contains the credible intervals for the requested statistics.}
#' \item{\code{mat}}{Contains a matrix of the posterior simulations for the requested statistic.}

#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(2, 6) distribution with 100 data points.
#' y <- rbeta(100, 2, 6)
#' duos_beta <- duos(y = y, k = 5, MH_N = 20000)
#'
#' # Get an estimate of the mean and its credible intervals
#' beta_mean <- duos_stat(duos_beta, stat = "mean")
#'
#' #Look at a histogram of the posterior simulations of the mean
#' hist(beta_mean$mat)
#' 
#' # Get an estimate of the variance and its credible intervals
#' beta_var <- duos_stat(duos_beta, stat = "var")
#'
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Normal(0, 1) distribution with 200 data points.
#' y <- rnorm(200, 0, 1)
#' duos_norm <- duos(y = y, k = 7, MH_N = 20000)
#'
#' # Get an estimate of the mean and credible intervals
#' mean_norm <- duos_stat(duos_norm, stat = "m")
#' 
#' #Get an estimate of the quantiles and their credible intervals
#' quant_norm <- duos_stat(duos_norm, stat = "q", p = c(0.1, 0.5, 0.9))

duos_stat <- function(duos_output,stat = "m", p=NA,burnin=NA,print=TRUE){
  
  # Get Cut-points
  C <- duos_output$C

  # Set defualt burnin
  if(is.na(burnin)){
    burnin <- ceiling(nrow(C)/2)
  }
  
  # Error checks 
  
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

  if(!(stat%in%c("mean","m", "var", "var", "quant", "q"))){
    stop("Please choose from the available statistics: 'mean', 'var', or 'quant'.")
  } 
 

  if(!is.na(p[1])&!(stat%in%c("q", "quantile", "quant"))){
    stop("p specified, but quantile not requested as the statistic.")
  }
  
  if(stat%in%c("q", "quantile", "quant")){
    if(max(p)>=1 | min(p)<=0){
      stop("Specified quantiles should be between 0 and 1.")
    }
  }
  
  # Get bin proportion parameters
  P <- duos_output$P

  # Number of cut-points
  k <<- ncol(C)

  # Get the original data
  y_orig <- duos_output$y
  # Get the max and min
  min_y <<- min(y_orig)
  max_y <<- max(y_orig)

  # Get the scale parameters
  scale_l <<- duos_output$scale_l
  scale_u <<- duos_output$scale_u
  
  #Calculate mean and variance
  stats_apply <- function(x){
    c_full <- c(0, x[1:k],1)
    p <- x[{k+1}:length(x)]
    c_full_lag <- Lag(c_full)[-1]
    c_full <- c_full[-1]

    if((min_y<0) | (max_y>1)){
      # Mean: needs to translated back to original scale
      m_scaled <- (sum((c_full+c_full_lag)*p)/2)
      m <- (sum((c_full+c_full_lag)*p)/2)*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l)

      c_full_lag3 <- c_full_lag^3
      c_full3 <- c_full^3

      v <- ((sum((c_full3-c_full_lag3)/3*(p/(c_full-c_full_lag))))-m_scaled^2)*(max_y+scale_u-(min_y-scale_l))^2
    }else{
      # mean
      m <- sum((c_full+c_full_lag)*p)/2

      c_full_lag3 <- c_full_lag^3
      c_full3 <- c_full^3

      v <- sum((c_full3-c_full_lag3)/3*(p/(c_full-c_full_lag)))-m^2

    }

    return(c(m,v))
  }

  # Calculates quantiles
  duos_quant <- function(x){

    c_full <- c(0, x[1:k],1)
    p <- x[{k+1}:length(x)]
    p_bounds <- c(0, cumsum(p))

    for(i in 1:{k+1}){
      if((q_loop>=p_bounds[i])&(q_loop<p_bounds[{i+1}])){
        if((min_y<0) | (max_y>1)){
        return((((q_loop-p_bounds[i])*(c_full[{i+1}]-c_full[i])/p[i])+c_full[i])*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l))
        }else{
          return((q_loop-p_bounds[i])*(c_full[{i+1}]-c_full[i])/p[i]+c_full[i])
        }
      }
    }
  }

  # Take burnin into account
  C_sub <- C[{burnin+1}:nrow(C),]
  P_sub <- P[{burnin+1}:nrow(C),]

  #Get matrix of pdf values
  stats_matrix <- apply(cbind(C_sub,P_sub), 1, stats_apply)

  # Create empty list to contain quantiles
  quantiles <- NA
  quantiles_cri <- matrix(NA, nrow=length(p), ncol=2)

  # Get quantile values
  if(!is.na(p[1])){
  for(i in 1:length(p)){
    q_loop <<- p[i]
    quant_matrix <- apply(cbind(C_sub,P_sub), 1, duos_quant)
    quantiles[i] <- mean(quant_matrix)
    quantiles_cri[i,] <- quantile(quant_matrix, c(0.025, 0.975))
  }
  }

  # Return results
  if(stat%in%c("mean", "m")){
    mean_list <- list(mean=mean(stats_matrix[1,]), cri=quantile(stats_matrix[1,],c(0.025, 0.975)), mat=stats_matrix[1,])
    if(print == TRUE){
    print(mean_list[1])
    print(mean_list[2])
    }
    return(mean_list)
  }else if (stat%in%c("var", "v")){
    var_list <- list(variance=mean(stats_matrix[2,]), cri=quantile(stats_matrix[2,],c(0.025, 0.975)),mat=stats_matrix[2,])
    if(print == TRUE){
    print(var_list[1])
    print(var_list[2])
    }
    return(var_list)
  }else{
    quant_list <- list(quantiles=quantiles, cri=quantiles_cri, mat = quant_matrix)
    if(print == TRUE){
    print(quant_list[1])
    print(quant_list[2])
    }
    return(quant_list)
  }

}
