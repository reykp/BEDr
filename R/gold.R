#' Bayesian Density Estimation Using \code{gold} (Gaussian process On a Log Density).
#'
#' Estimates a density from Bayesian techniques using a Gaussian process on a log density.
#' 
#' @usage 
#' gold(y, s1 = max(2.3 - 0.003 * n, 0.5), c1 = max([n/50], 1), s2 = NA, c2 = c1/2, N = 20000, graves = TRUE, scale_l = 0.00001, scale_u = 0.00001, poi = NA)
#'
#' @param y A numeric vector. \code{gold} estimates the density on this data.
#' @param s1 The standard deviation of the Gaussian prior. The default is set to a linear equaiton that decreases s1 as the sample size increases (see details).
#' @param c1 The correlation parameter of the Gaussian prior that controls the correlation in the covariance structure and smoothness in the density estimate. See details for the defaul value.
#' @param s2 The standard deviation of the Gaussian proposal distribution.
#' @param c2 The correlation parameter that controls the covariance structure of the Gaussian proposal distribution.
#' @param N The number of iterations to run in the algorithm. The default is 20,000.
#' @param graves An option to have the standard deviation of the proposal distribution (s1) chosen using an automatic step size selection (See Graves (2011)) (DEFAULT is FALSE).
#' @param scale_l A value >= 0 controlling the scaling based on the minimum data value. The default is 0.00001 (see details).
#' @param scale_u value >= 0 controlling the scaling based on the maximum data value. The default is 0.00001 (see details).
#' @param poi Points of interest at which to estimate the density. \code{gold_pdf} will produce estimates and credible intervals of the density at any point, but \code{mat} returned by \code{gold_pdf}  and \code{gold_cdf} only return the posterior draws if the density is estimated at that point.
#'
#' @export
#' @importFrom matrixcalc svd.inverse
#' @importFrom MASS mvrnorm
#'
#' @details
#'
#' The density being estimated takes the form below:
#'
#' \eqn{f(x) =}
#' \deqn{exp(g(x)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' where g(x) is an unknown log density.
#'
#' Given that g(x) is unknown, the integral in the normalizing constant is estimated using a weighted average, and the set of unknown parameters that receive a prior is g(x) at a finite set of points. This finite set of points includes each data point, as well as a grid to ensure good coverage.
#' 
#' This density operates on data between 0 and 1. Thus, if the input is not between 0 and 1, it is standardized. The formula for scaling is below:
#'
#'  \deqn{(y-(min(y)-scale_l))/(max(y)+scale_u-(min(y)-scale_l))}.
#'  
#'  Values of 0 for the scale parameters indicates the density will only be the edges of the data in \code{y}. The default for \code{scale_l} and \code{scale_u} is 0.00001 since some densities cannot be estimated at 0 and 1.
#'   
#' There is a default recommendation based on the size of \code{y}. The recommendation starts and 1 and then increases by 1 for each additional 50 data points:
#' 
#' Default: c1 = max(round(n/50), 1)
#' 
#' A default is also recommended for s1 in the form of: 1.1 - 0.002 * n. However, a minimum value of .05 is enforced, but the user can specify their own value for s1 to overide this.
#' @return
#'
#' \code{gold} returns a list containing the density estimate results.
#'
#' \item{\code{G}}{A matrix with the posterior draws for g(x) at a finite number of points. The number of rows is the number of iterations. The number of columns is the number of data points plus the grid points.}
#' \item{\code{widths}}{A vector with the widths around each x used for g(x). These are the widths used in the weighted average to estimate the integral in the normalizing constant.}
#' \item{\code{x}}{A vector with the points at which g(x) was estimated.}
#' \item{\code{y}}{A vector containing the data introduced to \code{gold} for density estimation.}
#' \item{\code{ar}}{The acceptance rate from the random walk proposals.}
#' \item{\code{prior}}{A vector of all prior and proposal parameter values.}
#' \item{\code{poi}}{The points of interest if \code{poi} is non-missing.}
#'
#' @references 
#' 
#' Graves, T. (2011). \emph{Automatic Step Size Selection in Random Walk Metropolis Algorithms}.
#' 
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'gold' on data sampled from a Beta(5, 1) distribution with 200 data points.
#' y <- rbeta(200, 5, 1)
#' # Use all defaults: this also uses Graves method to tune the standard deviation for the proposal
#' # distribution.
#' gold_beta <- gold(y = y)
#'
#' # Check trace plots
#' gold_mcmcplots(gold_beta)
#'
#' # Examine the estimate of the PDF
#' gold_plot(gold_beta, data = TRUE)
#'
#' # Examine the estimate of the CDF with the empirical CDF overlaied
#' gold_plot(gold_beta, type = "cdf", data = TRUE)
#' 
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'gold' on data sampled from a Normal(0, 1) distribution with 100 data points.
#' y <- rnorm(100, 0, 1)
#' # Specify all prior parameters: tune until achieve desired acceptance rate
#' gold_norm <- gold(y = y, s1 = 2, c1 = 1, s2 = 2, c2 = 0.5, N = 20000, graves = FALSE, scale_l = 2*sd(y), scale_u = 2*sd(y))
#'
#' # Check what s2 was chosen
#' gold_norm$prior
#' 
#' # Check acf plots
#' gold_mcmcplots(gold_norm, type = "acf")
#'
#' # Examine the estimate of the PDF with credible intervals
#' gold_plot(gold_norm, cri = TRUE)
#'
#' # Examine the estimate of the CDF using an interactive graph
#' gold_plot(gold_norm, type = "cdf", interact = TRUE)
#' 
#' ## --------------------------------------------------------------------------------
#' ## Bimodal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # Sample 400 random uniforms
#' u  <- runif(400)
#' y <- rep(NA,400)
#' # Sampling from the mixture
#' for(i in 1:400){
#'   if(u[i] < 0.3){
#'    y[i]  <-  rnorm(1, 0, 1)
#'   }else {
#'    y[i] <- rnorm(1, 4, 1)
#'   }
#' }
#' 
#' # First run 'gold' on data sampled from a bimodal distribution with 200 data points.
#' # Specify several points of interest
#' # USe graves, but use less correlation then default by making 'c1' large.
#' gold_bimodal <- gold(y = y, s1 = 2, c1 = 15, c2 = 7, graves = TRUE, poi = c(0, 6.5), scale_l = .5, scale_u = .5)
#'
#' # Check the running mean plots
#' gold_mcmcplots(gold_bimodal, type = "rm")
#' 
#' # Plot the PDF with the data and credible intervals
#' gold_plot(gold_bimodal, data = TRUE, cri = TRUE)
#'
#' # Plot the CDF
#' gold_plot(gold_bimodal, type = "cdf")
#' 
#' # Plot histograms at estimates of the CDF at the 'poi'
#' bimodal_cdf <- gold_cdf(x = c(0, 6.5), gold_bimodal)
#' bimodal_cdf$cdf
#' # x = 0
#' hist(bimodal_cdf$mat[,1])
#' # x = 6.5
#' hist(bimodal_cdf$mat[,2])

gold <- function(y, s1 = NA, c1 = NA, s2 = NA, c2 = NA, N = 20000, graves=TRUE, scale_l = .00001, scale_u = .00001, poi = NA){
  
  # Specify default for s1 if missing
  if(is.na(s1)){
    s1 = max(2.3 - 0.003 * length(y), 0.5)
  }
  
  # Specify default for c1 if missing
  if(is.na(c1)){
    c1 <- max(round(length(y)/50),1)
  }
  
  if(is.na(c2)){
    c2 <- c1/2
  }
  
  # All error checking ********************************************************************
  
  # No missing values in y
  if(length(which(is.na(y))) > 0){
    stop(" 'y' cannot contain missing values.")
  }
  
  # Error check for appropriate range of parameters
  if(s1<=0){
    stop(" 's1' should be > 0.")
  }

  if(c1 < 0){
    stop(" 'c1' should be >= 0.")
  }
  
  
  if(c2 < 0){
    stop(" 'c2' should be >= 0.")
  }
  
  if (scale_l < 0.0) {
    stop(" 'scale_l' should be >= 0.")
  }
  if (scale_u < 0.0) {         	# scale paramter needs to be greater than or equal to 0
    stop(" 'scale_u' should be >= 0.")
  }
  
  # Make sure numbers in 'poi' are unique, otherwise throws algorithm off  
  if(!is.na(poi[1])){
    total_length <- length(y)+length(poi)
    if(length(unique(poi))!=length(poi)){
      stop("The values in 'poi' are not unique.")
    }
    
    # Check for missing values
    if(length(which(is.na(poi))) > 0){
      stop(" 'poi' cannot contain missing values.")
    }
  }else{
    total_length <- length(y)
  }
  
  # Print a warning if choose graves, but also specify a value for 's2'
  if(graves == TRUE){
  if(!is.na(s2)){
    warning(" 'graves' == TRUE. Specified value for 's2' will be ignored.")
  }
  }
    
  # Make sure N is an integer 
  if(N/ceiling(N)<1){
    stop(" 'N' should be an integer value.")
  }
  
  # Make sure N is greater than 0
  if(N<=0){
    stop(" 'N' should be an integer greater than zero.")
  }
  
  
  #### end of error checking
  
  # Create vector to contain original data
  y_orig <- NA
  for(j in 1:length(y)){
    y_orig[j]<-y[j]
  }

  # Sort the data
  y_orig <- sort(y_orig)
  
  # Find the maximum and minimum
  # If points of interest are specified outside the range of the data, 
  # This needs to be incorporated into the max and min otherwise the standardized
  # valeus for 'poi' might be greater than 1 or less than 0
  if(!is.na(poi[1])){
    max_y <- max(y, poi)
    min_y <- min(y, poi)
  }else{
    max_y <- max(y)
    min_y <- min(y)
  }

  # Mark positions that are data
  names(y_orig) <- rep("data", length(y_orig))
  
  # Add points of interest to the data if there are any entered
  if(is.na(poi[1])){
      y_all_orig <- y_orig
  } else{
    names(poi) <- rep("poi", length(poi))
    y_all_orig <- sort(c(y_orig, poi))
  }
  
  # y_all becomes the scaled data between 0 and 1 to use in the algorithm
  y_all <- y_all_orig
  
  # Check if data is outside (0,1) range and scale it if it is
  if((max_y>1)|(min_y<0)){
    for(j in 1:length(y_all_orig)){
      y_all[j] <- (y_all_orig[j]-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l))
    }
  }

  # Sort the data again
  y_all <- sort(y_all)
  
  # Create grid:0, 0.01,..., 0.99, 1
  x_full <- 0:100
  x_full <- x_full/100

  # Will eventually remove end points, but just need for calculating bin widths
  # Mark positions that are the grid and the end points
  names(x_full) <- rep("grid", length(x_full))
  names(x_full)[1] <- "end"
  names(x_full)[length(x_full)] <- "end"
  
  # x_full contains the end points, x removes them
  x <- x_full[-1]
  x <- x[-length(x)]

  # Get summary of data to count how many times each data point appeared
  y_d <- data.frame(y_all_orig,y_all, names(y_all_orig))
  names(y_d) <- c("y_all_orig", "y_all", "names.y_all_orig.")
  y_d <- y_d %>% group_by(y_all_orig,y_all, names.y_all_orig.) %>% summarise(n=n()) 

  # Data to add to grid
  y_add <- y_d$y_all
  # Maintain names of data points: "data" and "poi"
  names(y_add) <- y_d$names.y_all_orig.

  # if contains "poi", remove rows since this data set should contain only the 'y' input
  if(length(which(y_d$names.y_all_orig.=="poi"))!=0){
    y_d <- y_d %>% filter(names.y_all_orig. != "poi")
  }
  
  # Add data to grid
  x <- c(x,y_add)
  # Sort the data
  x <- sort(x)
  # Sort the full data
  x_full <- sort(c(x_full,y_add))

  # calculate the distances between points
  pairwise_dist <- NA
  for(i in 1:(length(x_full)-1)){
    pairwise_dist[i] <- abs(x_full[{i+1}]-x_full[i])
  }

  # The grid is designed to have distances of .01 between them, if there are two data points within .01 of each other
  # remove grid point in between
  drop <- NA
  r_index <- 1
  check <- which(!(names(x)%in% c("data", "poi")))

  for(j in check){
    if(j>1&j<length(x)){
      if(names(x)[{j-1}]=="data"&names(x)[{j+1}]=="data"){
        distance <- x[{j+1}]-x[{j-1}]
        if(distance<.01&names(x)[j]=="grid"){
          drop[r_index] <- j
          r_index <- r_index+1
        }
      }
    }
  }
  
  # Just in case none to drop
  if(!is.na(drop[1])){
  x <- x[-drop]
  }
  x_full <- c(0,x,1)

  # Calculate the midpoints of the remaining points
  midpoints <- NA
  for (i in 1:(length(x)+1)){
    midpoints[i] <- (x_full[{i+1}]-x_full[{i}])/2+x_full[{i}]
  }

  # Calculate the widths between each pair of points
  widths <- NA
  for (i in 1:(length(x))){
    widths[i] <- abs(midpoints[{i+1}]-midpoints[{i}])
  }

  # Covariance matrix for prior distribution
  cov_prior <- matrix(0, nrow=length(x), ncol=length(x))
  for (i in 1:length(x)){
    for(j in 1:length(x)){
      d <- abs(x[i]-x[j])
      if(c1*d<.5){
        cov_prior[i,j] <- s1^2*(1-6*(c1*d)^2+6*(c1*d)^3)
      }else if (c1*d <1){
        cov_prior[i,j] <- s1^2*(2*(1-c1*d)^3)

      } else {
        cov_prior[i,j] <- 0

      }
    }
  }

  # Use SVD
  cov_prior_inv <- matrixcalc::svd.inverse(cov_prior)
  G <- matrix(0, nrow=N, ncol=length(x))
  
  # Starting values
  G[1,] <-rep(0, ncol(G))


  ###GRAVES################

  if(graves==TRUE){
  # Vector of values to establish trend in acceptance rate
  s2_graves <- c(.001,.01,.03,.06,.08,.1,.2,.3,.6,.8,1,3,5)
  # Number of step sizes to try
  n_step_size <- length(s2_graves)
  # The target is .35, slightly lower than .4 due to number of parameteres
  target <- .35
  # Store acceptance rates
  ar <- NA
  
  #Which which data points in x are the data
  which_x <- which(names(x)=="data")

  # Loop through each step size suggestion
  for (i in 1:n_step_size){
    
    ar_G <- NA
    # Covariate matrix
    cov_x <- matrix(0, nrow=length(x), ncol=length(x))
    for (l in 1:length(x)){
      for(j in 1:length(x)){
        cov_x[l,j] <- s2_graves[i]*exp(-c2*(x[l]-x[j])^2)
      }
    }

    for (l in 1:ncol(cov_x)){
      cov_x[l,l] <- cov_x[l,l]+.0000000000001
    }
    
    # Set up for multivariate simultions from normal
    n <- 1
    mu <- rep(0, length(x))
    tol = 1e-06
    p <- length(mu)
    eS <- svd(cov_x)
    ev <- eS$d
    if (!all(ev >= -tol * abs(ev[1L])))
      stop("'Sigma' is not positive definite")

    mvrnorm_calc <- eS$u %*% diag(sqrt(pmax(ev, 0)), p)

    # Store iterations for tuning purposes 
    tune<- matrix(0, nrow=201, ncol=ncol(cov_x))
    tune[1,] <- G[1,]

    for(j in 2:101){
      
      # Simulate length(x) values from standard normal
      rw <- matrix(rnorm(p * n), n)

      # Get proposal
      rw <- drop(mu) + mvrnorm_calc %*% t(rw)
      g_prop <- tune[{j-1},]+rw
      
      # Numerator and denominator for acceptance ratio on the log scale
      ap_num_log <- -length(y)*log(sum(widths*exp(g_prop)))+sum(g_prop[which_x])-.5*t(g_prop)%*%cov_prior_inv%*%g_prop
      ap_den_log <- -length(y)*log(sum(widths*exp(tune[{j-1},])))+sum(tune[{j-1},which_x])-.5*t(tune[{j-1},])%*%cov_prior_inv%*%tune[{j-1},]

      # Calculate acceptance ratio on the log scale
      ap_log <- ap_num_log-ap_den_log

      if (log(runif(1))<ap_log){
        tune[j,] <- g_prop
        ar_G[j] <- 1
      } else {
        tune[j,] <- tune[j-1,]
        ar_G[j] <- 0
      }
    }
    ar_G <- ar_G[-1]
    
    # Calculate the number that accepted
    ar[i] <- sum(ar_G)
  }

  # Relationship between log acceptance ratio and log of s2 is linear
  test <- glm(cbind(ar,100-ar)~log(s2_graves), family="binomial")
  # Back calculate what s2 leads to .35 acceptance rate
  s2_start <- exp((log(target/(1-target))-test$coefficients[1])/test$coefficients[2])

  names(s2_start) <- NULL
  # s2 should not be too huge or small, if it is, algorithm is not converging so 
  # pick more reasonable value by finding closest in suggestions to 
  if(s2_start>15 | s2_start<.0001){
    s2_start <- s2_graves[which.min(abs(ar - target*100)) ]
  }
  ##########
  s2 <- s2_start
}

# Calculate covariance matrix for proposals
cov_x <- matrix(0, nrow=length(x), ncol=length(x))
for (i in 1:length(x)){
  for(j in 1:length(x)){
    cov_x[i,j] <- s2*exp(-c2*(x[i]-x[j])^2)
  }
}

# Ensures has inverse
for (i in 1:ncol(cov_x)){
  cov_x[i,i] <- cov_x[i,i]+.0000000000001
}

####SVD for multivariate proposals#####################
n <- 1
mu <- rep(0, length(x))
tol = 1e-06
p <- length(mu)
eS <- svd(cov_x)
ev <- eS$d
if (!all(ev >= -tol * abs(ev[1L])))
  stop("'Sigma' is not positive definite")

mvrnorm_calc <- eS$u %*% diag(sqrt(pmax(ev, 0)), p)

# Collect acceptance rate
ar_G <- NA

# Get indices for which parameters are actual data verses part of the grid
which_x <- which(names(x)=="data")

# Metropolis-Hastings algorithm
for (i in 2:N){

    # Simulate length(x) values from standard normal
    rw <- matrix(rnorm(p * n), n)
    rw <- drop(mu) + mvrnorm_calc %*% t(rw)

    # Add to previous iteration
    g_prop <- G[{i-1},]+rw

    # Calculate log of numerator and denominator of acceptance probability
    ap_num_log <- -length(y)*log(sum(widths*exp(g_prop)))+sum(y_d$n*g_prop[which_x])-.5*t(g_prop)%*%cov_prior_inv%*%g_prop
    ap_den_log <- -length(y)*log(sum(widths*exp(G[{i-1},])))+sum(y_d$n*G[{i-1},which_x])-.5*t(G[{i-1},])%*%cov_prior_inv%*%G[{i-1},]

    # Calculate log of acceptance probability
    ap_log <- ap_num_log-ap_den_log

    if (log(runif(1))<ap_log){
      G[i,] <- g_prop #Acceptance proposal
      ar_G[i] <- 1
    } else {
      G[i,] <- G[i-1,] #Reject proposal
      ar_G[i] <- 0
    }

    if(i %% 1000 == 0){
     print(paste("Iteration:", i))
    }
    
  }


  # Acceptance Rates
  ar_G <- ar_G[-1]
  print("Acceptance Rate:")
  print(sum(ar_G)/length(ar_G))
  if(is.na(poi[1])){
    return(list(G=G,widths=widths,x=x,y=y_orig,ar=sum(ar_G)/length(ar_G), prior=c(s1 = s1, c1 = c1, s2 = s2, c2 = c2), scale_l = scale_l, scale_u = scale_u))
  }else{
      return(list(G=G,widths=widths,x=x,y=y_orig,ar=sum(ar_G)/length(ar_G), prior=c(s1 = s1, c1 = c1, s2 = s2, c2 = c2), scale_l = scale_l, scale_u = scale_u, poi = poi))
  }
}
