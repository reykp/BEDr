#' Bayesina Density Estimation Using GOLD (Gaussian process On a Log Density).
#'
#' Bayesian density estimation through a Gaussian process on a log density
#'
#' @param y a numeric vector containing the data on which to estimate a density
#' @param s1 the standard deviation of the Gaussian prior
#' @param c1 the correlation parameter of the Gaussian prior that controls the correlation in the covariance structure
#' @param s2 the standard deviation of the Gaussian proposal distribution
#' @param c2 the correlation paramtet that controls the covariance structure of the Gaussian proposal distribution
#' @param MH_N the number of iterations
#'
#' @export
#' @importFrom matrixcalc svd.inverse
#' @importFrom MASS mvrnorm



gold <- function(y,s1, c1,s2, c2,MH_N){
  #Create vector to contain original data
  y_orig <- NA
  for(j in 1:length(y)){
    y_orig[j]<-y[j]
  }

  y_orig <- sort(y_orig)
  max_y <- max(y)
  min_y <- min(y)

  #Check if data is outside (0,1) range
  if((max_y>1)|(min_y<0)){
    for(j in 1:length(y)){
      y[j] <- (y_orig[j]-(min_y-.00001))/(max_y+.00001-(min_y-.00001))
    }
  }

  y <- sort(y)
  #Create grid:0,0,01,...,.99,1
  x_full <- 0:100
  x_full <- x_full/100

  #Will eventually remove end points, but just need to calculating bin widths
  #also mark positions that are the grid
  names(x_full) <- rep("grid", length(x_full))
  names(x_full)[1] <- "end"
  names(x_full)[length(x_full)] <- "end"

  #Mark positions that are data
  names(y) <- rep("data", length(y))

  #x_full contains the end points, x removes them
  x <- x_full[-1]
  x <- x[-length(x)]

  #Get summary of data to count how many times each data point appeared
  y_d <- data.frame(y_orig,y)
  y_d <- y_d %>% group_by(y_orig,y) %>% summarise(n=n())

  y_add <- y_d$y
  names(y_add) <- rep("data", length(y_add))

  #Add in the grid
  x <- c(x,y_add)
  #Sort the data
  x <- sort(x)
  #Sort the full data
  x_full <- sort(c(x_full,y_add))

  #calculate the distances between points
  pairwise_dist <- NA
  for(i in 1:(length(x_full)-1)){
    pairwise_dist[i] <- abs(x_full[{i+1}]-x_full[i])
  }

  #The grid is designed to have distances of .01 between them, if there are two data points within .01 of each other
  #remove grid point in between
  drop <- NA
  r_index <- 1
  check <- which(names(x)!="data")

  for(j in check){
    if(j>1){
      if(names(x)[{j-1}]=="data"&names(x)[{j+1}]=="data"){
        distance <- x[{j+1}]-x[{j-1}]
        if(distance<.01){
          drop[r_index] <- j
          r_index <- r_index+1
        }
      }
    }
  }
  #Just in case none to drop
  if(!is.na(NA)){
  x <- x[-drop]
  }
  x_full <- c(0,x,1)


  #Calculate the midpoints of the remaining points
  midpoints <- NA
  for (i in 1:(length(x)+1)){
    midpoints[i] <- (x_full[{i+1}]-x_full[{i}])/2+x_full[{i}]
  }

  #Calculate the widths between each pair of points
  widths <- NA
  for (i in 1:(length(x))){
    widths[i] <- abs(midpoints[{i+1}]-midpoints[{i}])
  }

  #Covariance matirx of proposal distribution
  cov_x <- matrix(0, nrow=length(x), ncol=length(x))
  for (i in 1:length(x)){
    for(j in 1:length(x)){
      cov_x[i,j] <- s2^2*exp(-c2*(x[i]-x[j])^2)
    }
  }

  #For intvertability
  for (i in 1:ncol(cov_x)){
    cov_x[i,i] <- cov_x[i,i]+.0000000000001
  }

  ####SVD#####################
  n <- 1
  mu <- rep(0, length(x))
  tol = 1e-06
  p <- length(mu)
  eS <- svd(cov_x)
  ev <- eS$d
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")

  mvrnorm_calc <- eS$u %*% diag(sqrt(pmax(ev, 0)), p)
  ################################

  #Covariance matrix for prior distribution
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

  #Use SVD
  cov_prior_inv <- svd.inverse(cov_prior)
  G <- matrix(0, nrow=MH_N, ncol=length(x))

  #G[1,] <- mvrnorm(n = 1, mu=rep(0,length(x)), Sigma=cov_prior, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  G[1,] <-rep(0, ncol(G))

  #Collect acceptance rate
  ar_G <- NA

  #Get indeces for which parameters are actual data verses part of the grid
  which_x <- which(names(x)=="data")

  #Metropolis-Hastings algorithm
  for (i in 2:MH_N){

    #Simulate length(x) values from standard normal
    rw <- matrix(rnorm(p * n), n)
    rw <- drop(mu) + mvrnorm_calc %*% t(rw)

    #Add to previous iteration
    g_prop <- G[{i-1},]+rw

    #Calculate log of numerator and denominator of acceptance probability
    ap_num_log <- -length(y)*log(sum(widths*exp(g_prop)))+sum(y_d$n*g_prop[which_x])-.5*t(g_prop)%*%cov_prior_inv%*%g_prop
    ap_den_log <- -length(y)*log(sum(widths*exp(G[{i-1},])))+sum(y_d$n*G[{i-1},which_x])-.5*t(G[{i-1},])%*%cov_prior_inv%*%G[{i-1},]

    #Calculate log of acceptance probability
    ap_log <- ap_num_log-ap_den_log

    if (log(runif(1))<ap_log){
      G[i,] <- g_prop #Acceptance proposal
      ar_G[i] <- 1
    } else {
      G[i,] <- G[i-1,] #Reject proposal
      ar_G[i] <- 0
    }

  }


  #Acceptance Rates
  ar_G <- ar_G[-1]
  print("Acceptance Rate: ")
  print(sum(ar_G)/length(ar_G))
  #print(sum(ar_G)/length(ar_G))
  return(list(G,widths,x,y_orig,sum(ar_G)/length(ar_G)))
}
