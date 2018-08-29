#' Bayesina Density Estimation Using GOLD (Gaussian process On a Log Density).
#'
#' Estimates a density Bayesian techniques using a Gaussian process on a log density
#'
#' @param y A numberic vector to estimate the density on.
#' @param s1 The standard deviation of the Gaussian prior.
#' @param c1 The correlation parameter of the Gaussian prior that controls the correlation in the covariance structure.
#' @param s2 The standard deviation of the Gaussian proposal distribution.
#' @param c2 The correlation parameter that controls the covariance structure of the Gaussian proposal distribution.
#' @param MH_N the number of iterations to run in the algorithm.
#' @param graves An option to have the variance of the proposal distribution chosen using an automatice step size selection (insert reference) (DEFAULT is FALSE).
#'
#' @export
#' @importFrom matrixcalc svd.inverse
#' @importFrom MASS mvrnorm
#'
#' @details
#'
#' The density being estimated takes the form below:
#'
#' \deqn{f(x) =}
#' \deqn{exp(g(x)) / (\int_{0}^{1} exp(g(u)) du)}
#'
#' where g(x) is an unknown log density.
#'
#' Given that g(x) is unknown, the normalizing constant is estimated using a weighted average and the set of unknown paramters that recieve a prior is g(x) at a finite set of points.
#'
#' @return
#'
#' \code{gold} returns a list containing the density estimate results.
#'
#' \item{\code{G}}{A matrix with the posterior draws for g(x) at a finite number of points. The number of rows is the number of iterations.}
#' \item{\code{widths}}{A vector with the widths around each each x used for g(x).}
#' \item{\code{x}}{A vector with the points at which g(x) was estimated.}
#' \item{\code{y}}{A vector containing the data introduced to \code{gold} for density estimation.}
#' \item{\code{ar}}{The acceptance rate from the random walk proposals.}
#' \item{\code{s2}}{If graves==TRUE, a sixth element is returned that contains the standard deviation chosen by the graves method.}
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(5, 1) distribution wiht 200 data points.
#' # Based on the rule of thumb in the details, 7 cutpoints are used
#' y <- rbeta(200, 5, 1)
#' gold_beta <- gold(y, s1 = 1, c1 = 0.8, s2 = 0.4, c2 = 0.5, MH_N = 20000)
#'
#' #Check traceplots
#' gold_plot(gold_beta, type="pdf", data=TRUE)
#'
#' #Examine estimate of PDF
#' gold_plot(gold_beta, type="pdf", data=TRUE)
#'
#' #Examine estimate of CDF
#' gold_plot(gold_beta, type="cdf")
#'
#' #Find probability of being less than 0.4
#' gold_cdf(c(.4), gold_beta)$cdf


gold <- function(y, s1, c1, s2, c2, MH_N = 20000, graves=FALSE){
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
  if(!is.na(drop[1])){
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


  ###GRAVES################

  if(graves==TRUE){
  s2_graves <- c(.001,.01,.03,.06,.08,.1,.2,.3,.6,.8,1,3,5)
  n_step_size <- length(s2_graves)
  target <- .35
  ar <- NA
  which_x <- which(names(x)=="data")

  for (i in 1:n_step_size){
    ar_G <- NA
    #Covariate matrix
    cov_x <- matrix(0, nrow=length(x), ncol=length(x))
    for (l in 1:length(x)){
      for(j in 1:length(x)){
        cov_x[l,j] <- s2_graves[i]*exp(-c2*(x[l]-x[j])^2)
      }
    }

    for (l in 1:ncol(cov_x)){
      cov_x[l,l] <- cov_x[l,l]+.0000000000001
    }


    n <- 1
    mu <- rep(0, length(x))
    tol = 1e-06
    p <- length(mu)
    if (!all(dim(cov_x) == c(p, p)))
      stop("incompatible arguments")
    eS <- svd(cov_x)
    ev <- eS$d
    if (!all(ev >= -tol * abs(ev[1L])))
      stop("'Sigma' is not positive definite")

    mvrnorm_calc <- eS$u %*% diag(sqrt(pmax(ev, 0)), p)

    tune<- matrix(0, nrow=201, ncol=ncol(cov_x))
    tune[1,] <- G[1,]

    for(j in 2:101){
      #Simulate length(x) values from standard normal
      rw <- matrix(rnorm(p * n), n)

      rw <- drop(mu) + mvrnorm_calc %*% t(rw)
      g_prop <- tune[{j-1},]+rw
      #plot(x, exp(g_prop))

      ap_num_log <- -length(y)*log(sum(widths*exp(g_prop)))+sum(g_prop[which_x])-.5*t(g_prop)%*%cov_prior_inv%*%g_prop
      ap_den_log <- -length(y)*log(sum(widths*exp(tune[{j-1},])))+sum(tune[{j-1},which_x])-.5*t(tune[{j-1},])%*%cov_prior_inv%*%tune[{j-1},]

      ap_log <- ap_num_log-ap_den_log

      if (log(runif(1))<ap_log){
        tune[j,] <- g_prop
        ar_G[j] <- 1
      } else {
        tune[j,] <- tune[j-1,]
        ar_G[j] <- 0
      }
      #print(i)
    }
    ar_G <- ar_G[-1]

    ar[i] <- sum(ar_G)
  }

  test <- glm(cbind(ar,100-ar)~log(s2_graves), family="binomial")
  s2_start <- exp((log(target/(1-target))-test$coefficients[1])/test$coefficients[2])

  names(s2_start) <- NULL
  if(s2_start>15 | s2_start<.0001){
    s2_start <- s2_graves[which.min(abs(ar - 30)) ]
  }
  ##########
  s2 <- s2_start
}

cov_x <- matrix(0, nrow=length(x), ncol=length(x))
for (i in 1:length(x)){
  for(j in 1:length(x)){
    cov_x[i,j] <- s2*exp(-c2*(x[i]-x[j])^2)
  }
}

for (i in 1:ncol(cov_x)){
  cov_x[i,i] <- cov_x[i,i]+.0000000000001
}


####SVD#####################
n <- 1
mu <- rep(0, length(x))
tol = 1e-06
p <- length(mu)
if (!all(dim(cov_x) == c(p, p)))
  stop("incompatible arguments")
eS <- svd(cov_x)
ev <- eS$d
if (!all(ev >= -tol * abs(ev[1L])))
  stop("'Sigma' is not positive definite")

mvrnorm_calc <- eS$u %*% diag(sqrt(pmax(ev, 0)), p)



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
  print("Acceptance Rate:")
  print(sum(ar_G)/length(ar_G))
  #print(sum(ar_G)/length(ar_G))
  if(graves==TRUE){
      return(list(G=G,widths=widths,x=x,y=y_orig,ar=sum(ar_G)/length(ar_G), s2=s2))
  }else{
    return(list(G=G,widths=widths,x=x,y=y_orig,ar=sum(ar_G)/length(ar_G)))
  }
}
