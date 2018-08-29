#' Estimate of PDF from GOLD
#'
#' Calculates a posterior estimate of the PDF based on the output from \code{gold} for an individual
#' or vector of values.
#'
#' @param x A single value or vector of values at which to calculate the PDF. These values are to be entered on the scale of the data (i.e. values can fall outside of 0 and 1).
#' @param gold_output The list returned by \code{gold} containing the density estimate results.
#' @param burnin The desired burnin to discard from the results. If no values is entered, the default is half the number of iterations.
#' @param scale This value TRUE/FALSE indicates whether to return scaled or unscalled results IF the original data does not fall between 0 and 1. The default is FALSE (i.e. returns results on the original data scale).
#'
#' @export
#'
#' @details
#'
#' The function \code{gold_pdf} returns the posterior mean PDF. The PDF is calculated based on the following equation at each iteration:
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
#' \code{gold_pdf} returns a list of the PDF results from \code{gold}.
#'
#' \item{\code{pdf}}{A vector of the posterior mean PDF values at each value in \code{x}.}
#' \item{\code{cri}}{A matrix with 2 columns and rows equaling the length of \code{x} containing the 95\% credible interval for the PDF at each of the points in \code{x}.}
#' \item{\code{mat}}{A matrix containing the PDF values for each \code{x} at EACH itertation after the burnin is discarded. The number of columns is the length of \code{x}.}
#' \item{\code{x}}{A vector containing the values at which to estimate the PDF. If the data is not between 0 and 1 and scale=TRUE, the scaled version of \code{x} is returned.}
#'
#' @examples
#'
#' ## --------------------------------------------------------------------------------
#' ## Beta Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'duos' on data sampled from a Beta(2,5) distribution wiht 100 data points.
#' y <- rbeta(100, 2, 5)
#' gold_beta <- gold(y, s1 = 1, c1 = 1, s2 = 0.8, c2 = 0.8, MH_N = 20000)
#' #Calculate pdf at a variety of values
#' pdf_beta <- gold_pdf(x = c(.01, .25, .6, .9), gold_beta)
#'
#' #Examine the PDF at 'x'
#' pdf_beta$pdf
#'
#' #Examine the credibal intervals of the PDF at 'x'
#' pdf_beta$cri
#'
#' ## --------------------------------------------------------------------------------
#' ## Normal Distribution
#' ## --------------------------------------------------------------------------------
#'
#' # First run 'gold' on data sampled from a Normal(0,1) distribution with 200 data points.
#' y <- rnorm(200, 0, 1)
#' gold_norm <- gold(y, s1 = 1, c1 = 0.8, s2 = 1, c2 = 0.5, MH_N = 20000)
#' pdf_norm <- gold_pdf(x = c(-2, -1, 0, 0.8, 1.8), gold_norm)
#'
#' #Examine the PDF at 'x'
#' pdf_norm$pdf
#'
#' #Examine the credibal intervals of the CDF at 'x'
#' pdf_norm$cri
#'
#' Histogram of distribution of the PDF density estimate at 0.8
#' hist(pdf_norm$mat[, 4])


gold_pdf <- function(x, gold_output, burnin=NA,scale=FALSE){

  #Get parameters
  G <- gold_output$G
  #Get widths around points estimate parameters at
  widths <- gold_output$widths
  #Get points estimated density at
  x_gold <- gold_output$x
  #Get original data
  y_orig <- gold_output$y

  #Assign points to find pdf at to input
  input <- x

  #Find min and max of data
  min_y <- min(y_orig)
  max_y <- max(y_orig)

  #IF original data lies outside of 0 and 1, needs to be scaled
  if((min_y<0 | max_y>1)){
    input_scaled <- (x-(min_y-.00001))/(max_y+.00001-(min_y-.00001))
  }else{
    input_scaled <- input
  }

  #if no burnin assigned, use half of iterations as burnin
  if(is.na(burnin)){
    burnin <- nrow(G)/2
  }

  #Make sure burnin in less than interations
  if(burnin>nrow(G)){
    stop("The specified burnin is greater than the number of iterations.")
  }

  #Get iterations without burnin
  G_burnin <- G[burnin:nrow(G),]

  #Exponentiate for numerator
  G_burnin_exp <- exp(G_burnin)
  #Calculate normalizing constat
  nc <- (widths%*%t(G_burnin_exp))
  #nc <- matrix(nc, ncol=1)

  #Calcualte pdf at each iteration
  G_PDF <- G_burnin_exp
  for(i in 1:nrow(G_PDF)){
    G_PDF[i,] <- G_PDF[i,]/nc[i]
  }

  #Calculate posterior mean pdf
  pdf_y <- NA
  for (j in 1:ncol(G_PDF)){
    pdf_y[j] <- mean(G_PDF[,j])
  }

  #Calculate percentiles
  pdf_y_perc <- apply(G_PDF, 2, quantile, probs=c(.025, .975))

  #The next steps are to find estimates of densitity a locations that lie closest to the values in x
  G_PDF_return <- matrix(0, nrow=nrow(G_PDF), ncol=length(input))
  pdf_y_return <- NA
  pdf_y_perc_return <- matrix(0, nrow=length(input), ncol=2)

  x_pdfs <- NA
  for(i in 1:length(input)){
    #Check if any of the input outside of the max and min of x if needs to be scaled
    if((input_scaled[i]<1)&(input_scaled[i]>0)){
      diff1 <- input_scaled[i]-x_gold[max(which(x_gold<=input_scaled[i]))]
      if(!is.na(x_gold[max(which(x_gold<=input_scaled[i]))+1])){
        diff2 <- x_gold[max(which(x_gold<=input_scaled[i]))+1]-input_scaled[i]
      }else{
        diff2 <- 0
      }

      if(diff1<=diff2){
        x_pdfs[i] <- max(which(x_gold<=input_scaled[i]))
        G_PDF_return[,i] <- G_PDF[,x_pdfs[i]]
        pdf_y_return[i] <- pdf_y[x_pdfs[i]]
        pdf_y_perc_return[i,] <- t(pdf_y_perc)[x_pdfs[i],]
      }else{
        x_pdfs[i] <- max(which(x_gold<=input_scaled[i]))+1
        G_PDF_return[,i] <- G_PDF[,x_pdfs[i]]
        pdf_y_return[i] <- pdf_y[x_pdfs[i]]
        pdf_y_perc_return[i,] <- t(pdf_y_perc)[x_pdfs[i],]
      }
    }else{
      if(input[i]<min_y){
        G_PDF_return[,i] <- rep(0,nrow(G_PDF_return))
        pdf_y_return[i] <- 0
        pdf_y_perc_return[i,] <- c(0,0)
      }else{
        G_PDF_return[,i] <- rep(1,nrow(G_PDF_return))
        pdf_y_return[i] <- 1
        pdf_y_perc_return[i,] <- c(1,1)
      }
    }
  }


  if(scale==FALSE & (min_y<0 | max_y>1)){
    #To have pdf on original scale, need to transform estimate of pdf
    pdf_y_return <- pdf_y_return/(max_y+.00001-(min_y-.00001))
    pdf_y_perc_return[,1] <- pdf_y_perc_return[,1]/(max_y+.00001-(min_y-.00001))
    pdf_y_perc_return[,2] <- pdf_y_perc_return[,2]/(max_y+.00001-(min_y-.00001))
    G_PDF_return <- G_PDF_return/(max_y+.00001-(min_y-.00001))

    return(list(pdf=pdf_y_return, cri=pdf_y_perc_return, mat=G_PDF_return, x=input))

  }else{
    return(list(pdf=pdf_y_return, cri=pdf_y_perc_return, mat=G_PDF_return, x=input_scaled))
  }


}
