#' Estimate of PDF from GOLD
#'
#' Calculates a posterior estimate of the PDF based on the output from \code{gold} for an individual
#' or vector of values.
#'
#' @param x A single value or vector of values at which to calculate the PDF. These values are to be entered on the scale of the data (i.e. values can fall outside of 0 and 1).
#' @param gold_output The list returned by \code{gold} containing the density estimate results.
#' @param burnin The desired burnin to discard from the results. If no values is entered, the default is half the number of iterations.
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
#' \item{\code{x}}{A vector containing the values at which to estimate the PDF.}
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


gold_pdf <- function(x, gold_output, burnin = NA){

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
  if(!is.null(gold_output[["poi"]])){
    max_y <- max(y_orig, gold_output$poi)
    min_y <- min(y_orig, gold_output$poi)
  }else{
    min_y <- min(y_orig)
    max_y <- max(y_orig)
  }

  scale_l <- gold_output$scale_l
  scale_u <- gold_output$scale_u
  
  
  # Check to make sure can estimate p at x
  if(max(gold_output$y)>1 | min(gold_output$y)< 0){
    
    outside_range_u <- which(x > (max(gold_output$y)+gold_output$scale_u))
    outside_range_l <- which(x < (min(gold_output$y)-gold_output$scale_l))
    if((length(outside_range_u)>0) &(length(outside_range_l)==0)){
      print(x[outside_range_u])
      stop("The requested 'x' vector contains the above value(s) outside the range of max(y)+scale_u.")
    }
    
    if((length(outside_range_l)>0) &(length(outside_range_u)==0)){
      print(x[outside_range_l])
      stop("The requested 'x' vector contains the above value(s) outside the range of min(y)-scale_l.")
    }
    
    if((length(outside_range_l)>0) &(length(outside_range_u)>0)){
      print(x[c(outside_range_l, outside_range_u)])
      stop("The requested 'x' vector contains the above value(s) outside the range of min(y)-scale_l and max(y)+scale_u.")
    }
    
    
  }else{
    outside_range_u <- which(x > 1)
    outside_range_l <- which(x < 0)
    if((length(outside_range_u)>0)|(length(outside_range_l)>0)){
      print(x[c(outside_range_u, outside_range_l)])
      stop("The requested 'x' vector contains the above value(s) outside the range of (0, 1).")
    }
    
  }
  
  # new_upper <- ifelse(max_y>=0, max_y+2*sd(y), max_y-2*sd(y))
  # new_lower <- ifelse(min_y>=0, min_y+2*sd(y), min_y-2*sd(y))

  #IF original data lies outside of 0 and 1, needs to be scaled
  if((min_y<0 | max_y>1)){
    input_scaled <- (x-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l))
  }else{
    input_scaled <- input
  }

  # if((min_y<0 | max_y>1)){
  #   input_scaled <- (x-new_lower)/(new_upper-new_lower)
  # }else{
  #   input_scaled <- input
  # }

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

  #For extrapolating
  # 
  # ss_x <<- x_gold
  # 
  # 
  # ss_function <- function(x){
  #   ss_matrix <- smooth.spline(ss_x, x)
  #   return(predict(ss_matrix, x_ss)$y) 
  # }
  
  for(i in 1:length(input_scaled)){
    
    find_x <- which(x_gold==input_scaled[i]) 
    
    if((min_y<0|max_y>1) & (input[i]<(min_y-scale_l)|input[i]>(max_y+scale_u))){
          print(input[i])
          print("is out of the range of the data in 'y'.")
        }else if(length(find_x)>0){
      pdf_y_return[i] <- pdf_y[find_x]
      G_PDF_return[,i] <- G_PDF[,find_x]
      pdf_y_perc_return[i,] <- t(pdf_y_perc)[find_x,]
    }else{
      #x_ss <<- input_scaled[i]
      ss_indiv <- smooth.spline(x_gold, pdf_y)
      pdf_y_return[i] <- predict(ss_indiv, input_scaled[i])$y
      #G_PDF_return[,i] <- apply(G_PDF, 1, ss_function)
      G_PDF_return[,i] <- rep(NA, nrow(G_PDF_return))
      #pdf_y_perc_return[i,] <- quantile(G_PDF_return[,i], c(.025, .975))
      ss_indiv_q1 <- smooth.spline(x_gold, t(pdf_y_perc)[,1])
      pdf_y_perc_return[i,1] <- predict(ss_indiv_q1, input_scaled[i])$y
      ss_indiv_q2 <- smooth.spline(x_gold, t(pdf_y_perc)[,2])
      pdf_y_perc_return[i,2] <- predict(ss_indiv_q2, input_scaled[i])$y
    }
  }
  
  
  
  # x_pdfs <- NA
  # for(i in 1:length(input)){
  #   #Check if any of the input outside of the max and min of x if needs to be scaled
  #   if((input_scaled[i]<1)&(input_scaled[i]>0)){
  #     diff1 <- input_scaled[i]-x_gold[max(which(x_gold<=input_scaled[i]))]
  #     if(!is.na(x_gold[max(which(x_gold<=input_scaled[i]))+1])){
  #       diff2 <- x_gold[max(which(x_gold<=input_scaled[i]))+1]-input_scaled[i]
  #     }else{
  #       diff2 <- 0
  #     }
  # 
  #     if(diff1<=diff2){
  #       x_pdfs[i] <- max(which(x_gold<=input_scaled[i]))
  #       G_PDF_return[,i] <- G_PDF[,x_pdfs[i]]
  #       pdf_y_return[i] <- pdf_y[x_pdfs[i]]
  #       pdf_y_perc_return[i,] <- t(pdf_y_perc)[x_pdfs[i],]
  #     }else{
  #       x_pdfs[i] <- max(which(x_gold<=input_scaled[i]))+1
  #       G_PDF_return[,i] <- G_PDF[,x_pdfs[i]]
  #       pdf_y_return[i] <- pdf_y[x_pdfs[i]]
  #       pdf_y_perc_return[i,] <- t(pdf_y_perc)[x_pdfs[i],]
  #     }
  #   }else{
  #     if(input[i]<min_y){
  #       G_PDF_return[,i] <- rep(0,nrow(G_PDF_return))
  #       pdf_y_return[i] <- 0
  #       pdf_y_perc_return[i,] <- c(0,0)
  #     }else{
  #       G_PDF_return[,i] <- rep(1,nrow(G_PDF_return))
  #       pdf_y_return[i] <- 1
  #       pdf_y_perc_return[i,] <- c(1,1)
  #     }
  #   }
  # }

  if((min_y<0 | max_y>1)){
    #To have pdf on original scale, need to transform estimate of pdf
    pdf_y_return <- pdf_y_return/(max_y+scale_u-(min_y-scale_l))
    pdf_y_perc_return[,1] <- pdf_y_perc_return[,1]/(max_y+scale_u-(min_y-scale_l))
    pdf_y_perc_return[,2] <- pdf_y_perc_return[,2]/(max_y+scale_u-(min_y-scale_l))
    G_PDF_return <- G_PDF_return/(max_y+scale_u-(min_y-scale_l))

    # pdf_y_return <- pdf_y_return/(new_upper-new_lower)
    # pdf_y_perc_return[,1] <- pdf_y_perc_return[,1]/(new_upper-new_lower)
    # pdf_y_perc_return[,2] <- pdf_y_perc_return[,2]/(new_upper-new_lower)
    # G_PDF_return <- G_PDF_return/(new_upper-new_lower)

    return(list(pdf=pdf_y_return, cri=pdf_y_perc_return, mat=G_PDF_return, x=input))

  }else{
    return(list(pdf=pdf_y_return, cri=pdf_y_perc_return, mat=G_PDF_return, x=input_scaled))
  }


}
