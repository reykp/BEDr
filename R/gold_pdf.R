#' Bayesina Density Estimate of PDF from GOLD.
#'
#' Calculates an estimate of the density based on the output from \code{duos} for an individual
#' or vector of values
#'
#' @param x A single value or vector of values to calculate the density at
#' @param gold_output The list returned by \code{duos}
#' @param burnin The desired burnin to discard from including in the estimate
#' @param scale If the data is not already between 0 or 1, this determines if the pdf is returned on the (0,1) transformed data or the original data scale
#'
#' @export

gold_pdf <- function(x, gold_output, burnin=NA,scale=FALSE){


  G <- gold_output[[1]]
  widths <- gold_output[[2]]
  x_gold <- gold_output[[3]]
  y_orig <- gold_output[[4]]

  input <- x

  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if((min_y<0 | max_y>1)){
    input_scaled <- (x-(min_y-.00001))/(max_y+.00001-(min_y-.00001))
  }else{
    input_scaled <- input
  }

  if(is.na(burnin)){
    burnin <- nrow(G)/2
  }

  if(burnin>nrow(G)){
    stop("The specified burnin is greater than the number of iterations.")
  }

  G_burnin <- G[burnin:nrow(G),]

  #Exponentiate for numerator
  G_burnin_exp <- exp(G_burnin)
  #Calculate normalizing constat
  nc <- (widths%*%t(G_burnin_exp))

  #Calcualte pdf at each iteration
  G_PDF <- G_burnin_exp
  for (i in 1:nrow(G_burnin_exp)){
    G_PDF[i,] <- G_PDF[i,]/nc[i]
  }

  #Calculate posterior mean pdf
  pdf_y <- NA
  for (j in 1:ncol(G_PDF)){
    pdf_y[j] <- mean(G_PDF[,j])
  }

  #Calculate percentiles

  pdf_y_perc <- apply(G_PDF, 2, quantile, probs=c(.025, .975))

  G_PDF_return <- matrix(0, nrow=nrow(G_PDF), ncol=length(input))
  pdf_y_return <- NA
  pdf_y_perc_return <- matrix(0, nrow=length(input), ncol=2)

  x_pdfs <- NA
  for(i in 1:length(input)){
    if((input[i]<=max_y)&(input[i]>=min_y)){
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
    return(list(pdf_matrix=G_PDF_return, pdf_y=pdf_y_return, pdf_percentiles=pdf_y_perc_return, x=input))


  }else{
    return(list(pdf_matrix=G_PDF_return, pdf_y=pdf_y_return, pdf_percentiles=pdf_y_perc_return, x=input_scaled))
  }


}
