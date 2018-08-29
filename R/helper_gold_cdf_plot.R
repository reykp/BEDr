# Plot Probability Density from gold.

helper_gold_cdf_plot <- function(gold_output,burnin=NA, cri=FALSE, data=FALSE){

  if(is.na(burnin)){
    burnin <- nrow(gold_output[[1]])/2
  }



  y_orig <- gold_output$y
  x_gold <- gold_output$x
  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if(min_y<0 | max_y>1){
    input <- (x_gold)*(max_y+.00001-(min_y-.00001))+(min_y-.00001);
    duos_CDF <- gold_cdf(input,gold_output,burnin)
  }else{
    input <- x_gold
    duos_CDF <- gold_cdf(x_gold,gold_output,burnin)
  }



  #Get x
  x <- duos_CDF$x
  CDF <- duos_CDF$cdf

  #Get data
  plot_CDF <- data.frame(x,CDF)

  #Get credible intervals
  crdble <- data.frame(duos_CDF$cri)
  names(crdble) <- c("lower", "upper")
  crdble$x <- duos_CDF$x
  # if(min_y<0 | max_y>1){
  #   crdble$x <- crdble$x*(max_y+.00001-(min_y-.00001))+(min_y-.00001)
  #   crdble$lower <- crdble$lower/(max_y+.00001-(min_y-.00001))
  #   crdble$upper <- crdble$upper/(max_y+.00001-(min_y-.00001))
  # }

  data_y <- data.frame(gold_output$y)
  names(data_y) <- "data"

  g <- ggplot()+
    theme(axis.title = element_text(size = 12))+
    theme_bw()+expand_limits(y=0)
  if(data==TRUE){
    ECDF <- ecdf(data_y$data)
    ECDF_data <- data.frame(data_y$data, ECDF(data_y$data))
    names(ECDF_data) <- c("x", "ECDF")

    g <- g+geom_line(data=ECDF_data, aes(x,ECDF), color="black", size=.6)
  }
  if(cri==TRUE){
    g <- g+geom_line(data=crdble, aes(x,lower), color="red",size=.6)+
      geom_line(data=crdble, aes(x,upper), color="red", size=.6)
  }

  g+geom_line(data=plot_CDF, aes(x, CDF),color="blue", size=.8)+ylab("CDF Estimate")+
    xlab("X")
}
