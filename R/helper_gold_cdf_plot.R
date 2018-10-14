# Plot Probability Density from gold.

helper_gold_cdf_plot <- function(gold_output,burnin=NA, cri=FALSE, data=FALSE, interact){

  if(is.na(burnin)){
    burnin <- nrow(gold_output[[1]])/2
  }



  y_orig <- gold_output$y
  x_gold <- gold_output$x
  
  
  if(!is.null(gold_output[["poi"]])){
    max_y <- max(y, gold_output$poi)
    min_y <- min(y, gold_output$poi)
  }else{
    min_y <- min(y_orig)
    max_y <- max(y_orig)
  }
  
  scale_l <- gold_output$scale_l
  scale_u <- gold_output$scale_u
  
  if(min_y<0 | max_y>1){
    input <- (x_gold)*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l);
    gold_CDF <- gold_cdf(input,gold_output,burnin)
  }else{
    input <- x_gold
    gold_CDF <- gold_cdf(x_gold,gold_output,burnin)
  }



  #Get x
  x <- gold_CDF$x
  CDF <- gold_CDF$cdf

  #Get data
  plot_CDF <- data.frame(x,CDF)

  #Get credible intervals
  crdble <- data.frame(gold_CDF$cri)
  names(crdble) <- c("lower", "upper")
  crdble$x <- gold_CDF$x
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

  if(interact == TRUE){
    suppressMessages(plotly::ggplotly(g+geom_line(data=plot_CDF, aes(x, CDF),color="blue", size=.8)+ylab("CDF Estimate")+
      xlab("X")))
  }else{
    g+geom_line(data=plot_CDF, aes(x, CDF),color="blue", size=.8)+ylab("CDF Estimate")+
      xlab("X")
  }
    
  }

