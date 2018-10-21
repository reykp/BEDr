# Plot CDF from DUOS.



helper_duos_cdf_plot <- function(duos_output, estimate, burnin=NA, cri=FALSE, data=FALSE, interact){

  C <- duos_output$C
  if(is.na(burnin)){
    burnin <- ceiling(nrow(duos_output$C)/2)
  }

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
  
  scale_l <- duos_output$scale_l
  scale_u <- duos_output$scale_u
  

  y_orig <- duos_output$y
  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if(min_y<0 | max_y>1){
    input <- (1:999/1000)*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l);
    duos_CDF <- duos_cdf(input,duos_output,burnin, estimate = estimate)
  }else{
    input <- 1:999/1000
    duos_CDF <- duos_cdf(input,duos_output,burnin, estimate = estimate)
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


   data_y <- data.frame(duos_output$y)
    names(data_y) <- "data"
  

  g <- ggplot()+
    theme(axis.title = element_text(size = 12))+
    theme_bw()+expand_limits(y=0)


  if(data==TRUE){
    ECDF <- ecdf(data_y$data)
    ECDF_data <- data.frame(data_y$data, ECDF(data_y$data))
    names(ECDF_data) <- c("x", "ECDF")

    g <- g+geom_line(data=ECDF_data, aes(x,ECDF), color="black", size=.8)
  }
  if(cri==TRUE){
    g <- g+geom_line(data=crdble, aes(x,lower), color="red",size=.6)+
      geom_line(data=crdble, aes(x,upper), color="red", size=.6)
  }

  if(interact == TRUE){
    
    suppressMessages(plotly::ggplotly(g+geom_line(data=plot_CDF, aes(x, CDF),color="blue", size=.8)+ylab("CDF")+
                       xlab("X")))
    
  }else{
    g+geom_line(data=plot_CDF, aes(x, CDF),color="blue", size=.8)+ylab("CDF")+
      xlab("X")
  }
}
