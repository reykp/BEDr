# Plot CDF from DUOS.



helper_duos_cdf_plot <- function(duos_output,burnin=NA, cri=FALSE, data=FALSE){

  if(is.na(burnin)){
    burnin <- nrow(duos_output$C)/2
  }

  scale_l <- duos_output$scale_l
  scale_u <- duos_output$scale_u
  

  y_orig <- duos_output$y
  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if(min_y<0 | max_y>1){
    input <- (1:999/1000)*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l);
    duos_CDF <- duos_cdf(input,duos_output,burnin)
  }else{
    input <- 1:999/1000
    duos_CDF <- duos_cdf(input,duos_output,burnin)
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


  data_y <- data.frame(y_orig)
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
