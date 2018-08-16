# Plot Probability Density from DUOS.

helper_duos_pdf_plot <- function(duos_output,burnin=NA, cri=FALSE, data=FALSE){

  if(is.na(burnin)){
    burnin <- nrow(duos_output$C)/2
  }

  y_orig <- duos_output$y
  min_y <- min(y_orig)
  max_y <- max(y_orig)
  if(min_y<0 | max_y>1){
    x <- (1:999/1000)*(max_y+.00001-(min_y-.00001))+(min_y-.00001);
    duos_density <- duos_pdf(x,duos_output,burnin)
  }else{
    duos_density <- duos_pdf(1:999/1000,duos_output,burnin)
  }



  #Get x
  x <- duos_density$x
  PDF <- duos_density$pdf



  #Get data
  plot_density <- data.frame(x,PDF)

  #Get credible intervals
  crdble <- data.frame(duos_density$cri)
  names(crdble) <- c("lower", "upper")
  crdble$x <- duos_density$x


  data_y <- data.frame(duos_output$y)
  names(data_y) <- "data"

  g <- ggplot(data_y, aes(x=data))+
    theme(axis.title = element_text(size = 12))+
    theme_bw()+expand_limits(y=0)


  if(data==TRUE){
    g <- g+geom_histogram(aes(y=..density..), fill="grey", color="black")
  }
  if(cri==TRUE){
    g <- g+geom_line(data=crdble, aes(x,lower), color="red",size=.6)+
      geom_line(data=crdble, aes(x,upper), color="red", size=.6)
  }

  g+geom_line(data=plot_density, aes(x, PDF),color="blue", size=.8)+ylab("PDF Estimate")+
    xlab("X")
}
