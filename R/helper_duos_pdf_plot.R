# Plot Probability Density from DUOS.

helper_duos_pdf_plot <- function(duos_output,burnin=NA, cri=FALSE, data=FALSE, interact, scale){

  
  if(is.na(burnin)){
    burnin <- ceiling(nrow(duos_output$C)/2)
  }

  C <- duos_output$C
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
  
  y_orig <- duos_output$y
  min_y <- min(y_orig)
  max_y <- max(y_orig)
  
  scale_l <- duos_output$scale_l
  scale_u <- duos_output$scale_u
  
  if(min_y<0 | max_y>1){
    x <- (1:999/1000)*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l);
    duos_density <- duos_pdf(x,duos_output,burnin, scale)
  }else{
    duos_density <- duos_pdf(1:999/1000,duos_output,burnin, scale)
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

  if(scale == TRUE){
    data_y <- data.frame(duos_output$y)
    names(data_y) <- "data"
    data_y$data <- (data_y$data-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l))
  }else{
    data_y <- data.frame(duos_output$y)
    names(data_y) <- "data"
  }

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

  if(interact == TRUE){
    suppressMessages(plotly::ggplotly(g+geom_line(data=plot_density, aes(x, PDF),color="blue", size=.8)+ylab("PDF Estimate")+
                       xlab("X")))
    
  }else{
    g+geom_line(data=plot_density, aes(x, PDF),color="blue", size=.8)+ylab("PDF Estimate")+
      xlab("X")
  }
}
