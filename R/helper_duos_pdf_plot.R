# Plot Probability Density from DUOS.

helper_duos_pdf_plot <- function(duos_output,burnin=NA, cr_i=FALSE, data=FALSE){

  if(is.na(burnin)){
    burnin <- nrow(duos_output[[1]])/2
  }

  y_orig <- duos_output[[3]]
  min_y <- min(y_orig)
  max_y <- max(y_orig)
  if(min_y<0 | max_y>1){
    x <- (1:999/1000)*(max_y+.00001-(min_y-.00001))+(min_y-.00001);
    duos_density <- duos_pdf(x,duos_output,burnin)
  }else{
    duos_density <- duos_pdf(1:999/1000,duos_output,burnin)
  }



  #Get x
  x <- duos_density[[4]]
  PDF <- duos_density[[2]]



  #Get data
  plot_density <- data.frame(x,PDF)

  #Get credible intervals
  crdble <- data.frame(duos_density[[3]])
  names(crdble) <- c("lower", "upper")
  crdble$x <- duos_density[[4]]
  # if(min_y<0 | max_y>1){
  #   crdble$x <- crdble$x*(max_y+.00001-(min_y-.00001))+(min_y-.00001)
  #   crdble$lower <- crdble$lower/(max_y+.00001-(min_y-.00001))
  #   crdble$upper <- crdble$upper/(max_y+.00001-(min_y-.00001))
  # }

  data_y <- data.frame(duos_output[[3]])
  names(data_y) <- "data"

  g <- ggplot(data_y, aes(x=data))+
    theme(axis.title = element_text(size = 12))+
    theme_bw()+expand_limits(y=0)


  if(data==TRUE){
    g <- g+geom_histogram(aes(y=..density..), fill="grey", color="black")
  }
  if(cr_i==TRUE){
    g <- g+geom_line(data=crdble, aes(x,lower), color="red",size=.6)+
      geom_line(data=crdble, aes(x,upper), color="red", size=.6)
  }

  g+geom_line(data=plot_density, aes(x, PDF),color="blue", size=.8)+ylab("PDF Estimate")+
    xlab("X")
}
