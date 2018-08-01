# Plot Probability Density from gold.

helper_gold_pdf_plot <- function(gold_output,burnin=NA, cr_i=FALSE, data=FALSE){

  if(is.na(burnin)){
    burnin <- nrow(gold_output[[1]])/2
  }



  y_orig <- gold_output[[4]]
  x <- gold_output[[3]]
  min_y <- min(y_orig)
  max_y <- max(y_orig)

  if(min_y<0 | max_y>1){
    input <- x*(max_y+.00001-(min_y-.00001))+(min_y-.00001);
    gold_density <- gold_pdf(x,gold_output,burnin)
  }else{
    gold_density <- gold_pdf(x,gold_output,burnin)
  }



  #Get x
  x <- gold_density[[4]]
  PDF <- gold_density[[2]]

  #Get data
  plot_density <- data.frame(x,PDF)

  #Get credible intervals
  crdble <- data.frame(gold_density[[3]])
  names(crdble) <- c("lower", "upper")
  crdble$x <- gold_density[[4]]
  # if(min_y<0 | max_y>1){
  #   crdble$x <- crdble$x*(max_y+.00001-(min_y-.00001))+(min_y-.00001)
  #   crdble$lower <- crdble$lower/(max_y+.00001-(min_y-.00001))
  #   crdble$upper <- crdble$upper/(max_y+.00001-(min_y-.00001))
  # }

  data_y <- data.frame(gold_output[[4]])
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
