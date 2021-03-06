# Plot Probability Density from gold.

helper_gold_pdf_plot <- function(gold_output,burnin=NA, cri=FALSE, data=FALSE, interact){
  
  if(is.na(burnin)){
    burnin <- nrow(gold_output$G)/2
  }



  y_orig <- gold_output$y
  x_gold <- gold_output$x
  if(!is.null(gold_output[["poi"]])){
    max_y <- max(y_orig, gold_output$poi)
    min_y <- min(y_orig, gold_output$poi)
  }else{
    min_y <- min(y_orig)
    max_y <- max(y_orig)
  }
  
  scale_l <- gold_output$scale_l
  scale_u <- gold_output$scale_u

  if(min_y<0 | max_y>1){
    input <- x_gold*(max_y+scale_u-(min_y-scale_l))+(min_y-scale_l);
    gold_density <- gold_pdf(input,gold_output,burnin)
  }else{
    input <- x_gold
    gold_density <- gold_pdf(input,gold_output,burnin)
  }



  #Get x
  x <- gold_density$x
  PDF <- gold_density$pdf

  #Get data
  plot_density <- data.frame(x,PDF)

  #Get credible intervals
  crdble <- data.frame(gold_density$cri)
  names(crdble) <- c("lower", "upper")
  crdble$x <- gold_density$x
  # if(min_y<0 | max_y>1){
  #   crdble$x <- crdble$x*(max_y+.00001-(min_y-.00001))+(min_y-.00001)
  #   crdble$lower <- crdble$lower/(max_y+.00001-(min_y-.00001))
  #   crdble$upper <- crdble$upper/(max_y+.00001-(min_y-.00001))
  # }

    data_y <- data.frame(gold_output$y)
    names(data_y) <- "X"
  

  g <- ggplot(data_y, aes(x=X))+
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
    if(data==FALSE){
      suppressMessages(plotly::ggplotly(g+geom_line(data=plot_density, aes(x, PDF),color="blue", size=.8)+ylab("Density")+
      xlab("X"), tooltip = c("X", "PDF", "lower", "upper")))
    }else{
      suppressMessages(plotly::ggplotly(g+geom_line(data=plot_density, aes(x, PDF),color="blue", size=.8)+ylab("Density")+
                                          xlab("X")))
    }
  }else{
    g+geom_line(data=plot_density, aes(x, PDF),color="blue", size=.8)+ylab("Density")+
      xlab("X")
  }
}
