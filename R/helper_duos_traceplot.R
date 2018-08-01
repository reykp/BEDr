# Plot CDF from DUOS.



helper_duos_traceplot <- function(duos_output,burnin=NA){

  if(is.na(burnin)){
    burnin <- nrow(duos_output[[1]])/2
  }

  duos_CDF <- duos_cdf(1:999/1000,duos_output,burnin)

  #Get x
  x <- duos_CDF[[4]]
  CDF <- duos_CDF[[2]]

  #Get data
  plot_CDF <- data.frame(x,CDF)

  #Get credible intervals
  crdble <- data.frame(duos_CDF[[3]])
  names(crdble) <- c("lower", "upper")
  crdble$x <- 1:999/1000

  data_y <- data.frame(duos_CDF[[5]])
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
  if(cr_i==TRUE){
    g <- g+geom_line(data=crdble, aes(x,lower), color="red",size=.6)+
      geom_line(data=crdble, aes(x,upper), color="red", size=.6)
  }

  g+geom_line(data=plot_CDF, aes(x, CDF),color="blue", size=.8)+ylab("CDF Estimate")
}
