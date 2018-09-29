// [[Rcpp::plugins(cpp11)]]
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;


//' @title duos (Distribution of Uniform Order Statistics)
//' @description
//' Runs Gibbs sampling algorithm to perform Bayesian density estimation through \code{duos}.
//' 
//' @usage
//' 
//' duos(y, k = ceiling(n/50) + 3, N = 20000, alpha = 1, scale_l = 0.00001, scale_u = 0.00001, start = NA)
//' 
//' @name duos
//' @param y A numeric vector. \code{duos} estimates the density on this data.
//' @param k The number of cut-points for the density estimate. There is a default based on the size of \code{y} (see details).
//' @param N The number of iterations to run in the algorithm. The default is 20,000.
//' @param alpha The parameter values for the Dirichlet prior. The parameter is constant and set to 1 by default.
//' @param scale_l A value >= 0 controlling the scaling based on the minimum data value. The default is 0.00001 (see details).
//' @param scale_u A value >= 0 controlling the scaling based on the maximum data value. The default is 0.00001 (see details).
//' @param start A vector of starting values of length k. This is useful for creating output to use in Gelman and Rubin's convergence diagnostics (see examples). If NA, reasonable starting values are selected by the algorithm. 
//' @examples
//'
//' @export
//'
//' @details
//'
//' The density being estimated takes the form below:
//' 
//' \eqn{f(x) =}
//' \deqn{(\pi_1) / (\gamma_1) , 0 \le x < \gamma_1}
//' \deqn{(\pi_2) / (\gamma_2-\gamma_1) , \gamma_1 \le x  < \gamma_2}
//' \deqn{(\pi_3) / (\gamma_3-\gamma_2) , \gamma_2 \le x  < \gamma_3}
//' \deqn{...  , ... \le x  < ...}
//' \deqn{(\pi_k) / (\gamma_k-\gamma_{k-1}) , \gamma_{k-1} \le x  < \gamma_k} 
//' \deqn{(\pi_{k+1}) / (1-\gamma_k) , \gamma_k \le x  < 1}
//'
//' where \eqn{\gamma_1 < \gamma_2 < ... < \gamma_k  is in (0,1) and \pi_1 + \pi_2 + ... + \pi_{k+1} = 1}
//'
//' The prior on the \eqn{\gamma's} is: draw k values from $unif(0,1)$ and order: $0 < \gamma_1 < \gamma_2 < ... < \gamma_{k-1} < \gamma_k < 1$.
//' The prior on the \eqn{\pi's} is Dirichlet(\eqn{\alpha}).
//' 
//' This density operates on data between 0 and 1. Thus, if the input is not between 0 and 1, it is standardized. The formula for scaling is below:
//' 
//' \deqn{(y-(min(y)-scale_l))/(max(y)+scale_u-(min(y)-scale_l))}.
//' 
//' Values of 0 for the scale parameters indicates the density will only be the edges of the data in \code{y}. The default for \code{scale_l} and \code{scale_u} is 0.0001.
//'
//' The recommended number of cut-points starts at 3, and then for each increment of 50, adds a cut-point.
//'
//' Default: k = floor(n/50)+3
//' 
//' @return
//'
//' \code{duos} returns a list of density estimate results.
//'
//' \item{\code{C}}{A matrix containing the posterior draws for the cut-point parameters. The number of columns is \code{k}, and the number of rows is \code{N}.}
//' \item{\code{P}}{A matrix containing the posterior draws for the bin proportion parameters. The number of columns is \code{k}+1, and the number of rows is \code{N}.}
//' \item{\code{y}}{A vector containing the data introduced to \code{duos} for density estimation.}
//' \item{\code{k}}{The number of cut-points.}
//' \item{\code{alpha}}{The parameter for the prior on the bin proportions.}
//' \item{\code{scale_l}}{The scaling parameter for the lower end of the data.}
//' \item{\code{scale_u}}{The scaling parameter for the upper end of the data.}
//'
//' @examples
//' ## --------------------------------------------------------------------------------
//' ## Beta Distribution
//' ## --------------------------------------------------------------------------------
//'
//' # Run 'duos' on data sampled from a Beta(5, 1) distribution with 100 data points.
//' # All defaults are used.
//' y <- rbeta(100, 5, 1)
//' duos_beta <- duos(y)
//'
//' #Examine estimate of PDF
//' duos_plot(duos_beta, type = "pdf", data = TRUE)
//'
//' #Examine estimate of CDF with credible intervals
//' duos_plot(duos_beta, type = "cdf")
//'
//' #Find probability of being less than 0.1
//' duos_cdf(c(.1), duos_beta)$cdf
//' 
//' ## --------------------------------------------------------------------------------
//' ## Claw Distribution
//' ## --------------------------------------------------------------------------------
//' # Run 'duos' on data sampled from a 'Claw' distribution with 200 data points.
//' u <- runif(200)
//' 
//' # Variable to store the samples from the mixture distribution 
//' y <- rep(NA, 200)
//' # Sampling from the mixture
//' for(i in 1:200){
//' if(u[i] < 0.5){
//'      y[i] <- rnorm(1, 0, 1)
//'      }else if(u[i] < 0.6){
//'      y[i] <- rnorm(1, -1, 0.1)
//'      }else if (u[i] < 0.7){
//'      y[i] <- rnorm(1, -0.5, 0.1)
//'      }else if (u[i] < 0.8){
//'      y[i] <- rnorm(1, 0, 0.1)
//'      }else if (u[i] <- 0.9){
//'      y[i] <- rnorm(1, 0.5, 0.1)
//'      }else{
//'      y[i] <- rnorm(1, 1, 0.1)
//'      }
//'      }
//'      
//' # Choose 12 cutpoints
//' duos_claw <- duos(y = y, k = 12, N = 20000)
//'
//' # Examine estimate of PDF
//' duos_plot(duos_claw, type = "pdf", data = TRUE)
//'
//' # Examine estimate of CDF
//' duos_plot(duos_claw, type = "cdf")
//'
//' # Find probability of being greater than 90
//' 1-duos_cdf(c(0.90), duos_claw)$cdf
//' 
//' ## --------------------------------------------------------------------------------
//' ## Normal Distribution
//' ## --------------------------------------------------------------------------------
//' 
//' # Run 'duos' on data sampled from a Normal(0, 1) distribution with 50 data points.
//' y <- rnorm(50, 0, 1)
//' # Run for 15,000 iterations and change the 'scale' values to 1.5*sd(y)
//' duos_norm <- duos(y = y, scale_l = 1.5*sd(y), scale_u = 1.5*sd(y), N = 15000)
//'
//' # Examine estimate of PDF
//' duos_plot(duos_norm, type = "pdf", data = TRUE)
//'
//' # Examine estimate of CDF
//' duos_plot(duos_norm, type = "cdf", cri = TRUE)
//'
//' # Find probability of being greater than 2.5
//' 1-duos_cdf(c(2.5), duos_norm)$cdf
//'
//' ## --------------------------------------------------------------------------------
//' ## Uniform Distribution
//' ## --------------------------------------------------------------------------------
//' 
//' # Below is an example of using the 'start' variable to calculate the Gelman and Rubin diagnostic
//' 
//' # First run 'duos' on data sampled from a Uniform(0, 1) distribution with 50 data points.
//' y <- runif(50)
//' 
//' # Use 4 cutpoints
//' # Create three runs with diverse starting values
//' duos_unif1 <- duos(y = y, k = 4, N = 20000, start = runif(4, 0, 1/3))
//' duos_unif2 <- duos(y = y, k = 4, N = 20000, start = runif(4, 1/3, 2/3))
//' duos_unif3 <- duos(y = y, k = 4, N = 20000, start = runif(4, 2/3, 1))
//' 
//' # Load package coda
//' library(coda)
//' 
//' # Turn matrices of cut-points into mcmc objects
//' C1 <- mcmc(duos_unif1$C)
//' C2 <- mcmc(duos_unif2$C)
//' C3 <- mcmc(duos_unif3$C)
//' 
//' # Turn into mcmc list
//' all_chains <- mcmc.list(C1, C2, C3)
//' 
//' #Get Gelman and Rubin diagnostic
//' gelman.diag(all_chains)[[1]][,1]
//' 
//' #Examine estimate of PDF's from each run
//' duos_plot(duos_unif1, type="pdf", data=TRUE)
//' duos_plot(duos_unif2, type="pdf", data=TRUE)
//' duos_plot(duos_unif3, type="pdf", data=TRUE)


// [[Rcpp::export]]
List duos(NumericVector y, double k=NA_REAL, double N=20000, NumericVector alpha=1, double scale_l=0.00001, double scale_u=0.00001, NumericVector start = NA_REAL){
  
  // Set default for cut-points if a values was not entered ********************************************************************//
  
  if(k!=k){
    k = round(y.size()/50)+3;
  }
  
  NumericVector alpha_output(k+1); 
  
  
  // All error checking ********************************************************************//
  
  
  //check if y contains any missing values
  if(std::find(y.begin(), y.end(), NA)!=y.end()){
    throw std::range_error("The input vector 'y' contains missing values. Please remove the missing values.");
  }
  
  
  // Check starting values to make sure values are within the appropriate range and length
  if(start[0]==start[0]){
    if(k!=start.size()){
      throw std::range_error("The length of the vector in start should be the number of cut-points: k.");
    }
    double max_start = *std::max_element(start.begin(), start.end());
    double min_start = *std::min_element(start.begin(), start.end());
    
    
    if(max_start>=1){
      throw std::range_error("The maximum start value is outside the appropriate range: (0, 1).");
    }
    
    if(min_start>=1){
      throw std::range_error("The minimum start value is outside the appropriate range: (0, 1).");
    }
    
  }
  
  // Scale paramters needs to be greater than or equal to 0
  if (scale_l < 0.0) {
    throw std::range_error("Scale parameter should be >= 0.");
  }
  if (scale_u < 0.0) {         	// scale paramter needs to be greater than or equal to 0
    throw std::range_error("Scale parameter should be >= 0.");
  }
  
  // Make sure k is an integer since allowed it to be double otherwise duos does not throw an error
  
  if(k/ceil(k)<1){
    throw std::range_error("k should be an integer value greater than zero.");
  }
  
  // Throw an error is k is less than 1
  if(k<=0){
    throw std::range_error("k should be an integer value greater than zero.");
  }
  
  // Make sure N is an integer 
  if(N/ceil(N)<1){
    throw std::range_error("N should be an integer value.");
  }
  
  // Make sure N is greater than 0
  if(N<=0){
    throw std::range_error("N should be an integer greater than zero.");
  }
  
  // Set default values for alpha
  if((alpha[0]==1)&(alpha.size()==1)){
    for(int j=0; j<k+1; j++){
      alpha_output[j] = 1;
    }   
  }else{
    for(int j=0; j<k+1; j++){
      alpha_output[j] = alpha[j];
    }
  }
  
  // Alpha has to be greater than zero
  for(int j=0; j<alpha_output.size(); j++){
    if(alpha_output[j]<=0){
      throw std::range_error("alpha should be greater then 0.");
    }
  }
  
  //Throw error if alpha not the right length
  if(alpha_output.size()!=(k+1)){
    throw std::range_error("alpha does not contain the correct number of parameters. It should have the length: k+1.");
  }
  
  
  // **************************************************************************************//
  
  // Begin set up for algorithm *********************************************************//
  
  // Find max and min of data
  
  double max_y = *std::max_element(y.begin(), y.end());
  double min_y = *std::min_element(y.begin(), y.end());
  
  // Create a new vector to analyze: contains either same data as y or scaled data
  NumericVector y_analyize(y.size());
  
  // Check if data is outside (0,1) range in order to scale if necessary
  if((max_y>1)|(min_y<0)){
    // Scale the data
    for(int j=0; j<y.size(); j++){
      y_analyize[j]=(y[j]-(min_y-scale_l))/(max_y+scale_u-(min_y-scale_l));
    }}else{
      // Do not scale
      for(int j=0; j<y.size(); j++){
        y_analyize[j]=y[j];
      }
    }
    
    
    // The data needs to be sorted for the algorithm
    std::sort(y_analyize.begin(), y_analyize.end());
    
    // Declare 2D arrays to store cut-points and bin proportion parameters
    NumericMatrix C(N,k);
    NumericMatrix P(N,k+1);
    // temp
    
    
    
    // Initial values for cutpoints and bin proportions
    NumericVector c_init(k);
    NumericVector p_init(k+1);
    // The full version is created to contain 0 and 1 in addition to the ordered cut-points
    NumericVector c_init_full(k+2);
    
    // Sum of gamma simulations to normalize to sum to one
    double sum_p=0;
    
    // Check if have starting values
    if(start[0]!=start[0]){
      //Initialize values
      for(int i = 0; i < k; i++) {
        c_init[i] = runif(1, 0, 1)[0];
      }
      //Sort from small to large
      std::sort(c_init.begin(), c_init.end());
    }else{
      // Initialize values to the provided start values
      for(int i = 0; i < k; i++) {
        c_init[i] = start[i];
      }
      // Sort from small to large
      std::sort(c_init.begin(), c_init.end());
      
    }
    
    // Assign c_init to first row of matrix of cut-point interations
    for(int i = 0; i < k; i++) {
      C(0,i) = c_init[i];
    }
    
    // Create c_init_full: with 0 and 1 added to the ends
    c_init_full[0] = 0;
    for(int i=1; i<k+1; i++){
      c_init_full[i]=c_init[i-1];
    }
    c_init_full[k+1] = 1;
    
    
    // Create initial p vector by simulating from a gamma distribution
    for(int i=0; i<k+1; i++){
      p_init[i] = rgamma(1, alpha_output[i]+sum((y_analyize>=c_init_full[i]) & (y_analyize<c_init_full[i+1])),1)[0];
    }
    
    // Normalize initial p
    for(int i=0; i<k+1; i++){
      sum_p+=p_init[i];
    }
    // Add initial bin proportions to the matrix storing their iterations
    for(int i=0; i<k+1; i++){
      P(0,i) = p_init[i]/sum_p;
      
    }
    
    // Gibbs algorithm
    for(int i=1; i<N;i++){
      
      // Vectors to contain previous runs of parameters
      NumericVector Q_C(k);
      NumericVector Q_P(k+1);
      NumericVector Q_C_extra(k+2);
      
      // Index for cutpoints
      IntegerVector k_index(k);
      
      // Start by initializing vectors with previous run
      Q_C = C(i-1,_);
      Q_P = P(i-1,_);
      
      // Create vector to store previous run of cut-points with 0 and 1 added to the ends
      Q_C_extra[0]=0;
      for(int j=1; j<k+1; j++){
        Q_C_extra[j]=Q_C[j-1];
      }
      Q_C_extra[k+1]=1;
      
      // Randomize order in which to sample cut-points
      for(int j=0;j<k;j++){
        k_index[j]=j;
      }
      std::random_shuffle(k_index.begin(), k_index.end());
      
      // Loop through the cut-points in the random order provided
      for(int j=0; j<k; j++){
        // The cut-point to sample
        int w=k_index[j];
        
        // Count the number of data points between the cut-point above and below the current cut-point 
        // being sampled
        double m = sum(((y_analyize>=Q_C_extra[w]) & (y_analyize<Q_C_extra[w+2])));
        
        // Check if there are any data points between the two cut-points conditioning on
        if(m>0){
          // Vector for data between two cut-points
          NumericVector xm(m+2);
          int btw_index=1;
          
          // Get data between two bordering cut-points
          for(int l=0; l<y_analyize.size(); l++){
            if((y_analyize[l]>=Q_C_extra[w]) & (y_analyize[l]<Q_C_extra[w+2])){
              xm[btw_index]=y_analyize[l];
              btw_index+=1;
            }
          }
          
          // Add the bordering cut-points to the vector
          xm[0]=Q_C_extra[w];
          xm[xm.size()-1]=Q_C_extra[w+2];
          
          // Initialize vector to contains maximums on each interval
          NumericVector c_max(m+1);
          // Initialize vector to contain choices
          NumericVector c_max_choices1(2);
          
          // Initialize value that will contain ratio to find maximum
          double c_max_maxes1;
          
          // Get choices for first interval in the vector
          c_max_choices1[0]=xm[0];
          c_max_choices1[1]=xm[1];
          
          
          // The maximum from solving the derivative is Q_C_extra[w] which is xm[1]
          
          // option1/option2
          c_max_maxes1= std::pow((Q_C_extra[w+2]-c_max_choices1[0]),-m)/(std::pow((Q_C_extra[w+2]-c_max_choices1[1]),(-m)));
          
          // Find maximum by comparing ratio to 1
          if (c_max_maxes1>1){
            c_max[0] = c_max_choices1[0];
          } else{
            c_max[0] = c_max_choices1[1];
            
          }
          
          // As long as there is more than one data value between the bording cut-points, loop through
          // them
          if(m>1){
            for(int q=1; q<m; q++){
              // Set up vector to contain choices for the maximum
              NumericVector c_max_choices2(3);
              // Set up vector to contain ratios to find the maximum
              NumericVector c_max_maxes2(3);
              
              // First and second choice are the borders: third choice is from solving the derivative
              c_max_choices2[0]=xm[q];
              c_max_choices2[1]=xm[q+1];
              c_max_choices2[2]=(((m-q)*Q_C_extra[w]+(q)*Q_C_extra[w+2])/m);
              
              // option1/option2
              c_max_maxes2[0] = exp((-q)*log(c_max_choices2[0]-Q_C_extra[w])+(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[0])-(-q)*log(c_max_choices2[1]-Q_C_extra[w])-(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[1]));
              // Option1/option3
              c_max_maxes2[1] = exp((-q)*log(c_max_choices2[0]-Q_C_extra[w])+(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[0])-(-q)*log(c_max_choices2[2]-Q_C_extra[w])-(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[2]));
              // Option2/option3
              c_max_maxes2[2] = exp((-q)*log(c_max_choices2[1]-Q_C_extra[w])+(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[1])-(-q)*log(c_max_choices2[2]-Q_C_extra[w])-(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[2]));
              
              // Use ratios to find maximum
              if ((c_max_choices2[2] >= xm[q]) & (c_max_choices2[2] < xm[q+1])){
                
                if ((c_max_maxes2[0] > 1) & (c_max_maxes2[1] > 1)) {
                  
                  c_max[q] = c_max_choices2[0];
                } else if ((c_max_maxes2[0] < 1) & (c_max_maxes2[2] > 1)) {
                  
                  c_max[q] = c_max_choices2[1];
                } else {
                  c_max[q] = c_max_choices2[2];
                }
              }else if ( c_max_maxes2[0] > 1) {
                
                c_max[q] = c_max_choices2[0];
              } else {
                c_max[q] = c_max_choices2[1];
              }
            } // loop through m-1 middle intervals looking for maximums
            
          } // if m>1
          
          // Find maximum in last bin
          // Set up vector for choices
          NumericVector c_max_choices3(2);
          // Set up value to contain ratio
          double c_max_maxes3;
          // Set up the choices
          c_max_choices3[0]=xm[m];
          c_max_choices3[1] = xm[m+1];
          
          // option1/option2
          c_max_maxes3 = std::pow((c_max_choices3[0]-Q_C_extra[w]),-m)/(std::pow((c_max_choices3[1]-Q_C_extra[w]),-m));
          
          // Compare the ratio to one to find the maximum
          if (c_max_maxes3 > 1){
            
            c_max[m] = c_max_choices3[0];
          } else {
            c_max[m] = c_max_choices3[1];
          }
          
          // Set up for calculating the normalizing constant
          // Alternatives are calculating these values on the log and a multiple scale 
          NumericVector nk_indiv_log(m+1);
          NumericVector nk_indiv(m+1);
          NumericVector nk_indiv_mult(m+1);
          
          // Find the inegral over each bin 
          for (int v=0; v<m+1;v++) {
            nk_indiv_log[v] = (v)*log(Q_P[w])-(v)*log(c_max[v]-Q_C_extra[w])+(m-v)*log(Q_P[w+1])-(m-v)*log(Q_C_extra[w+2]-c_max[v])+log(xm[v+1]-xm[v]);
            nk_indiv[v] = std::pow((Q_P[w]/(c_max[v]-Q_C_extra[w])),v)*std::pow((Q_P[w+1]/(Q_C_extra[w+2]-c_max[v])),m-v)*(xm[v+1]-xm[v]);
            nk_indiv_mult[v] = std::pow((10*Q_P[w]/(c_max[v]-Q_C_extra[w])),v)*std::pow((10*Q_P[w+1]/(Q_C_extra[w+2]-c_max[v])),m-v)*(xm[v+1]-xm[v]);
          }
          
          NumericVector sum_across(m);
          double sum_across_sum1=0;
          double nk_log;
          double nk=0;
          double nk_mult=0;
          
          
          // Normalizing constant for bounding function on log scale
          for(int l=0; l<m; l++){
            sum_across[l] = std::pow(exp(1),nk_indiv_log[l+1]-nk_indiv_log[0]);
          }
          for(int l=0;l<sum_across.size(); l++){
            sum_across_sum1+=sum_across[l];
          }
          nk_log = nk_indiv_log[0]+log(1+sum_across_sum1) ;
          
          // Normalizing constant for bounding function on true scale
          for(int l=0; l<nk_indiv.size(); l++){
            nk+=nk_indiv[l];
          }
          
          // Normalizing constant for bounding function on multiple scale
          for(int l=0; l<nk_indiv.size(); l++){
            nk_mult+=nk_indiv_mult[l];
          }
          
          //Find probabilites of each interval (these sum to one)
          std::vector< std::pair <double,int> > discrete_prob_log(m+1);
          std::vector< std::pair <double,int> > discrete_prob(m+1);
          std::vector< std::pair <double,int> > discrete_prob_mult(m+1);
          
          // Find probabilities to use for calculating the cumulative probability
          // Give each value in vector a number indicating with interval
          for(int l=0; l<m+1; l++){
            discrete_prob_log[l].first=nk_indiv_log[l]-nk_log;
            discrete_prob_log[l].second=l;
            discrete_prob[l].first=nk_indiv[l]/nk;
            discrete_prob[l].second=l;
            discrete_prob_mult[l].first=nk_indiv_mult[l]/nk_mult;
            discrete_prob_mult[l].second=l;
          }
          
          // Sort in desceding order
          sort(discrete_prob_log.begin(), discrete_prob_log.end());
          reverse(discrete_prob_log.begin(), discrete_prob_log.end());
          
          std::sort(nk_indiv_log.begin(), nk_indiv_log.end());
          std::reverse(nk_indiv_log.begin(), nk_indiv_log.end());
          
          sort(discrete_prob.begin(), discrete_prob.end());
          reverse(discrete_prob.begin(), discrete_prob.end());
          
          sort(discrete_prob_mult.begin(), discrete_prob_mult.end());
          reverse(discrete_prob_mult.begin(), discrete_prob_mult.end());
          
          // Calculate cumulative probabilites
          
          std::vector< std::pair <double,int> > cum_prob_log(m+1);
          std::vector< std::pair <double,int> > cum_prob_mult(m+1);
          std::vector< std::pair <double,int> > cum_prob(m+1);
          
          // Set up the first value
          cum_prob[0].first=discrete_prob[0].first;
          cum_prob[0].second=discrete_prob[0].second;
          cum_prob_log[0].first=discrete_prob_log[0].first;
          cum_prob_log[0].second=discrete_prob_log[0].second;
          cum_prob_mult[0].first=discrete_prob_mult[0].first;
          cum_prob_mult[0].second=discrete_prob_mult[0].second;
          
          double sum_discrete_prob=0;
          double sum_across_sum=0;
          
          // Total sum of probabilites
          for(int l=0; l<m+1; l++){
            sum_discrete_prob+=discrete_prob[l].first;
          }
          
          // Calculate cumulative probabilities
          for (int l=1; l<m+1; l++){
            if(sum_discrete_prob!=0){
              cum_prob[l].first=cum_prob[l-1].first+discrete_prob[l].first;
              cum_prob[l].second=discrete_prob[l].second;
            }
            cum_prob_mult[l].first=cum_prob_mult[l-1].first+discrete_prob_mult[l].first;
            cum_prob_mult[l].second=discrete_prob_mult[l].second;
            
            for(int q=1; q<l+1; q++){
              sum_across_sum += std::pow(exp(1),nk_indiv_log[q]-nk_indiv_log[0]);
            }
            cum_prob_log[l].first=nk_indiv_log[0]+log(1+sum_across_sum)-nk_log;
            cum_prob_log[l].second=discrete_prob_log[l].second;
            sum_across_sum=0;
            
          }
          
          // Create an indicator that tells whether or not rejection sampling found a value
          int found=0;
          
          // Rejeciton Sampling
          while (found==0){
            // All of the following is to sample bin from discrete uniform
            double unif_rand;
            double unif_rand_log;
            double ci_new;
            double new_prob;
            double unif_prop;
            int interval;
            
            // Check if all probabilities too small to handle
            if (sum_discrete_prob!=0){
              unif_rand = runif(1,0,1)[0];
              int select =0;
              while(cum_prob[select].first<unif_rand){
                select+=1;
              }
              interval= cum_prob[select].second;
              
            }else if ((nk_log!=nk_log) | (std::isinf(nk_log))){
              // Sample on multiple scale to make values larger
              unif_rand = runif(1,0,1)[0];
              int select =0;
              while(cum_prob_mult[select].first<unif_rand){
                select+=1;
              }
              interval= cum_prob_mult[select].second;
              
            }else{
              // Sample on log scale
              unif_rand_log = log(runif(1,0,1)[0]);
              
              int select =0;
              while(cum_prob_log[select].first<unif_rand_log){
                select+=1;
              }
              
              interval = cum_prob_log[select].second;
              
            }
            
            // Propose a new value usin the selected value
            ci_new = runif(1, xm[interval], xm[interval+1])[0];
            
            // Compute ratio of bounding and original function
            new_prob = (-interval)*log(ci_new-Q_C_extra[w])+(-m+interval)*log(Q_C_extra[w+2]-ci_new)+(interval)*log(c_max[interval]-Q_C_extra[w])+(m-interval)*log(Q_C_extra[w+2]-c_max[interval]);
            unif_prop = log(runif(1,0,1)[0]);
            
            // Check if has large enough probability
            if (unif_prop<new_prob){
              found=1;
              Q_C[w]=ci_new;
              Q_C_extra[w+1]=ci_new;
            }
            
          } // while loop
          
        }else{ // If there is no data between the two cut-points, distribution is just uniform
          // over the bounding cut-points
          Q_C[w] = runif(1, Q_C_extra[w], Q_C_extra[w+2])[0];
          Q_C_extra[w+1]=Q_C[w];
        }// if m>0 condition
        
        // Assign the new
        C(i,w)=Q_C[w];
        
        // Create c_init_full
        c_init_full[0] = 0;
        for(int l=1; l<k+1; l++){
          c_init_full[l]=Q_C[l-1];
        }
        c_init_full[k+1] = 1;
        
        // Create initial p vector
        for(int l=0; l<k+1; l++){
          p_init[l] = rgamma(1, alpha_output[l]+sum((y_analyize>=c_init_full[l]) & (y_analyize<c_init_full[l+1])),1)[0];
        }
        
        sum_p =0;
        // Normalize initial p
        for(int l=0; l<k+1; l++){
          sum_p+=p_init[l];
        }
        for(int l=0; l<k+1; l++){
          Q_P[l] = p_init[l]/sum_p;
        }
        
        
      } // loop through cut points
      
      for(int j=0; j<k+1; j++){
        P(i,j) = Q_P[j];
      }
      
      // Output iteration    
      if(i % 1000 == 0){
        Rcout  << "Iteration: " << i << std::endl;
      }
      
    } // Gibbs iterations
    
    
    
    // Return list containing C and P matrix
    return List::create(
      _["C"] = C,
      _["P"] = P,
      _["y"] = y,
      _["k"] = k,
      _["alpha"] = alpha_output,
      _["scale_l"] = scale_l,
      _["scale_u"] = scale_u);
}
