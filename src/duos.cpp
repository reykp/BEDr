// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List duos(NumericVector ys,int MH_N=20000,int k=12, double alpha=1){
  //Declare list to return: will contain two matrices
  List C_P(2);
  //sort data
  std::sort(ys.begin(), ys.end());

  //Declare 2D arrays to store parameters
  NumericMatrix C(MH_N,k);
  NumericMatrix P(MH_N,k+1);
  NumericVector c_init(k);
  NumericVector p_init(k+1);
  NumericVector c_init_full(k+2);
  double sum_p=0;

  //Initialize values
  for(int i = 0; i < k; i++) {
    c_init[i] = runif(1, 0, 1)[0];
  }
  //Sort from small to large
  std::sort (c_init.begin(), c_init.end());

  //assign c_init to first row of matrix
  for(int i = 0; i < k; i++) {
    C(0,i) = c_init[i];
  }

  //Create c_init_full
  c_init_full[0] = 0;
  for(int i=1; i<k+1; i++){
    c_init_full[i]=c_init[i-1];
  }
  c_init_full[k+1] = 1;


  //Create initial p vector
  for(int i=0; i<k+1; i++){
    p_init[i] = rgamma(1, alpha+sum((ys>=c_init_full[i]) & (ys<c_init_full[i+1])),1)[0];

  }

  //Normalize initial p
  for(int i=0; i<k+1; i++){
    sum_p+=p_init[i];
  }
  for(int i=0; i<k+1; i++){
    P(0,i) = p_init[i]/sum_p;

  }

  //Gibbs algorithm
  for(int i=1; i<MH_N;i++){

    NumericVector Q_C(k);
    NumericVector Q_P(k+1);
    NumericVector Q_C_extra(k+2);
    IntegerVector k_index(k);

    //Start by initializing vectors with previous run to
    Q_C = C(i-1,_);
    Q_P = P(i-1,_);


    //Calculate vector of c's to use for Gibbs with c0<-0 and ck+1<-1
    Q_C_extra[0]=0;
    for(int j=1; j<k+1; j++){
      Q_C_extra[j]=Q_C[j-1];
    }
    Q_C_extra[k+1]=1;


    //Randomize order in which to sample cut-points
    for(int j=0;j<k;j++){
      k_index[j]=j;
    }
    std::random_shuffle(k_index.begin(), k_index.end());

    for(int j=0; j<k; j++){
      //cut-point to sample
      int w=k_index[j];

      double m = sum((ys>=Q_C_extra[w]) & (ys<Q_C_extra[w+2]));

      if(m>0){
        //Vector for data between two cut-points
        NumericVector xm(m+2);
        int btw_index=1;

        //Get data between two end points
        for(int l=0; l<ys.size(); l++){
          if((ys[l]>=Q_C_extra[w]) & (ys[l]<Q_C_extra[w+2])){
            xm[btw_index]=ys[l];
            btw_index+=1;
          }
        }
        //Add C[w-1] and C[w+1] to vector
        xm[0]=Q_C_extra[w];
        xm[xm.size()-1]=Q_C_extra[w+2];

        //Initialize vector to contains maximums on each interval
        NumericVector c_max(m+1);
        //Initialize vector to contain minimums on each interval
        NumericVector c_max_choices1(2);

        double c_max_maxes1;

        //Get for first interval in the vector
        c_max_choices1[0]=xm[0];
        c_max_choices1[1]=xm[1];


        //The maximum from derivatives is Q_C_extra[w] which is xm[1]
        //Pairs (c(1/2), c(1/3), c(2/3))
        //Create to contain ratios of the pdf's

        //option1/option2
        c_max_maxes1= pow((Q_C_extra[w+2]-c_max_choices1[0]),-m)/(pow((Q_C_extra[w+2]-c_max_choices1[1]),(-m)));


        if (c_max_maxes1>1){
          c_max[0] = c_max_choices1[0];

        } else{
          c_max[0] = c_max_choices1[1];

        }

        if(m>1){

          for(int q=1; q<m; q++){

            NumericVector c_max_choices2(3);
            NumericVector c_max_maxes2(3);

            c_max_choices2[0]=xm[q];
            c_max_choices2[1]=xm[q+1];
            c_max_choices2[2]=(((m-q)*Q_C_extra[w]+(q)*Q_C_extra[w+2])/m);

            //option1/option2
            c_max_maxes2[0] = exp((-q)*log(c_max_choices2[0]-Q_C_extra[w])+(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[0])-(-q)*log(c_max_choices2[1]-Q_C_extra[w])-(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[1]));
            //Option1/option3
            c_max_maxes2[1] = exp((-q)*log(c_max_choices2[0]-Q_C_extra[w])+(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[0])-(-q)*log(c_max_choices2[2]-Q_C_extra[w])-(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[2]));
            //Option2/option3
            c_max_maxes2[2] = exp((-q)*log(c_max_choices2[1]-Q_C_extra[w])+(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[1])-(-q)*log(c_max_choices2[2]-Q_C_extra[w])-(-m+q)*log(Q_C_extra[w+2]-c_max_choices2[2]));


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
          } //loop through m-1 middle intervals looking for maximums

        } //if m>1

        NumericVector c_max_choices3(2);
        double c_max_maxes3;
        c_max_choices3[0]=xm[m];
        c_max_choices3[1] = xm[m+1];

        //Pairs (c(1/2), c(1/3), c(2/3))
        //option1/option2
        c_max_maxes3 = pow((c_max_choices3[0]-Q_C_extra[w]),-m)/(pow((c_max_choices3[1]-Q_C_extra[w]),-m));

        if (c_max_maxes3 > 1){

          c_max[m] = c_max_choices3[0];
        } else {
          c_max[m] = c_max_choices3[1];
        }


        //normalizing constant
        NumericVector nk_indiv_log(m+1);
        NumericVector nk_indiv(m+1);
        NumericVector nk_indiv_mult(m+1);

        for (int v=0; v<m+1;v++) {
          nk_indiv_log[v] = (v)*log(Q_P[w])-(v)*log(c_max[v]-Q_C_extra[w])+(m-v)*log(Q_P[w+1])-(m-v)*log(Q_C_extra[w+2]-c_max[v])+log(xm[v+1]-xm[v]);

          nk_indiv[v] = pow((Q_P[w]/(c_max[v]-Q_C_extra[w])),v)*pow((Q_P[w+1]/(Q_C_extra[w+2]-c_max[v])),m-v)*(xm[v+1]-xm[v]);
          //nk_indiv[v] = ys.size()*(pow((Q_P[w]/(c_max[v]-Q_C_extra[w])),v/ys.size())*pow((Q_P[w+1]/(Q_C_extra[w+2]-c_max[v])),(m-v)/ys.size()))*(xm[v+1]-xm[v]);

          nk_indiv_mult[v] = pow((10*Q_P[w]/(c_max[v]-Q_C_extra[w])),v)*pow((10*Q_P[w+1]/(Q_C_extra[w+2]-c_max[v])),m-v)*(xm[v+1]-xm[v]);

        }

        //Normalizing constant for bounding function
        NumericVector sum_across(m);
        double sum_across_sum1=0;
        double nk_log;
        double nk=0;
        double nk_mult=0;

        for(int l=0; l<m; l++){
          sum_across[l] = pow(exp(1),nk_indiv_log[l+1]-nk_indiv_log[0]);

        }
        for(int l=0;l<sum_across.size(); l++){
          sum_across_sum1+=sum_across[l];
        }
        nk_log = nk_indiv_log[0]+log(1+sum_across_sum1) ;

        //Normalizing constant for bounding function
        for(int l=0; l<nk_indiv.size(); l++){
          nk+=nk_indiv[l];
        }
        //nk = sum(nk_indiv);

        for(int l=0; l<nk_indiv.size(); l++){
          nk_mult+=nk_indiv_mult[l];
        }
        //nk_mult = sum(nk_indiv_mult);


        //Find probabilites of each interval (these sum to one)
        std::vector< std::pair <double,int> > discrete_prob_log(m+1);
        std::vector< std::pair <double,int> > discrete_prob(m+1);
        std::vector< std::pair <double,int> > discrete_prob_mult(m+1);

        for(int l=0; l<m+1; l++){
          //discrete_prob_log.push_back(std::make_pair((nk_indiv_log[l]-nk_log),l) );
          discrete_prob_log[l].first=nk_indiv_log[l]-nk_log;
          discrete_prob_log[l].second=l;
          discrete_prob.push_back(std::make_pair(nk_indiv[l]/nk,l) );
          discrete_prob_mult.push_back(std::make_pair(nk_indiv_mult[l]/nk_mult,l) );
        }


        //Give each value in vector a number indicating with interval
        //discrete_prob_log.names() = CharacterVector::create(vector_names);


        //Sort in desceding order
        sort(discrete_prob_log.begin(), discrete_prob_log.end());
        reverse(discrete_prob_log.begin(), discrete_prob_log.end());

        std::sort(nk_indiv_log.begin(), nk_indiv_log.end());
        std::reverse(nk_indiv_log.begin(), nk_indiv_log.end());

        sort(discrete_prob.begin(), discrete_prob.end());
        reverse(discrete_prob.begin(), discrete_prob.end());

        sort(discrete_prob_mult.begin(), discrete_prob_mult.end());
        reverse(discrete_prob_mult.begin(), discrete_prob_mult.end());

        //Calculat cumulative probabilites

        std::vector< std::pair <double,int> > cum_prob_log(m+1);
        std::vector< std::pair <double,int> > cum_prob_mult(m+1);
        std::vector< std::pair <double,int> > cum_prob(m+1);


        cum_prob[0].first=discrete_prob[0].first;
        cum_prob[0].second=discrete_prob[0].second;
        cum_prob_log[0].first=discrete_prob_log[0].first;
        cum_prob_log[0].second=discrete_prob_log[0].second;
        cum_prob_mult[0].first=discrete_prob_mult[0].first;
        cum_prob_mult[0].second=discrete_prob_mult[0].second;



        double sum_discrete_prob=0;
        double sum_across_sum=0;

        for(int l=0; l<m+1; l++){
          sum_discrete_prob+=discrete_prob[l].first;

        }
        for (int l=1; l<m+1; l++){

          if(sum_discrete_prob!=0){
            cum_prob[l].first=cum_prob[l-1].first+discrete_prob[l].first;
            cum_prob[l].second=discrete_prob[l].second;
          }

          cum_prob_mult[l].first=cum_prob_mult[l-1].first+discrete_prob_mult[l].first;
          cum_prob_mult[l].second=discrete_prob_mult[l].second;

          for(int q=1; q<l+1; q++){
            sum_across_sum += pow(exp(1),nk_indiv_log[q]-nk_indiv_log[0]);
          }

          cum_prob_log[l].first=nk_indiv_log[0]+log(1+sum_across_sum)-nk_log;
          cum_prob_log[l].second=discrete_prob_log[l].second;
          sum_across_sum=0;

        }


        int found=0;


        while (found==0){
          double unif_rand;
          double unif_rand_log;
          double ci_new;
          double new_prob;
          double unif_prop;
          int interval;

          if (sum_discrete_prob!=0){
            unif_rand = runif(1,0,1)[0];
            int select =0;
            while(cum_prob[select].first<unif_rand){
              select+=1;
            }

            interval= cum_prob[select].second;

          }else if ((nk_log!=nk_log) | (std::isinf(nk_log))){
            unif_rand = runif(1,0,1)[0];
            int select =0;
            while(cum_prob_mult[select].first<unif_rand){
              select+=1;
            }

            interval= cum_prob_mult[select].second;

          }else{
            unif_rand_log = log(runif(1,0,1)[0]);
            int select =0;
            while(cum_prob_log[select].first<unif_rand_log){
              select+=1;
            }

            interval = cum_prob_log[select].second;

          }

          ci_new = runif(1, xm[interval], xm[interval+1])[0];


          new_prob = (-interval)*log(ci_new-Q_C_extra[w])+(-m+interval)*log(Q_C_extra[w+2]-ci_new)+(interval)*log(c_max[interval]-Q_C_extra[w])+(m-interval)*log(Q_C_extra[w+2]-c_max[interval]);
          unif_prop = log(runif(1,0,1)[0]);

          //if ((unif_prop<new_prob)&(ci_new>Q_C_extra[w])&(ci_new<Q_C_extra[w+2])){
          if (unif_prop<new_prob){
            found=1;
            Q_C[w]=ci_new;
            Q_C_extra[w+1]=ci_new;
          }

        } //while loop





      }else{
        Q_C[w] = runif(1, Q_C_extra[w], Q_C_extra[w+2])[0];
        Q_C_extra[w+1]=Q_C[w];
      }// if m>0 condition

      //assign new value
      C(i,w)=Q_C[w];

      //Create c_init_full
      c_init_full[0] = 0;
      for(int l=1; l<k+1; l++){
        c_init_full[l]=Q_C[l-1];
      }
      c_init_full[k+1] = 1;

      //Create initial p vector
      for(int l=0; l<k+1; l++){
        p_init[l] = rgamma(1, alpha+sum((ys>=c_init_full[l]) & (ys<c_init_full[l+1])),1)[0];
      }

      sum_p =0;
      //Normalize initial p
      for(int l=0; l<k+1; l++){
        sum_p+=p_init[l];
      }
      for(int l=0; l<k+1; l++){
        Q_P[l] = p_init[l]/sum_p;
      }


    } //loop through cut points

    for(int j=0; j<k+1; j++){
      P(i,j) = Q_P[j];
    }


  } //Gibbs iterations

  //Assign C and P to list
  C_P[0] =C;
  C_P[1] = P;

  //Return list containing C and P matrix
  return C_P;
}
