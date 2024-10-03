#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]



// [[Rcpp::export]]
arma::mat rMVNormCpp(int n, arma::vec Mean, arma::mat Var) {
  int ncols = Var.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(Mean, 1, n).t() + Y * arma::chol(Var);
}



//[[Rcpp::export]]
arma::mat sim_thetacpp(int S, arma::vec lambda, int n_sources, 
                           int n_tracers, bool solo){
  
  arma::mat theta(S, (n_sources + n_tracers));
  
  arma::vec mean = lambda.subvec(0, n_sources-1);
  
  // for(int i = 0; i<n_sources; i++){
  //   mean(i) = lambda(i);
  // }
  
  
  arma::mat chol_prec(n_sources, n_sources);
  
  
  int count = 0;
  for(int i = 0; i< n_sources; i++){ 
    for(int j = 0; j<n_sources; j++){
      if (i <= j){
        count +=1;
        chol_prec((i),(j)) = lambda(n_sources -1 +count);
        
        
      }
      else{
        chol_prec(i,j) = 0;
      }
      
    }
  }
  
  arma::mat normmat(S, n_sources);
  
  arma::mat prec(n_sources, n_sources);
  prec = chol_prec.t() * chol_prec;

  
  arma::mat solve_prec(n_sources, n_sources);
  
  arma::mat b = arma::eye(prec.n_rows, prec.n_cols);
  
  solve_prec = solve(prec, b);
  
  // solvearma returns the transpose of the inverse of the precision matrix
  //so need to transpose it back
  arma::mat var(n_sources, n_sources);
  for(int i=0; i<n_sources; i++){
    for(int j=0; j<n_sources; j++){
      var(i,j) = solve_prec(j,i);
    }
  }
  
  
  normmat = rMVNormCpp(S, mean, var);
  
  
  arma::vec solovec(S);
  
  for(int s = 0; s<S; s++){
    solovec(s) = 1000;
  }
  

  
  int lambda_offset = n_sources + (n_sources * (n_sources + 1)) / 2;
  
  if (!solo) {
    // When solo == FALSE
    for (int i = 0; i < n_tracers; i++) {
      // Generate gamma random numbers and convert to arma::vec
      Rcpp::NumericVector gamma_vals = Rcpp::rgamma(S, lambda(lambda_offset + i), 1 / lambda(lambda_offset + n_tracers + i));
      theta.col(i + n_sources) = arma::vec(gamma_vals.begin(), gamma_vals.size(), false);
    }
  } else {
    // When solo == TRUE
    for (int i = 0; i < n_tracers; i++) {
      // Generate gamma random numbers and convert to arma::vec
      Rcpp::NumericVector gamma_vals = Rcpp::rgamma(S, lambda(lambda_offset + i), 1 / lambda(lambda_offset + n_tracers + i));
      arma::vec gamma_vec(gamma_vals.begin(), gamma_vals.size(), false);
      theta.col(i + n_sources) = solovec + 0.00001 * gamma_vec;
    }
  }
  
  
  
  for(int i=0; i<n_sources; i++){
    theta.col(i) = normmat.col(i);
  }
  
  
  
  return theta;
  
}



// This function takes theta and calculates the proportions
//[[Rcpp::export]]
arma::vec hfn(arma::vec theta, int n_sources) {
  // Compute exp of each theta and store in exptheta
  arma::vec theta_sub = theta.head(n_sources);
  arma::vec exptheta = arma::exp(theta_sub);
  
  // Compute sum of all exptheta
  double sumexptheta = arma::sum(exptheta);
  
  // Compute p as element-wise division of exptheta by sumexptheta
  arma::vec p = exptheta / sumexptheta;
  
  return p;
}


//[[Rcpp::export]]
double hcpp(int n_sources, int n_isotopes,
            double beta_prior,
            arma::mat concentrationmeans, arma::mat sourcemeans,
            arma::mat correctionmeans,
            arma::mat corrsds, arma::mat sourcesds, 
            arma::vec theta, arma::mat y ){
  
  double x = 0;
  
  arma::vec p(n_sources);
  
  p = hfn(theta, n_sources);
  
  double ly = y.n_rows;
  
  // Set prior values for hyperparameters using vectorization
  arma::vec prior_means = arma::zeros<arma::vec>(n_sources);  // Vector of zeros
  arma::vec prior_sd = arma::ones<arma::vec>(n_sources);      // Vector of ones
  
  // Set the prior values for c_0 and d_0
  arma::vec c_0 = arma::ones<arma::vec>(n_isotopes);          // Vector of ones
  arma::vec d_0 = arma::ones<arma::vec>(n_isotopes) * beta_prior;  // All values set to beta_prior
  
  
  if(n_isotopes == 2){
    
    
    // This is to get dnorm(y[,1], sum(p*q) etc)
    double mutop1 = 0;
    double mubtm1 = 0;
    double mu1 = 0;
    double sigmasq1 = 0;
    double sigmatopsq1 = 0;
    double sigmabtmsq1 = 0;
    
    
    
    // Calculate numerator and denominator of mu
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_isotopes; j++){
        mutop1 +=  p(i)*concentrationmeans(i,j) * (sourcemeans(i,0) + correctionmeans(i,0));
        mubtm1 += p(i) * concentrationmeans(i,j);
      }
    }
    
    // Same for sigma
    for(int i=0; i<n_sources; i++){
      for(int j =0; j<n_isotopes; j++){
        sigmatopsq1 += pow(p(i),2) * pow(concentrationmeans(i,j),2) * (pow(sourcesds(i,0),2) +
          pow(corrsds(i,0),2));
        sigmabtmsq1 += pow(p(i),2) * pow(concentrationmeans(i,j),2);
      }
    }
    
    //Calculate mu and sd
    mu1 = mutop1/mubtm1;
    sigmasq1 = sigmatopsq1/sigmabtmsq1;
    double sigma1 = pow(sigmasq1 + 1/theta((n_sources)), 0.5);
    
    
    // This is to get dnorm(y[,2], sum(p*q) etc)
    double mutop2 = 0;
    double mubtm2 = 0;
    double mu2 = 0;
    double sigmasq2 = 0;
    double sigmatopsq2 = 0;
    double sigmabtmsq2 = 0;
    for(int i=0; i<n_sources; i++){
      for(int j =0; j<n_isotopes; j++){
        mutop2 += p(i) * concentrationmeans(i,j) * (sourcemeans(i,1) + correctionmeans(i,1));
        mubtm2 += p(i) * concentrationmeans(i,j);
      }
    }
    
    
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_isotopes; j++){
        sigmatopsq2 += pow(p(i),2) * pow(concentrationmeans(i,j),2) * (pow(sourcesds(i,1),2) +
          pow(corrsds(i,1),2));
        sigmabtmsq2 += pow(p(i),2) * pow(concentrationmeans(i,j),2);
      }
    }
    
    mu2 = mutop2/mubtm2;
    sigmasq2 = sigmatopsq2/sigmabtmsq2;
    
    double sigma2 = pow(sigmasq2 + 1/theta((1+n_sources)), 0.5);
    
    double yminusmu1 = 0;
    double yminusmu2 = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu1 += pow((y(i,0) - mu1),2);
      yminusmu2 +=  pow((y(i,1) - mu2),2);
    }
    
    
    
    // This is log(dnorm(y, p*q, p^2*q^2 etc) for y1 and y2
    
    x = - ly * log(sigma1) - 0.5 * ly * log(2 * M_PI)
      - 0.5 * yminusmu1 * 1/(pow(sigma1,2))
      - ly * log(sigma2) - 0.5 * ly * log(2 * M_PI)
      - 0.5 * yminusmu2 * 1/(pow(sigma2,2));
      
      
  }  else{
    //This is for just one isotope!!
    
    // This is to get dnorm(y[,1], sum(p*q) etc)
    double mutop = 0;
    double mubtm = 0;
    double mu = 0;
    double sigmasq = 0;
    double sigmatopsq = 0;
    double sigmabtmsq = 0;
    
    // calculate mu numerator and denominator
    for(int i=0; i<n_sources; i++){
      mutop += p(i) * concentrationmeans(i,0) * (sourcemeans(i,0) + correctionmeans(i,0));
      mubtm += p(i) * concentrationmeans(i,0);
    }
    
    
    
    
    // same for sigma
    for(int i=0; i<n_sources; i++){
      sigmatopsq += pow(p(i),2) * pow(concentrationmeans(i,0),2) * (pow(sourcesds(i,0),2) +
        pow(corrsds(i,0),2));
      sigmabtmsq += pow(p(i),2) * pow(concentrationmeans(i,0),2);
    }
    
    
    
    //Calculate mu and sd
    mu = mutop/mubtm;
    sigmasq = sigmatopsq/sigmabtmsq;
    double sigma = pow(sigmasq + 1/theta(n_sources), 0.5);
    
    
    
    
    
    double yminusmu = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu += pow((y(i,0) - mu),2);
    }
    
    
    
    
    
    // This has y
    
    x = - ly * log(sigma) - 0.5 * ly * log(2 * M_PI) - 0.5 * yminusmu* 1/(pow(sigma,2));
    
    
    
    
  }
  
  double thetanorm = 0;
  
  
  
  
  for(int i = 0; i<n_sources; i++){
    thetanorm +=  - n_sources * log(prior_sd(i)) - 0.5 * log(2 * M_PI) - (pow((theta(i) - prior_means(i)), 2)
                                                                            * 1/(2 * pow(prior_sd(i), 2)));
  }
  
  double gammaprior = 0;
  for (int i=0; i <(n_isotopes); i++){
    //theta should be log_theta?? I think
    gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +(c_0(i) - 1) *theta((i+n_sources)) -
      d_0(i) * theta((i+n_sources));
    
  }
  
  double totx = x + gammaprior + thetanorm;
  
  return totx;
  
}



//[[Rcpp::export]]
double log_q_cpp(arma::vec theta, arma::vec lambda, 
                 int n_sources, int n_tracers){
  
  // Create chol_prec as a lower triangular matrix
  arma::mat chol_prec = arma::zeros<arma::mat>(n_sources, n_sources);
  int count = 0;
  
  // Fill the upper triangle of chol_prec using lambda values
  for (int i = 0; i < n_sources; ++i) {
  for (int j = 0; j < n_sources; ++j) {
      if(i<=j){
      count += 1;
      chol_prec(i, j) = lambda(n_sources - 1 + count);
      }
    }
  }
  
  // Calculate precision matrix (prec)
  arma::mat prec = chol_prec.t() * chol_prec;
  
  arma::mat b = arma::eye(chol_prec.n_rows, chol_prec.n_cols);

  
  // Solve for inverse of prec using Armadillo's efficient solve function
  arma::mat solve_prec = solve(chol_prec, b); // inv_sympd is faster for symmetric positive definite matrices
  
  arma::mat var(n_sources, n_sources);
  
  var = solve_prec.t();
  
  
  // arma::mat y(1, n_sources);
  // 
  // arma::vec mean(n_sources);
  
  arma::rowvec yminusmean = (theta.head(n_sources) - lambda.head(n_sources)).t();
  
  arma::mat Z = yminusmean * chol_prec.t();
  arma::colvec tZ = Z.t();
  
  arma::mat ZtZmat = Z * tZ;
  
  double ZtZ = ZtZmat(0,0);
  

  double prod_sig;
  
  arma::vec sig(n_sources);
  
  
  prod_sig =  (arma::sum(arma::log(arma::diagvec(chol_prec.t()))));
  double thetanorm = 0;
  
  thetanorm = -  (n_sources/2) * log(2 * M_PI) + prod_sig - 0.5 * ZtZ;
  
  
  double gamman = 0;
  for (int i=0; i <(n_tracers); i++){
    gamman += lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) * 
      log(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i)) 
    - log(tgamma(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)))) 
    +(lambda(((n_sources + (n_sources * (n_sources + 1)) / 2) + i)) - 1) * log(theta((i+n_sources))) - 
    lambda(((n_sources + (n_sources * (n_sources + 1)) / 2)) + n_tracers + i) * theta((i+n_sources));
  }
  
  double x = thetanorm + gamman;
  
  return x;
}




// [[Rcpp::export]]
arma::vec delta_lqltcpp(arma::vec lambda, arma::vec theta,
                        double eps, int n_sources, int n_tracers) {
  double k = lambda.n_elem;
  arma::vec ans(k);
  arma::vec d = arma::zeros(k);  // Initialize as zeros

  arma::vec lambdaplusd(k);
  arma::vec lambdaminusd(k);

  for (int i = 0; i < k; i++) {

    for (int j = 0; j<k; j++){
      d(j) = 0;
    }
    d(i) = eps;


    for (int j = 0; j<k; j++){
      lambdaplusd(j) = lambda(j) + d(j);
      lambdaminusd(j) = lambda(j) - d(j);
    }

    ans(i) = (log_q_cpp(theta, lambdaplusd, n_sources, n_tracers) -
      log_q_cpp(theta, lambdaminusd, n_sources, n_tracers)) / (2 * eps);

    d(i) = 0;  // Reset d(i) back to zero
  }

  return ans;
}

// // [[Rcpp::export]]
// arma::vec delta_lqltcpp(arma::vec lambda, arma::vec theta,
//                                double eps, int n_sources, int n_tracers) {
//   // eps = 0.001;
//   double k = lambda.n_elem;
//   arma::vec ans(k);
//   arma::vec d(k);
//   arma::vec lambdaplusd(k);
//   arma::vec lambdaminusd(k);
//
//
//   for(int i = 0; i<k; i++){
//
//     for (int j = 0; j<k; j++){
//       d(j) = 0;
//     }
//     d(i) = eps;
//
//
//     for (int j = 0; j<k; j++){
//       lambdaplusd(j) = lambda(j) + d(j);
//       lambdaminusd(j) = lambda(j) - d(j);
//     }
//     ans(i) = (log_q_cpp(theta, lambdaplusd, n_sources, n_tracers) -
//       log_q_cpp(theta, lambdaminusd, n_sources, n_tracers))/(2 * eps);
//   }
//   return  ans;
// }

// [[Rcpp::export]]
double h_lambdacpp(int n_sources, int n_isotopes,
                   double beta_prior,
                   arma::mat concentrationmeans, arma::mat sourcemeans,
                   arma::mat correctionmeans,
                   arma::mat corrsds, arma::mat sourcesds,
                   arma::vec theta, arma::mat y,
                   arma::vec lambda) {

  return hcpp(n_sources, n_isotopes, beta_prior, concentrationmeans, sourcemeans, correctionmeans,
              corrsds, sourcesds, theta, y) - log_q_cpp(theta, lambda, n_sources, n_isotopes);
}


// [[Rcpp::export]]
arma::mat cov_mat_cpp(arma::mat x, arma::mat y) {
  int xcol = x.n_cols;
  int ycol = y.n_cols;
  int xrow = x.n_rows;
  int yrow = y.n_rows;

  arma::vec meanx(xcol);
  arma::vec meany(ycol);
  arma::mat covmat(xcol, ycol);


  for(int i = 0; i<xcol; i++){
    meanx(i) = mean(x.col(i));
  }
  for(int i = 0; i<ycol; i++){
    meany(i) = mean(y.col(i));
  }

  arma::mat xminusmean(xrow, xcol);
  arma::mat yminusmean(yrow, ycol);

  for(int j = 0; j<xcol; j++){
    for(int i=0; i<xrow; i++){
      xminusmean(i,j) = x(i,j) - meanx(j);
    }
  }

  for(int j = 0; j<ycol; j++){
    for(int i =0; i<yrow; i++){
      yminusmean(i,j) = y(i,j) - meany(j);
    }
  }

  arma::mat sumxy(xcol, ycol);

  // arma::vec xcol(x.ncol());
  // arma::vec ycol(y.ncol());

  for(int i = 0; i<xcol; i++){
    for(int j=0; j<ycol; j++){
      for(int n =0; n<xrow; n++){

        sumxy(i,j) += xminusmean(n,i) * yminusmean(n,j);
      }
    }}



  for(int i=0; i<xcol; i++){
    for(int j = 0; j<ycol; j++){
      covmat(i,j) = sumxy(i,j)/(xrow-1);
    }
  }

  return covmat;
}


// [[Rcpp::export]]
arma::vec nabla_LB_cpp(arma::vec lambda, arma::mat theta, int n_sources, int n_tracers, double beta_prior,
                           arma::mat concentrationmeans, arma::mat sourcemeans,
                           arma::mat correctionmeans,
                           arma::mat corrsds, arma::mat sourcesds, arma::mat y,
                           arma::vec c){

  int thetanrow = theta.n_rows;
  int lambdalength = lambda.n_elem;

  arma::mat big_c(thetanrow, c.n_elem);

  //working
  arma::mat big_delta_lqlt(thetanrow, lambdalength);
  arma::mat big_h_lambda_rep(lambdalength, thetanrow);
  arma::mat big_h_lambda_rep_transpose(thetanrow, lambdalength);

  arma::vec big_h_lambda(thetanrow);
  arma::vec big_h_lambda_transpose(thetanrow);

  for(int i = 0; i <thetanrow; i++){
  big_delta_lqlt.row(i) = delta_lqltcpp(lambda, theta.row(i).t(), 0.01, n_sources, n_tracers).t();
  }

  for(int i =0; i<thetanrow; i++){
    big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers, beta_prior,
                 concentrationmeans, sourcemeans,
                 correctionmeans,
                 corrsds,sourcesds, theta.row(i).t(), y,
                 lambda);
  }



  for(int i =0; i<lambdalength; i++){
    big_h_lambda_rep.row(i) = big_h_lambda.t();
  }

  for(int i=0; i<lambdalength; i++){
    for (int j=0; j < theta.n_rows; j++){
      big_h_lambda_rep_transpose(j,i) = big_h_lambda_rep(i,j);
    }}


  for(int i =0; i<thetanrow; i++){
    big_c.row(i) = c.t();
  }



  arma::mat big_h_minus_c(thetanrow, lambdalength);
  //arma::mat big_h_minus_c_t(lambda.length(), theta.nrow());

  for (int i = 0; i<thetanrow; i++){
    for(int j = 0; j<lambdalength; j++){
      big_h_minus_c(i,j) = big_h_lambda_rep_transpose(i,j) - big_c(i,j);
    }
  }

  //big_h_minus_c_t = transpose(big_h_minus_c);

  arma::mat ansmat(big_delta_lqlt.n_rows, big_h_minus_c.n_cols);

  for (int i = 0; i < big_delta_lqlt.n_rows; i++)
  {
    for (int j = 0; j < big_delta_lqlt.n_cols; j++) {


      ansmat(i,j) = big_delta_lqlt(i,j) * big_h_minus_c(i,j);


    }
  }

  arma::vec ans(ansmat.n_cols);
  for(int i = 0; i<ansmat.n_cols; i++){

    ans(i) = mean(ansmat.col(i));

  }

  return ans;
}




// [[Rcpp::export]]
arma::mat control_var_cpp(arma::vec lambda,
                          arma::mat theta,
                          int n_sources, int n_tracers,
                          double beta_prior,
                          arma::mat concentrationmeans,
                          arma::mat sourcemeans,
                          arma::mat correctionmeans,
                          arma::mat corrsds,
                          arma::mat sourcesds,
                          arma::mat y){

  int S = theta.n_rows;           // Number of rows in theta
  int lambdallength = lambda.n_elem;  // Length of lambda vector

  arma::mat big_delta_lqlt(S, lambdallength); // Matrix of size S x lambdallength
  arma::mat big_h_lambda_rep(S, lambdallength); // Matrix of size S x lambdallength for replication
  arma::vec big_h_lambda(S); // Vector of size S


  // Fill big_delta_lqlt matrix using delta_lqltcpp
  for(int i = 0; i < S; i++) {
    big_delta_lqlt.row(i) = delta_lqltcpp(lambda, theta.row(i).t(), 0.01, n_sources, n_tracers).t();
  }

 // Fill big_h_lambda vector using h_lambdacpp
 for(int i = 0; i < S; i++) {
   big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers, beta_prior,
                concentrationmeans, sourcemeans,
                correctionmeans, corrsds, sourcesds,
                theta.row(i).t(), y, lambda);
 }

 // Replicate big_h_lambda as columns in big_h_lambda_rep
 big_h_lambda_rep.each_col() = big_h_lambda;

 // Compute big_nabla matrix
 arma::mat big_nabla(S, lambdallength);
 for (int i = 0; i < S; i++) {
   for (int j = 0; j < lambdallength; j++) {
     big_nabla(i, j) = big_delta_lqlt(i, j) * big_h_lambda_rep(i, j);
   }
 }

 // Compute variance of big_delta_lqlt columns
 arma::vec var_big_delta_lqlt(lambdallength);
 for (int i = 0; i < lambdallength; i++) {
   var_big_delta_lqlt(i) = var(big_delta_lqlt.col(i));
 }

 // Compute covariance matrix
 arma::mat covmat = cov_mat_cpp(big_nabla, big_delta_lqlt);

 // Extract diagonal of covariance matrix
 arma::vec diag = covmat.diag();

 // Compute the final answer
 arma::vec ans(lambdallength);
 for (int i = 0; i < lambdallength; i++) {
   ans(i) = diag(i) / var_big_delta_lqlt(i);
 }

 return ans;

}



// [[Rcpp::export]]
double LB_lambda_cpp(arma::mat theta, arma::vec lambda, int n_sources, int n_isotopes,
                     double beta_prior,
                     arma::mat concentrationmeans, arma::mat sourcemeans,
                     arma::mat correctionmeans,
                     arma::mat corrsds, arma::mat sourcesds, arma::mat y){
  int S = theta.n_rows;

  arma::vec hlambdaapply(S);

  for(int i = 0; i <S; i++){
    hlambdaapply(i) = h_lambdacpp(n_sources, n_isotopes, beta_prior, concentrationmeans, sourcemeans,
                 correctionmeans, corrsds, sourcesds, theta.row(i).t(), y, lambda);
  }

  double ans = mean(hlambdaapply);

  return ans;


}



// [[Rcpp::export]]
arma::vec run_VB_cpp(arma::vec lambdastart,
                         int n_sources,
                         int n_tracers,
                         double beta_prior,
                         arma::mat concentrationmeans,
                         arma::mat sourcemeans,
                         arma::mat correctionmeans,
                         arma::mat corrsds,
                         arma::mat sourcesds,
                         arma::mat y,
                         int S,
                         int P,
                         double beta_1,
                         double beta_2,
                         int tau,
                         double eps_0,
                         int t_W,
                         bool solo
){



  int lsl = lambdastart.n_elem;

  arma::mat theta(S, (n_sources + n_tracers));

  theta = sim_thetacpp(S, lambdastart, n_sources, n_tracers, solo);

  arma::vec c(lsl);

  c = control_var_cpp(lambdastart, theta, n_sources, n_tracers, beta_prior,
                      concentrationmeans,
                      sourcemeans, correctionmeans,
                      corrsds, sourcesds, y);

  arma::vec g_0(lsl);

  arma::vec c_0(lsl);
  for(int i=0; i<lsl; i++){
    c_0(i) = 0;
  }

  g_0 = nabla_LB_cpp(lambdastart, theta,
                     n_sources, n_tracers,
                     beta_prior,
                     concentrationmeans, sourcemeans,
                     correctionmeans, corrsds,
                     sourcesds, y, c_0);

  arma::vec nu_0(lsl);

  for(int i = 0; i<lsl; i++){
    nu_0(i) = pow(g_0(i),2);
  }


  arma::vec g_bar(lsl);
  g_bar = g_0;

  arma::vec nu_bar(lsl);
  nu_bar = nu_0;

  arma::vec g_t(lsl);
  arma::vec nu_t(lsl);

  double patience = 0;
  bool stop = FALSE;
  double max_LB_bar = R_NegInf;
  double alpha_t = 0;
  double t = 0;
  arma::vec LB(t_W+1);

  for(int i=0; i<t_W; i++){
    LB(i) = NA_REAL;
  }

  arma::vec lambda(lsl);

  for(int i = 0; i<lsl; i++){
    lambda(i) = lambdastart(i);
  }

   while(stop == FALSE){

     theta = sim_thetacpp(S, lambda, n_sources, n_tracers, solo);

     g_t = nabla_LB_cpp(lambda, theta, n_sources, n_tracers, beta_prior, concentrationmeans,
                        sourcemeans, correctionmeans, corrsds, sourcesds,
                        y, c);

     c = control_var_cpp(lambda, theta,n_sources,n_tracers, beta_prior,
                         concentrationmeans, sourcemeans,
                         correctionmeans,
                         corrsds,sourcesds, y);

    for(int i=0; i<lsl; i++){
      nu_t(i) = pow(g_t(i),2);
    }

    for(int i=0; i<lsl; i++){
      g_bar(i) = (beta_1 * g_bar(i)) + ((1-beta_1) * g_t(i));
      nu_bar(i) = (beta_2 * nu_bar(i)) + ((1-beta_2) * nu_t(i));
    }

    arma::vec alpha_min(2);

    for(int i=0; i<2; i++){
      alpha_min(0) = eps_0;
      alpha_min(1) = eps_0 * (tau/(t+1));
    }

    alpha_t = arma::min(alpha_min);
  //
  //
  //
  //
    //# Update lambda
    for(int i = 0; i<lsl; i++){
      //lambda(i) = lambda(i) + alpha_t * 1/pow(nu_bar(i), 0.5);
      lambda(i) = lambda(i) + alpha_t * (g_bar(i)/(pow(nu_bar(i), 0.5)));
    }

    Rcout << "Iteration : " << t << "\n";


    //////////// This was written by Ahmed

    int r = t;
    int inn = 0;
    while(1){
      inn++;
      r = r/10;
      if(r == 0) break;
    }

    for(int j = 0 ; j < (13+inn) ;j++){
      Rcout<<"\b";
    }



     //# Compute the moving average LB if out of warm-up
     if(t<=t_W){



       LB(t) = LB_lambda_cpp(theta, lambda, n_sources, n_tracers,
          beta_prior,
          concentrationmeans, sourcemeans,
          correctionmeans,
          corrsds,sourcesds, y);
     }
     else{
       for (int i = 0; i<(t_W-1); i++){
         LB(i) = LB(i+1);
       }



      LB(t_W) = LB_lambda_cpp(theta, lambda, n_sources, n_tracers,
         beta_prior,
         concentrationmeans, sourcemeans,
         correctionmeans,
         corrsds,sourcesds, y);


      double LB_bar = mean(LB);

      arma::vec maxbar(2);

      for(int i=0; i<2; i++){
        maxbar(0) = max_LB_bar;
        maxbar(1) = LB_bar;
      }

      max_LB_bar = arma::max(maxbar);

      if(LB_bar>= max_LB_bar){
        patience = 0;
      } else{
        patience = patience +1;
      }

    }

    //////////////////////
    if(patience>P){
      stop = TRUE;
    }
    t = t + 1;
  }

  return lambda;

}

