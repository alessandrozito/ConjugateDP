#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// In this file, we run the simulations for our ConjugateDP prior.

// Multivariate T distribution - Very useful
// [[Rcpp::export]]
double logdT(arma::vec x, double nu, arma::vec mu, arma::mat Sigma) {
  // Log density for the multivatiate t distribution
  double p = Sigma.n_cols;  // Dimension of the vectors
  double pi = 3.141592653589793238463 ;
  // Log determinant of sigma
  double logdetSigma = arma::log_det_sympd(Sigma);
  // Mahalanobis distance
  arma::mat Sigma_inv = arma::inv_sympd(Sigma);
  arma::vec x_center= x - mu;
  arma::mat dist(1,1);
  dist = x_center.t() * Sigma_inv * x_center;
  // Value of the log density
  double log_dt = lgamma((nu + p)/2) - lgamma(nu/2) -.5*logdetSigma -.5*p*log(nu*pi) -.5*(nu + p)*log(1 + dist(0,0)/nu);
  return log_dt;
}


//[[Rcpp::export]]
arma::uword get_id_cl(arma::vec cluster_names, arma::uword cl){
  // Get the cluster names
  arma::uvec id_cl = find(cluster_names == cl);
  return id_cl(0);
}

// Functions to pack the matrix into a vector for better storage
// [[Rcpp::export]]
arma::mat vec2mat(arma::vec x, int p) {
  arma::mat y(x);
  y.reshape(p, p);
  return y;
}

// Functions to unpack the vector into a matrix for posterior computations
// [[Rcpp::export]]
arma::vec mat2vec(arma::mat M) {
  arma::vec y = M.as_col();
  //y.reshape(nrow, ncol);
  return y;
}


// Sample the Stirling-gamma
double rStirgamma(double a, double b, double m){
  // Obtaining namespace of ConjugateDP package
  Environment pkg = Environment::namespace_env("ConjugateDP");
  // Picking up rStirgamma() function
  Function f = pkg["rStirgamma"];
  SEXP sample = f(Named("a")= a, Named("b")=b, Named("m")=m, Named("nsamples")= 1);
  return *REAL(sample);
}


// Remove observation from the posterior parameters
// [[Rcpp::export]]
void RemoveObservation(arma::uword id_cl,
                       arma::vec y_obs,
                       arma::vec& cl_counts,
                       arma::vec& cluster_names,
                       arma::vec& nu,
                       arma::vec& kappa,
                       arma::mat& y_bar,
                       arma::mat& mu,
                       arma::mat& S,
                       arma::vec& mu0,
                       double& kappa0){

  // Substract 1 from the counts
  cl_counts(id_cl) -= 1;
  // Check if the cluster needs to be removed
  if (cl_counts(id_cl) == 0) {
    // Remove the entire cluster from the current state
    // Shed matrices
    y_bar.shed_col(id_cl);
    mu.shed_col(id_cl);
    S.shed_col(id_cl);
    // Shed vectors
    arma::uvec sel = arma::find(cl_counts > 0);
    cl_counts = cl_counts(sel);
    cluster_names = cluster_names(sel);
    nu = nu(sel);
    kappa = kappa(sel);
  } else {
    // Calculate useful quantities
    double n = cl_counts(id_cl) + 1;
    double n_old =  cl_counts(id_cl);
    double kappa_old = kappa(id_cl) - 1;
    arma::mat S_temp = vec2mat(S.col(id_cl), y_obs.n_elem);
    arma::vec y_bar_old = (n * y_bar.col(id_cl) - y_obs)/n_old;
    arma::vec y_minus_mu0 = y_bar.col(id_cl) - mu0;
    arma::vec yold_minus_mu0 = y_bar_old - mu0;
    // Decrease the hyperparameters
    nu(id_cl) -= 1;
    mu.col(id_cl) = (kappa(id_cl) / kappa_old) * mu.col(id_cl) - y_obs / kappa_old;
    S_temp = S_temp +
      n * y_bar.col(id_cl) * y_bar.col(id_cl).t() -
      n_old * y_bar_old * y_bar_old.t() -
      (kappa0 * n / (kappa0 + n)) * y_minus_mu0 * y_minus_mu0.t() +
      (kappa0 * n_old / (kappa0 + n_old)) * yold_minus_mu0 * yold_minus_mu0.t() -
      y_obs * y_obs.t();
    S.col(id_cl) = mat2vec(S_temp);
    // Decrease the hyperparameters
    kappa(id_cl) -= 1;
    // Remove from the mean
    y_bar.col(id_cl) = y_bar_old;
  }
}

// Add observation to the posterior
// [[Rcpp::export]]
void AddObservation(double& sampled_cl,
                    arma::vec& y_obs,
                    arma::vec& cl_counts,
                    arma::vec& cluster_names,
                    arma::vec& nu,
                    arma::vec& kappa,
                    arma::mat& y_bar,
                    arma::mat& mu,
                    arma::mat& S,
                    arma::vec& mu0,
                    arma::mat& S0,
                    double& kappa0,
                    double& nu0){

  // The sampled cluster is new. We need to include it everywhere
  double m = max(cluster_names);
  int K = cluster_names.n_elem;

  // New cluster
  if (sampled_cl > m) {
    // Cluster names and counts
    cluster_names.resize(K + 1);
    cluster_names(K) = sampled_cl;
    cl_counts.resize(K + 1);
    cl_counts(K) = 1;
    // nu
    nu.resize(K + 1);
    nu(K) = nu0 + 1;
    // kappa0
    kappa.resize(K + 1);
    kappa(K) = kappa0 + 1;
    // y_bar
    y_bar = arma::join_horiz(y_bar, y_obs);
    // mu
    mu = arma::join_horiz(mu, (mu0 * kappa0 + y_obs)/(kappa0 + 1));
    // S
    arma::mat S_temp;
    arma::mat sigma2 = y_obs.t() * y_obs;
    S_temp = S0 + sigma2(0,0) - y_obs * y_obs.t() +
      kappa0/(kappa0 + 1) * (y_obs - mu0) * (y_obs - mu0).t();
    S = arma::join_horiz(S, mat2vec(S_temp));
  } else {
    // New cluster
    arma::uword id_cl = get_id_cl(cluster_names, sampled_cl);
    double n = cl_counts(id_cl);
    double n_new = n + 1;
    double kappa_new = kappa(id_cl) + 1;
    arma::vec y_bar_new = (n * y_bar.col(id_cl) + y_obs)/n_new;
    nu(id_cl) += 1;
    // Posterior for mu
    mu.col(id_cl) =  kappa(id_cl)/kappa_new * (mu.col(id_cl) + y_obs/kappa(id_cl));
    // Posterior for S
    arma::vec y_minus_mu0 = y_bar.col(id_cl) - mu0;
    arma::vec ynew_minus_mu0 = y_bar_new - mu0;
    arma::mat S_temp = vec2mat(S.col(id_cl), y_obs.n_elem);
    S_temp = S_temp +
     n * y_bar.col(id_cl) * y_bar.col(id_cl).t() -
     n_new * y_bar_new * y_bar_new.t() -
     (kappa0 * n / (kappa0 + n)) * y_minus_mu0 * y_minus_mu0.t() +
     (kappa0 * n_new / (kappa0 + n_new)) * ynew_minus_mu0 * ynew_minus_mu0.t() +
     y_obs * y_obs.t();
    S.col(id_cl) = mat2vec(S_temp);
    // Update the rest
    y_bar.col(id_cl) = y_bar_new;
    cl_counts(id_cl) += 1;
    kappa(id_cl) += 1;
  }

}

// Sample the cluster
// [[Rcpp::export]]
double sampleCluster(arma::vec& y_obs,
                   int& p,
                   arma::vec cluster_names,
                   double& alpha,
                   arma::vec& cl_counts,
                   arma::vec& nu, arma::vec& kappa,
                   arma::mat& y_bar,
                   arma::mat& mu, arma::mat& S,
                   arma::mat& S0, arma::vec& mu0, double& kappa0, double& nu0) {
  int K = cluster_names.n_elem;
  arma::vec probs(K + 1);
  cluster_names.resize(K + 1);
  cluster_names(K) = max(cluster_names) + 1;

  // Calculate the probability for the existing clusters
  for (int j = 0; j < K; j++) {
    probs(j) = log(cl_counts(j)) +
      logdT(y_obs, nu(j) - p + 1, mu.col(j), vec2mat(S.col(j), 2)/(nu(j) - p + 1) * (kappa(j) + 1)/kappa(j));
  }
  // Calculate the probabilities for the new cluster
  probs(K) = log(alpha) +
    logdT(y_obs, nu0 - p + 1, mu0, S0/(nu0 - p + 1) * (kappa0 + 1)/kappa0);
  // Renormalize probabilities
  probs = exp(probs - max(probs));
  probs = probs/sum(probs);
  // Sample
  double sampled_cl = RcppArmadillo::sample(cluster_names, 1L, false, probs)[0];
  return sampled_cl;
}


// Gibbs sampler for Conjugate multivariate normal
// [[Rcpp::export]]
List ConjugateNormalMixture_cpp(arma::mat y, arma::vec clusters,
                                List pars, List prior, int iter,
                                double alpha = 1, bool random_alpha = true,
                                double a = 5, double b = 1, bool verbose = true) {
  int n = y.n_rows;
  int p = y.n_cols;

  y = y.t(); // Transpose the y for now
  // Unpack the prior parameters
  double nu0 = prior["nu0"];
  double kappa0 = prior["kappa0"];
  arma::vec mu0 = prior["mu0"];
  arma::mat S0 = prior["S0"];

  // Unpack the posterior parameters for each cluster, as initialized in R
  arma::vec cl_counts = pars["cl_counts"];
  arma::vec nu = pars["nu"];
  arma::vec kappa = pars["kappa"];
  arma::mat y_bar = pars["y_bar"];
  arma::mat mu = pars["mu"];
  arma::mat S = pars["S"];
  arma::vec cluster_names = pars["cluster_names"];

  // Storage quantities
  arma::mat CLUSTERS(iter, n);
  arma::vec ALPHA(iter);

  // Define useful quantities now
  arma::uword id_cl;
  arma::uword cl;
  arma::vec y_obs(p);
  double sampled_cl;
  int R_show = iter/10;
  if (R_show == 0) {
    R_show = 1;
  }
  // Begin sampling now
  for (int r = 0; r < iter; r++) {
   if(((r+1)%R_show==0) and (verbose == true)) {
     int state = round(100 * r/iter);
     Rprintf("Iteration: %i [%i%%] \n", r+1,  state);
   }
    // Gibbs sampling - Algorithm 3 of Neal (2001)
    for (int i = 0; i < n; i++) {
      // Extract cluster info
      cl = clusters(i);
      id_cl = get_id_cl(cluster_names, cl);
      y_obs = y.col(i);
      // Remove observation from the posterior
      RemoveObservation(id_cl = id_cl, y_obs = y_obs, cl_counts = cl_counts,
                        cluster_names = cluster_names, nu = nu, kappa = kappa,
                        y_bar = y_bar, mu = mu, S = S, mu0 = mu0, kappa0 = kappa0);

      // Sample the cluster
      sampled_cl = sampleCluster(y_obs = y_obs, p = p,
                                cluster_names = cluster_names,
                                alpha = alpha,
                                cl_counts = cl_counts,
                                nu = nu, kappa = kappa,
                                y_bar = y_bar,
                                mu = mu, S = S,
                                S0 = S0, mu0 = mu0, kappa0 = kappa0,  nu0 = nu0);
      CLUSTERS(r, i) = sampled_cl;
      //Rcout << sampled_cl << "\n";
      //sampled_cl = 11;
      // Update the posterior
      AddObservation(sampled_cl =  sampled_cl,
                    y_obs =  y_obs, cl_counts = cl_counts,
                    cluster_names = cluster_names, nu =  nu,
                    kappa = kappa, y_bar= y_bar, mu = mu,
                    S = S, mu0 = mu0,
                    S0 = S0, kappa0 = kappa0, nu0 = nu0);
    }

    // Sample alpha through the Stirling-gamma
    if(random_alpha == true) {
      alpha = rStirgamma(a + cluster_names.n_elem, b + 1, n);
    }
    ALPHA(r) = alpha;

  }
  return(List::create(_["Cluster_chain"] = CLUSTERS,
                      _["alpha"] = ALPHA));

}


/*** R
set.seed(15)
out <- ConjugateNormalMixture_cpp(y = y, clusters = clusters, pars = pars,
                                  prior = prior, iter = 1, verbose = FALSE)
out$Cluster_chain
out$alpha
*/


