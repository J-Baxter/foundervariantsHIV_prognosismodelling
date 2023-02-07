#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// C++ version of final model
// Erlang shape fixed at 4 for infectious categories and 2 for recovery categories

// [[Rcpp::export]]
List SIR_cpp_model(NumericVector times,NumericVector Y,
                   List parms){
  

  // Declare Parms
  double gamma_cd4 = parms["gamma_cd4"];
  double gamma_cd8 = parms["gamma_cd8"];
  double a_0 = parms['a_0'];
  double mu = parms['mu'];
  double x_s = parms['x_s'];
  double mu_a = parms['mu_a'];
  double p_a = parms['p_a'];
  double beta = parms['beta'];
  double alpha = parms['alpha'];
  double alpha_latent = parms['alpha_latent'];
  double f_latent = parms['f_latent'];
  double p = parms['p'];
  double d = parms['d'];
  double z_0 = parms['z_0'];
  double bc = parms['bc'];
  double sigma = parms['sigma'];
  double N_PB = parms['N_PB'];
  double theta_cd4 = parms['theta_cd4'];
  double theta_cd8 = parms['theta_cd8'];
  double e_cd4 = parms['e_cd4'];
  double e_cd8 = parms['e_cd8'];
  
  double v = beta*y/c;
  double a_4 = a_0*k_4*(x4/(x4+x_s));
  double a_8 = a_0*k_8*(x8/(x8+x_s));
  double epsilon = y/(y+y_t);

  double dx4 = gamma_cd4 - a_4 + 2*p_a*mu_a*x_4a + mu*x4 - ;
  double dx8 = gamma_cd8 - a_8 + 2*p_a*mu_a*x_8a + mu*x8;
  double dx_4a = a_4 - mu_a*x_4a - beta*v*x_4a;
  double dx_8a = a_8 - mu_a*x_8a;
  double dy = (1-f_latent)*beta*v*x_4a + alpha_latent*y_l - (alpha + sigma*epsilon*z)*y;
  double dy_l = f_l*beta*v*x_4a - alpha_latent*y_l;
  double dz = d(z_0 - z) + p*epsilon*z;
  
  NumericVector res_vec = NumericVector::create(dx4, dx8, dx_4a, dx_8a, dy, dy_l, dz);
  
  List res = List::create(res_vec);
  
  
  return res;
}
