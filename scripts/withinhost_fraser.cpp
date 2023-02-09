#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// C++ version of final model
// Erlang shape fixed at 4 for infectious categories and 2 for recovery categories

// [[Rcpp::export]]
List Fraser_cpp_model(NumericVector times,NumericVector state,
                   List parms){
  
  // Declare State Variables
  double cd4 = state["cd4"]; // Uninfected quiescent CD4 cells
  double cd8 = state["cd8"]; // Uninfected quiescent CD8 cells
  double cd4_activated = state["cd4_activated"]; // Uninfected activated CD4 cells
  double cd8_activated = state["cd8_activated"]; // Uninfected activated CD8 cells
  double infectious_active = state["infectious_active"]; //Actively infectious cells
  double infectious_latent = state["infectious_latent"]; //Latent infectious cells
  double v = state["v"]; //free virus
  double z = state["z"]; //anti-HIV CTLs
  double k_cd4 = state["k_cd4"];
  double k_cd8 = state["k_cd8"];
  double pk_cd4 = state["k_cd4"];
  double pk_cd8 = state["k_cd8"];
 
  

  // Declare Parms
  double lambda_cd4 = parms["lambda_cd4"]; //daily thymic production of new cd4
  double lambda_cd8 = parms["lambda_cd8"]; //daily thymic production of new cd8
  double a_0 = parms["a_0"]; //average rate of T cell activation per antigenic exposure
  double mu = parms["mu"]; //daily rate of non-antigen-driven homeostaatic T cell division
  double x_s = parms["x_s"]; //relative t cell pool size, below which t cell activation fails due to exhaustion
  double mu_a = parms["mu_a"]; // activated t cell division rate
  double p_a = parms["p_a"]; // average probability of an activated t cell successfully dividing in an HIV negative control
  double beta = parms["beta"]; // average per virion infection rate of an activated CD4+ T cell
  double alpha = parms["alpha"]; //death rate of productively infected cell (in the absence of CTL)
  double alpha_latent = parms["alpha_latent"]; //rate of reactivation of latent infected cells
  double frac_latent = parms["frac_latent"]; //proportion of successful infections that result in latency
  double p = parms["p"]; //maximum proliferation rate of anti-HIV CTLs
  double d = parms["d"]; //death rate of anti-HIV CTLs
  double z_0 = parms["z_0"]; //pre-infection frequency of anti-HIV CTLs
  double infectious_0 = parms["infectious_0"];
  double bc_ratio = parms["b"]; // ratio of viral production rate in productively infected cells and viral lifetime
  double sigma = parms["sigma"]; // max rat of CTL killing of HIV infected cells
  //double N_PB = parms["N_PB"];
  double theta_cd4 = parms["theta_cd4"]; // average clearance rate in antigenic exposure model
  double theta_cd8 = parms["theta_cd8"]; // average exposure rate in antigenic exposure model
  double e_cd4 = parms["e_cd4"];
  double e_cd8 = parms["e_cd8"];

  
  //Antigenic Exposure and T Cell Activation
  double a_cd4 = a_0*k_cd4*(cd4/(cd4+x_s));
  double a_cd8 = a_0*k_cd8*(cd8/(cd8+x_s));
  
  double dpk_cd4 = theta_cd4*(k_cd4+1)*(pk_cd4+1) - (e_cd4 + theta_cd4 * k_cd4)*pk_cd4 + e_cd4*(pk_cd4-1); //Ask Katie
  double dpk_cd8 = theta_cd8*(k_cd8+1)*(pk_cd8+1) - (e_cd8 + theta_cd8 * k_cd8)*pk_cd8 + e_cd8*(pk_cd8-1);
  double epsilon = infectious_active/(infectious_active+ infectious_0);
  double homeostasis = lambda_cd4 + lambda_cd8 + mu + 2*a_0*(2*p_a + 1);

  // T Cell Dynamics
  double dcd4 = lambda_cd4-a_cd4 + 2*p_a*mu_a*cd4_activated + mu*cd4 - homeostasis*cd4*(cd4 + cd8);
  double dcd8 = lambda_cd8-a_cd8 + 2*p_a*mu_a*cd8_activated + mu*cd8 - homeostasis*cd8*(cd4 + cd8);
  double dcd4_activated = a_cd4 - mu_a*cd4_activated - beta*v*cd4_activated;
  double dcd8_activated = a_cd8 - mu_a*cd8_activated;
  
  //HIV Viral Dynamics and CTL Control of HIV
  double dinfectious_active = (1-frac_latent)*beta*v*cd4_activated + alpha_latent*infectious_latent - (alpha + sigma*epsilon*z)*infectious_active;
  double dinfectious_latent = frac_latent*beta*v*cd4_activated - alpha_latent*infectious_latent;
  double dz = d*(z_0 - z) + p*epsilon*z;
  double dv = infectious_active*bc_ratio;
  
  NumericVector res_vec = NumericVector::create(dcd4, dcd8, dcd4_activated, dcd8_activated, dinfectious_active, dinfectious_latent, dz, dv, k_cd4, k_cd8, dpk_cd4, dpk_cd8);
  
  List res = List::create(res_vec);
  
  
  return res;
}
