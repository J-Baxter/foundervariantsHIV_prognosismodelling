#include <Rcpp.h>

using namespace Rcpp;

// C++ version of populationdata_models.R
// for fixed spVL, per virion prob of infection and proportion of exposure with potential for infection this
// function calculates i) the probability of infection, ii) the probability of multiple founder variants

// [[Rcpp::export]]

List populationmodel_fixedVL_cpp(double sp_ViralLoad = 10^5, 
                                 double PerVirionProbability = 4.715*1e-8, 
                                 double PropExposuresInfective = 0.029){
  
  //this is to speed up the code, but can be updated based on phylogenetic results
  
  int nSims = 10000;
  
  //rename for reducing code length
  f <- PropExposuresInfective; //fraction of exposure that can lead to infection
    

  //parameters from the Thomson paper for ....
  double alpha = -3.55;
  double sigma = 0.78/(sqrt(1 - ((2*alpha^2)/(pi*(1 + alpha^2)))));
  double mu = 4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha^2)));
  double Vp = 87173000;
  double Va = 24004000;
  double kappa = 1;
  double np = round(kappa*Vp);
  double na = round(kappa*Va);
  double Dmax = 25.4;
  double Dfifty = 3058;
  double Dk = 0.41;
  double nc = round(sp_ViralLoad);
  double taup = 0.24;
  double taua = 0.75;
  double tauc = Dmax*(Dfifty^Dk)/(nc^Dk + (Dfifty^Dk));
  double maximumTime = taup + taua + tauc;
  
  double totalOfTimeVals = taup + tauc + taua;
  
  //NOW CALCULATE FULL DISTRIBUTION OF VIRIONS NUMBER
  double m = np*PerVirionProbability;
  double s = np*PerVirionProbability*(1-PerVirionProbability);
  double iter = 1;
  

  primary_prob_fulldist = R::dbinom(0:ceiling(m+iter*s), np, PerVirionProbability);

    
    
  List output =list(probTransmissionPerSexAct = probTransmissionPerSexAct, variant_distribution = variant_distribution);
  
  return output;
}
