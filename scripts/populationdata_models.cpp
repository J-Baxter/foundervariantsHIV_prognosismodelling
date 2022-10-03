#include <Rcpp.h>

using namespace Rcpp;

// C++ version of populationdata_models.R
// for fixed spVL, per virion prob of infection and proportion of exposure with potential for infection this
// function calculates i) the probability of infection, ii) the probability of multiple founder variants

// [[Rcpp::export]]

List populationmodel_fixedVL_cpp(double sp_ViralLoad = 10**5, 
                                 double PerVirionProbability = 4.715*1e-8, 
                                 double PropExposuresInfective = 0.029){
  
  //this is to speed up the code, but can be updated based on phylogenetic results
  
  int nSims = 10000;
  
  //rename for reducing code length
  f = PropExposuresInfective; //fraction of exposure that can lead to infection
  p = PerVirionProbability;
    

  //parameters from the Thomson paper for ....
  double alpha = -3.55;
  double sigma = 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))));
  double mu = 4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)));
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
  double tauc = Dmax*(Dfifty**Dk)/(nc**Dk + (Dfifty**Dk));
  double maximumTime = taup + taua + tauc;
  
  double totalOfTimeVals = taup + tauc + taua;
  
  //NOW CALCULATE FULL DISTRIBUTION OF VIRIONS NUMBER
  double m = np*p;
  double s = np*p*(1-p);
  double iter = 1;

  
  double probNoTransmissionPerSexAct = 0;
  double probTransmissionPerSexAct = 0;
  double logBitToAddOn = 0;

  probNoTransmissionPerSexAct = probNoTransmissionPerSexAct + ((1 - f) + f*((taup/(taup + tauc + taua))*(1 - p)**np + (tauc/(taup + tauc + taua))*(1 - p)**nc + (taua/(taup + tauc + taua))*(1 - p)**na));
  probTransmissionPerSexAct = 1 - probNoTransmissionPerSexAct;
  
  double probTransmitnparticles[10**7] = { 0 };
  double threshold = 10**(-6);
  double n = 1;

  
  while (sum((probTransmitnparticles)/(probTransmissionPerSexAct)) < (1 - threshold)){
    
    double p0 = 0;
    
    if (n <= np){
      double logBitToAddOn = log(taup/(taup + tauc + taua)) + n*log(p) + (np - n)*log(1-p);
      
      if (n > np){
  
        for (z = (np - n + 1):np){
          logBitToAddOn = logBitToAddOn + log(z);}
        
        for (z = 1:n){
          logBitToAddOn = logBitToAddOn - log(z);}
      }
      
      p0 = p0 + exp(logBitToAddOn)*f;
    }
    
    if (n <= nc){
      logBitToAddOn = log(tauc/(taup + tauc + taua)) + n*log(p) + (nc - n)*log(1-p);
      
      if (nc > n){
        
        for (z = (nc-n+1):nc){
          logBitToAddOn = logBitToAddOn + log(z);}
        
        for (z = 1:n){
          logBitToAddOn = logBitToAddOn - log(z);}
      }
      
      p0 = p0 + exp(logBitToAddOn)*f
    }
    
    if (n <= na){
      logBitToAddOn = log(tauc/(taup + tauc + taua)) + n*log(p) + (na - n)*log(1-p);
      
      if (na > n){
        
        for (z = (na-n+1):na){
          logBitToAddOn = logBitToAddOn + log(z);}//
        
        for (z = 1:n){
          logBitToAddOn = logBitToAddOn - log(z);}//
        
      }
      
      p0 = p0 + exp(logBitToAddOn)*f
    }
    
    probTransmitnparticles(n) = p0; //should be parenthesis for vector?
    n = n + 1;
  }
  
  probTransmitnparticles = probTransmitnparticles(1:(n-1));
  probTransmitnparticles = probTransmitnparticles./sum(probTransmitnparticles);
  
  
  //Now, the probability of transmitting N variants vs N
  
  nparticlesConsidered = length(probTransmitnparticles);
  maximumTime = max(taup + tauc + taua);
  nTimeSteps = 1000;
  timeWindowEdges = [0:maximumTime/nTimeSteps:maximumTime];
  timeVals = zeros(nTimeSteps,1);
  probTransmitNvariantsGivenTimeTAndTransmitNparticles = zeros(nparticlesConsidered, nTimeSteps, nparticlesConsidered);
  timeSinceInfectionBeingCalculated = 0;
  
  for (timeV = 1:nTimeSteps){
    timeVals(timeV) = (timeWindowEdges(timeV) + timeWindowEdges(timeV + 1))/2;
    
    if (timeVals(timeV) > timeSinceInfectionBeingCalculated){
      timeSinceInfectionBeingCalculated = timeSinceInfectionBeingCalculated + 1;
      
    }
    
    probDist = zeros(maxnvariants,1); //
    
    for (i = 1:maxnvariants){
      probDist(i) = gampdf(i, 0.417, timeVals(timeV)/0.563);//vector []
      
    }
    
    probDist = probDist/sum(probDist);
    
    nSims = 100000;
    
    for (simNo = 1:nSims){
      variantIIndicator = zeros(maxnvariants,1);
      
      for(nparticles= 1:nparticlesConsidered){
        variantPicked = 1;
        variantPickedCount = 0;
        r = rand();
        
        while (r>variantPickedCount){
          variantPickedCount = variantPickedCount + probDist(variantPicked);//
          variantPicked = variantPicked + 1;
        }
        
        variantIIndicator(variantPicked -1) = 1;
        
        nvariantsTransferred = sum(variantIIndicator);
        probTransmitNvariantsGivenTimeTAndTransmitNparticles(nvariantsTransferred, timeV, nparticles) = probTransmitNvariantsGivenTimeTAndTransmitNparticles(nvariantsTransferred, timeV, nparticles) + 1;
      }
    }
  }
  
  probTransmitNvariantsGivenTimeTAndTransmitNparticles = probTransmitNvariantsGivenTimeTAndTransmitNparticles./nSims;
  probTransmitMvariants = zeros(min(maxnvariants, nparticlesConsidered),1); // double probTransmitnparticles[10**7] = { 0 };
  
  double integralPrimary = 0;
  double integralChronic = 0;
  double integralPreAids = 0;
  
  for (M = 1:min(maxnvariants, nparticlesConsidered)){
    
    for (j = 1:nTimesteps){
      
      if(timeVals(j) <= taup){
        
        for (nVal = M:min(np, nparticlesConsidered)){
          binomialBit = 1;
          
          if(np > nVal){
            logBinBit = nVal*log(p) + (np - nVal)*log(1-p);
            
            for (z = (np - nVal + 1):np){
              logBinBit = logBinBit + log(z);
            }
            
            for (z = 1:nVal){
              logBinBit = logBinBit - log(z);
            }
            
            binomialBit = exp(logBinBit);
          }
          
          integralPrimary = integralPrimary + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc + taua))*binomialBit;
        }
      }
    }
    
   
    for (j = 1:nTimesteps){
      
      if((timeVals(j) > taup) && (timeVals(j) <= (taup + tauc))){
        
        for (nVal = M:min(nc, ncarticlesConsidered)){
          binomialBit = 1;
          
          if(nc > nVal){
            logBinBit = nVal*log(p) + (nc - nVal)*log(1-p);
            
            for (z = (nc - nVal + 1):nc){
              logBinBit = logBinBit + log(z);
            }
            
            for (z = 1:nVal){
              logBinBit = logBinBit - log(z);
            }
            
            binomialBit = exp(logBinBit);
          }
          
          integralChronic = integralChronic + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitncarticles(M,j,nVal)*(1/(taup + tauc + taua))*binomialBit;
        }
      }
    }
    
    
    for (j = 1:nTimesteps){
      
      if((timeVals(j) > taup) && (timeVals(j) <= (taup + tauc + taua))){
        
        for (nVal = M:min(na, naarticlesConsidered)){
          binomialBit = 1;
          
          if(na > nVal){
            logBinBit = nVal*log(p) + (na - nVal)*log(1-p);
            
            for (z = (na - nVal + 1):na){
              logBinBit = logBinBit + log(z);
            }
            
            for (z = 1:nVal){
              logBinBit = logBinBit - log(z);
            }
            
            binomialBit = exp(logBinBit);
          }
          
          integralPreAids = integralPreAids + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitnaarticles(M,j,nVal)*(1/(taup +  taua))*binomialBit;
        }
      }
    }
    
    probTransmitMvariants(M) = probTransmitMvariants(M) + (integralPrimary + integralChronic + integralPreAids);
  }
  
  probTransmitMvariants = probTransmitMvariants/sum(probTransmitMvariants);
  //probMultivariantTransmission = 1 - probTransmitMvariants(1);
  
  
  output = list(probTransmissionPerSexAct, probTransmitMvariants);
  
  return output;
}
