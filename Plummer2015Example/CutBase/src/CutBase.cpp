#include <Rcpp.h>

using namespace Rcpp;

NumericVector Npart = {111, 71, 162, 188, 145, 215, 166, 37, 173,143, 229, 696, 93};
NumericVector Npop = {26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218};

inline double ljc(double x){  
  int i;
  double s=0;
  for(i=1;i<=x;i++) s+=log(i);
  return s;
}

double sinpx(double phi, double Z, double Npart){
  double result;
  double shapeOne = 1+Z;
  double shapeTwo = 1+Npart-Z;
  result = (shapeOne-1)*log(phi)+(shapeTwo-1)*log(1-phi)-log(R::beta(shapeOne,shapeTwo));
  
  return result;
}


double sinpy(double phi, double Y, double Npop, double thetaOne, double thetaTwo){
  double lambda;
  lambda = Npop*0.001*exp(thetaOne+thetaTwo*phi);
  
  return Y*log(lambda)-lambda-ljc(Y);
}



//[[Rcpp::export]]

double px(NumericVector phi, NumericVector Z){
     return sum(mapply(phi, Z, Npart, sinpx));
}


//[[Rcpp::export]]

double py(NumericVector Y, NumericVector theta, NumericVector phi){
  double result=log(1)-log(200);
  for(int i = 0; i < 13; i++)
  {
    result += sinpy(phi[i],Y[i],Npop[i],theta[0],theta[1]);
  }
  
  return(result);
}
