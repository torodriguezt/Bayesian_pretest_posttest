functions{
  real bivbin_lpdf(vector theta, real a1,real a2,real a3){
    real res;
    res= lgamma((a1+a2+a3))-
    lgamma(a2)-lgamma(a3)-lgamma(a1)+
    (a2-1)*log(theta[1])+(a3+a1-1)*log(1-theta[1])+
    (a3-1)*log(theta[2])+(a2+a1-1)*log(1-theta[2])-
    (a1+a2+a3)*log(1-theta[1]*theta[2]);
    return res; 
}
}

data {
  // dimension parametros y número de muestras
  int<lower=0> P;
  // variable 
  int X[P];
  real alpha1;
  real alpha2;
  real alpha3;
 
  int n[P];
  }
  

// parametros del modelo
parameters {
  
  // vector de parametros

  vector<lower=0,upper=1> [P] Theta;
}

model {

 // target += bivbin_lpdf(Theta|alpha1,alpha2,alpha3);
Theta ~ bivbin(alpha1,alpha2,alpha3);
//  for (i in 1:P){
//  X[i] ~ binomial(n[i],Theta[i]);
//  }
X ~ binomial(n,Theta);
}
