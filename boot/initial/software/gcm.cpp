#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Lrel);        // length at release (tags)
  DATA_VECTOR(Lrec);        // length at recapture (tags)
  DATA_VECTOR(liberty);     // time at liberty (tags)
  DATA_VECTOR(Aoto);        // age (otoliths)
  DATA_VECTOR(Loto);        // length (otoliths)
  DATA_SCALAR(Lshort);      // length where sd(length) is sigma_1
  DATA_SCALAR(Llong);       // length where sd(length) is sigma_2
  DATA_SCALAR(Afix);        // age where Lfix starting point can be provided
  int Ntags = Lrel.size();  // number of tags
  int Noto = Loto.size();   // number of otoliths

  PARAMETER(log_Linf);        // asymptotic length
  PARAMETER(log_rmax);        // initial slope
  PARAMETER(log_k);           // growth coefficient
  PARAMETER(log_Lfix);        // length at Afix
  PARAMETER(log_sigma_1);     // sd(length) at Lshort
  PARAMETER(log_sigma_2);     // sd(length) at Llong
  PARAMETER_VECTOR(log_age);  // estimated age at release

  Type Linf = exp(log_Linf);
  Type rmax = exp(log_rmax);
  Type k = exp(log_k);
  Type Lfix = exp(log_Lfix);
  Type L0 = Lfix - rmax * Afix;
  Type t50 = log(exp(k * (Linf-L0) / rmax) - Type(1.0)) / k;
  Type sigma_1 = exp(log_sigma_1);
  Type sigma_2 = exp(log_sigma_2);
  Type sigma_slope = (sigma_2 - sigma_1) / (Llong - Lshort);  // s = a + b*age
  Type sigma_intercept = sigma_1 - Lshort * sigma_slope;

  vector<Type> age(Ntags);
  vector<Type> Lrel_hat(Ntags);
  vector<Type> Lrec_hat(Ntags);
  vector<Type> Loto_hat(Noto);
  vector<Type> sigma_Lrel(Ntags);  // sigma for each Lrel datapoint
  vector<Type> sigma_Lrec(Ntags);  // sigma for each Lrec datapoint
  vector<Type> sigma_Loto(Noto);   // sigma for each Loto datapoint
  vector<Type> nll_Lrel(Ntags);    // nll for each Lrel datapoint
  vector<Type> nll_Lrec(Ntags);    // nll for each Lrec datapoint
  vector<Type> nll_Loto(Noto);     // nll for each Loto datapoint
  vector<Type> age_seq(10*365+1);  // ages to calculate growth curve
  vector<Type> curve(10*365+1);    // growth curve for ages 0->10 (confint)

  age=exp(log_age);
  Lrel_hat = L0 + rmax *
    ((log(Type(1.0) + exp(-k*t50)) -
      log(Type(1.0) + exp(k*(age-t50)))) / k + age);
  Lrec_hat = L0 + rmax *
    ((log(Type(1.0) + exp(-k*t50)) -
      log(Type(1.0) + exp(k*((age+liberty)-t50)))) / k + age+liberty);
  Loto_hat = L0 + rmax *
    ((log(Type(1.0) + exp(-k*t50)) -
      log(Type(1.0) + exp(k*(Aoto-t50)))) / k + Aoto);
  sigma_Lrel = sigma_intercept + sigma_slope * Lrel_hat;
  sigma_Lrec = sigma_intercept + sigma_slope * Lrec_hat;
  sigma_Loto = sigma_intercept + sigma_slope * Loto_hat;

  Type nll = 0.0;
  nll_Lrel = -dnorm(Lrel, Lrel_hat, sigma_Lrel, true);
  nll_Lrec = -dnorm(Lrec, Lrec_hat, sigma_Lrec, true);
  nll_Loto = -dnorm(Loto, Loto_hat, sigma_Loto, true);
  nll = sum(nll_Lrel) + sum(nll_Lrec) + sum(nll_Loto);

  for(int i=0; i<=10*365; i++)
    age_seq(i) = Type(i) / Type(365.0);
  curve = L0 + rmax *
    ((log(exp(-k*t50) + Type(1.0)) - log(exp(k*(age_seq-t50)) +
                                         Type(1.0))) / k + age_seq);

  ADREPORT(curve);
  REPORT(L0);
  REPORT(rmax);
  REPORT(k);
  REPORT(t50);
  REPORT(Linf);
  REPORT(Afix);
  REPORT(Lfix);
  REPORT(age);
  REPORT(liberty);
  REPORT(Lrel);
  REPORT(Lrec);
  REPORT(Aoto);
  REPORT(Loto);
  REPORT(Lrel_hat);
  REPORT(Lrec_hat);
  REPORT(Loto_hat);
  REPORT(Lshort);
  REPORT(Llong);
  REPORT(sigma_1);
  REPORT(sigma_2);
  REPORT(sigma_Lrel);
  REPORT(sigma_Lrec);
  REPORT(sigma_Loto);
  REPORT(nll_Lrel);
  REPORT(nll_Lrec);
  REPORT(nll_Loto);
  REPORT(curve);

  return nll;
}
