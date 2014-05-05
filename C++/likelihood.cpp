#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include "likelihood.h"
#include "model.h"


//////////////////////////////////
////  Prior parameter values  ////
//////////////////////////////////

#define NUMSPLINES 7
#define MAX_R 3
#define MAX_PREV 0.95

const double ancBiasPrior[2] = {0, 1};
const double logiota_unif_prior[] = {log(1e-13), log(0.0025)};
const double invGammaParameter = 0.001;

const double tau2_sample_rate = 0.5;
const double u0_sample_mu = 1.0;      // mean for prior sampling first spline coef 
const double u0_sample_sigma = 2.0;   // sd for sampling first spline coef


///////////////////////////////////////
////  Declare data for likelihood  ////
///////////////////////////////////////

const double ancLogitPrev[] = {-4.8719780, -4.2914736, -3.6969050, -3.1148216, -2.5022585,
                               -2.1492642, -1.8012415, -1.5856273, -1.2196389, -1.2425065,
                               -1.1254595, -1.1093076, -1.0201407, -0.9494274, -0.8712224,
                               -0.8377921, -0.8905323, -0.8760355, -0.8808581, -0.8760355,
                               -0.8377921, -0.8712224};
const double ancLogitVar[] = {0.0174009254, 0.0084540456, 0.0053964641, 0.0032544019, 0.0016378424,
                              0.0010016838, 0.0009123278, 0.0016000000, 0.0020186588, 0.0011394100,
                              0.0009205834, 0.0011694585, 0.0007565123, 0.0007092568, 0.0006018407,
                              0.0006458887, 0.0003913861, 0.0003867049, 0.0003882484, 0.0003398774,
                              0.0003295350, 0.0003385354};
const size_t ancIdx[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                         30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                         40, 41}; // note -1 vs R version to change from 1 based to 0 based indexing
const size_t numANCYears = 22;

const double hhsurvLogitPrev[] = {-1.688296, -1.643422, -1.592731, -1.463058};
const double hhsurvLogitVar[] = {0.005099825, 0.002752585, 0.002773331, 0.002178224};
const size_t hhsurvIdx[] = {32, 35, 38, 42}; // note -1 vs R version to change from 1 based to 0 based indexing
const size_t numHHSurvYears = 4;

/////////////////////
////  Functions  ////
/////////////////////

double logit(const double x) { return log(x / (1.0 - x)); };

double log_invgamma_pdf(const double x, const double shape, const double rate)
{
  return (shape * log(rate) - gsl_sf_lngamma(shape) - (shape + 1) * log(x) - rate / x);
}

double ll(struct modprev * out)
{

  // ll for ANC data
  double ll_S2 = 0.0, ll_dbar = 0.0, ll_d2bar = 0.0, ll_anc = 0.0;
  for(size_t i = 0; i < numANCYears; i++){
    double p_i = out->ANCprev[ancIdx[i]];
    if(isnan(p_i) || p_i <= 0.0 || p_i > MAX_PREV)
      return -INFINITY;
    double logitPi = logit(p_i);
    ll_S2 += 1.0/ancLogitVar[i];
    ll_dbar += (ancLogitPrev[i] - logitPi)/ancLogitVar[i];
    ll_d2bar += pow(ancLogitPrev[i] - logitPi, 2.0)/ancLogitVar[i];
    ll_anc += log(ancLogitVar[i]);
  }
  ll_S2 = 1.0/ll_S2;
  ll_dbar *= ll_S2;
  ll_d2bar *= ll_S2;

  ll_anc = log(ll_S2)/2.0 - ll_anc/2.0 + (pow(ll_dbar, 2.0) - ll_d2bar)/ll_S2 + log(gsl_cdf_gaussian_P(ancBiasPrior[1] - ll_dbar, sqrt(ll_S2)) - gsl_cdf_gaussian_P(ancBiasPrior[0] - ll_dbar, sqrt(ll_S2)));

  // ll for HSRC survey data
  double ll_hhsurv = 0.0;
  for(size_t i = 0; i < numHHSurvYears; i++){
    double p_i = out->a15to49prev[hhsurvIdx[i]];
    if(isnan(p_i) || p_i <= 0.0 || p_i > MAX_PREV)
      return -INFINITY;
    ll_hhsurv += -log(hhsurvLogitVar[i])/2.0 - pow(hhsurvLogitPrev[i] - logit(p_i), 2.0)/(2.0*hhsurvLogitVar[i]);
  }

  return ll_anc + ll_hhsurv;
}

double likelihood(const gsl_vector * theta)
{

  double iota = exp(gsl_vector_get(theta, 0));

  double * rVec = fnGenRVec(gsl_vector_const_ptr(theta, 2), NUMSPLINES, MAX_R); // mallocs rVec if successful
  if(rVec == NULL)
    return 0.0;
  
  struct modprev out;
  fnSpectrumPrev(iota, rVec, &out);
  fnFreeRVec(rVec);

  return(exp(ll(&out)));
}

double prior(const gsl_vector * theta)
{
  double tau2 = exp(gsl_vector_get(theta, 1));
  double tau = sqrt(tau2);
  
  double lprior = log(gsl_ran_flat_pdf(gsl_vector_get(theta, 0), 
				       logiota_unif_prior[0], 
				       logiota_unif_prior[1]));

  lprior += log_invgamma_pdf(tau2, invGammaParameter, invGammaParameter);

  for(size_t i = 4; i < (2+NUMSPLINES); i++)
    lprior += log(gsl_ran_gaussian_pdf(gsl_vector_get(theta, i), tau));
  return exp(lprior);
}

void sample_prior(gsl_rng * r, size_t numSamples, gsl_matrix * storeSamples)
{

  for(size_t i = 0; i < numSamples; i++){
    gsl_matrix_set(storeSamples, i, 0, 
		   gsl_ran_flat(r, logiota_unif_prior[0], logiota_unif_prior[1]));
    double tau2 = gsl_ran_exponential(r, 1.0/tau2_sample_rate);
    gsl_matrix_set(storeSamples, i, 1, log(tau2));
    double tau = sqrt(tau2);
    gsl_matrix_set(storeSamples, i, 2, 
		   gsl_ran_gaussian(r, u0_sample_sigma) + u0_sample_mu);
    for(size_t j = 3; j < (2+NUMSPLINES); j++)
      gsl_matrix_set(storeSamples, i, j, gsl_ran_gaussian(r, tau));
  }

  return;
}
