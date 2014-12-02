#include <stddef.h>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#include "states.h"
#include "parameters.h"
#include "mvndstpack.h"


void fnEffectiveAgePrev(const states y, const double relinfectART, double eff_ageprev[NG][AG])
{

  // calculate the HIV prevalence adjusted for relative infectiousness of ART by age / sex

  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++){
      double Xhivn = y.X[g][a][0][0];
      double Xhivp_noart = 0.0;
      double Xart = 0.0;
      for(size_t m = 1; m < y.hiv_idx; m++){
        Xhivp_noart += y.X[g][a][m][0];
        for(size_t u = 1; u < y.art_idx; u++){
          Xart += y.X[g][a][m][u];
        }
      }
      eff_ageprev[g][a] = (Xhivp_noart + relinfectART * Xart) / (Xhivp_noart + Xart + Xhivn);
    }

  return;
}



void fnAgeInc(const states y, const struct parameters * p, double age_inc[NG][AG])
{

  if(p->incmod == INC_INCRR){
    ////////////////////////////////////////////
    ////  fixed incidence rate ratio model  ////
    ////////////////////////////////////////////

    double Xhivn = 0.0, Xhivp_noart = 0.0, Xart = 0.0;
    for(size_t g = 0; g < NG; g++)
      for(size_t a = IDX_15TO49; a < IDX_15TO49+AG_15TO49; a++){
        Xhivn += y.X[g][a][0][0];
        for(size_t m = 1; m < y.hiv_idx; m++){
          Xhivp_noart += y.X[g][a][m][0];
          for(size_t u = 1; u < y.art_idx; u++){
            Xart += y.X[g][a][m][u];
          }
        }
      }
    double Xtot = Xhivn + Xhivp_noart + Xart;

    double inc_rate_15to49 = p->rvec[p->ts] * (Xhivp_noart + p->relinfect_art * Xart)/Xtot + p->iota_ts;

    // printf("ts %lu, Xhivn %f, Xhivp_noart %f, Xart %f, incrate %f\n", p->ts, Xhivn, Xhivp_noart, Xart, inc_rate_15to49);

    double inc_rr[NG][AG]; // incidence rate ratio by age and sex
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
	inc_rr[g][a] = p->agesex_incrr[p->ts][g][a];

    double Xhivn_incrr = 0;
    for(size_t g = 0; g < NG; g++)
      for(size_t a = IDX_15TO49; a < IDX_15TO49+AG_15TO49; a++){
        Xhivn_incrr += inc_rr[g][a] * y.X[g][a][0][0];
      }

    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
        age_inc[g][a] = inc_rate_15to49 * inc_rr[g][a] * Xhivn / Xhivn_incrr;

    /*
    // if(p->iota_ts > 0){
    if(p->year_idx == 6){
      printf("iota_ts %f\n", p->iota_ts);
      printf("year_idx %lu\n", p->year_idx);
      for(size_t i = 0; i < AG; i++)
	printf("%f %f\n", age_inc[0][i], age_inc[1][i]);
      // printf("%f %f\n", p->inc_agerat[MALE][i], p->inc_agerat[FEMALE][i]);
    }
    */
    

  } else if(p->incmod == INC_RMAT) {
    /////////////////////////////////////////////////
    ////  age-specific force of infection model  ////
    /////////////////////////////////////////////////

    double eff_ageprev[NG][AG];
    fnEffectiveAgePrev(y, p->relinfect_art, eff_ageprev);

    // add initial infection pulse
    if(p->iota_ts > 0){
      for(size_t i = 0; i < NG; i++)
        for(size_t j = 0; j < AG; j++)
          eff_ageprev[i][j] += p->iota_ts;
    }

    // initialise age_inc to 0
    for(size_t i = 0; i < NG; i++)
      for(size_t j = 0; j < AG; j++)
        age_inc[i][j] = 0.0;

    // matrix multiplication of eff_ageprev %*% Rmat
    for(size_t i = 0; i < AG; i++)  // i - susceptible age index
      for(size_t j = 0; j < AG; j++){  // j - infected age index
        age_inc[MALE][i] +=  eff_ageprev[FEMALE][j] * p->rmat[i][j];
        age_inc[FEMALE][i] +=  eff_ageprev[MALE][j] * p->rmat[j][i];
      }

    for(size_t i = 0; i < NG; i++)
      for(size_t j = 0; j < AG; j++)
        age_inc[i][j] *= p->rvec[p->ts];
  }

  return;
}


void createRmat(double m_mean, double m_sd, double f_mean, double f_sd, double corr, double rmat[AG][AG])
{

  if(m_mean <= 0 || m_sd <= 0.1 || f_mean <= 0 || f_sd <= 0.1){
    for(size_t i = 0; i < AG; i++)
      for(size_t j = 0; j < AG; j++)
	rmat[i][j] = 0.0;
    return;
  }

  // gamma distribution parameters (GSL -- shape/scale parameterisation)
  double m_shape = m_mean*m_mean / (m_sd * m_sd);
  double m_scale = m_sd * m_sd / m_mean;
  double f_shape = f_mean*f_mean / (f_sd * f_sd);
  double f_scale = f_sd * f_sd / f_mean;

  double m_quant[AG-3], f_quant[AG-3];
  for(size_t i = 0; i < AG-3; i++){

    m_quant[i] = gsl_cdf_gamma_P(AG_SPAN*i, m_shape, m_scale);
    f_quant[i] = gsl_cdf_gamma_P(AG_SPAN*i, f_shape, f_scale);

    m_quant[i] = gsl_cdf_ugaussian_Pinv(m_quant[i]);
    f_quant[i] = gsl_cdf_ugaussian_Pinv(f_quant[i]);
  }

  // arguments for mvndst
  double lower[2], upper[2];
  int n[1] = {2}, infin[2], maxpts[1] = {25000};
  double abseps[1] = {0.001}, releps[1] = {0.0};

  // pass for outputs of mvndst
  double error[1], value[1];
  int inform[1];

  for(size_t i = 0; i < AG; i++)
    for(size_t j = 0; j < AG; j++){

      if(i < IDX_15PLUS || j < IDX_15PLUS)
        rmat[i][j] = 0.0;
      else {
        size_t midx = i - IDX_15PLUS;
        size_t fidx = j - IDX_15PLUS;

        if(m_quant[midx] == GSL_NEGINF){
          infin[0] = 0;
          lower[0] = 0.0;
          upper[0] = m_quant[midx+1];
        } else if(i == AG-1){
          infin[0] = 1;
          lower[0] = m_quant[midx];
          upper[0] = 0.0;
        } else {
          infin[0] = 2;
          lower[0] = m_quant[midx];
          upper[0] = m_quant[midx+1];
        }

        if(f_quant[fidx] == GSL_NEGINF){
          infin[1] = 0;
          lower[1] = 0.0;
          upper[1] = f_quant[fidx+1];
        } else if(j == AG-1){
          infin[1] = 1;
          lower[1] = f_quant[fidx];
          upper[1] = 0.0;
        } else {
          infin[1] = 2;
          lower[1] = f_quant[fidx];
          upper[1] = f_quant[fidx+1];
        }

        mvndst_(n, lower, upper, infin, &corr, maxpts, abseps, releps, error, value, inform);
        rmat[i][j] = *value;
      }
    }

  return;
}
