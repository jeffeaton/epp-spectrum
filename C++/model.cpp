#include <stddef.h>
#include <math.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>

#include "states.h"
#include "model.h"
#include "parameters.h"

#include <stdio.h>

/////////////////////////////
////  Declare functions  ////
/////////////////////////////

states rk4(states y, double dt, const struct parameters *);
states euler(states y, const double dt, const struct parameters *);
states grad(states y, const struct parameters *);
void RecordFullOutput(states * y, const size_t outIdx, const size_t numOutDates, double * Xout);
void PrepareHIV(states * y);
void PrepareART(states * y);


////////////////////////////
////  Define functions  ////
////////////////////////////

void fnSpectrum(const double iota, const double * rVec, const size_t numOutDates, double * Xout)
{
  // set parameters
  struct parameters param;

  // initialise population
  states current;
  current.hiv_idx = 1;
  current.art_idx = 1;
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < current.hiv_idx; m++)
        for(size_t u = 0; u < current.art_idx; u++)
          current.X[g][a][m][u] = 0.0;

  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      current.X[g][a][0][0] = init_pop[g][a];

  // simulate the model
  for(size_t ts = 0; ts < PROJ_STEPS; ts++){

    param.year_idx = (size_t) ts*dt;
    param.r = rVec[ts];

    // determine desired number on ART
    double frac_yr = fmod(ts*dt, 1.0) + dt;
    param.art_ts = frac_yr * artnum_15plus[param.year_idx] + (1.0 - frac_yr) * ((param.year_idx==0)?0.0:artnum_15plus[param.year_idx-1]);

    if(ts == ts0){
      param.iota = iota;
      PrepareHIV(&current);
    } else {
      param.iota = 0;
    }

    if(param.art_ts > 0 && current.art_idx != TS){
      PrepareART(&current);
    }

    current = euler(current, dt, &param);

    // record the outputs (midyear)
    if(ts % ((size_t) (1.0/dt)) + 1 == (1.0/(2*dt))){
      size_t out_idx = (size_t) ts * dt;
      RecordFullOutput(&current, out_idx, numOutDates, Xout);
    }

  }

  return;
}


/* numerical solvers */
states rk4(states y, double dt, const struct parameters * param)
{
  states k[4], out;
  k[0] = grad(y, param);
  k[1] = grad(y + k[0]*(.5*dt), param);
  k[2] = grad(y + k[1]*(.5*dt), param);
  k[3] = grad(y + k[2]*dt, param);

  out = y + (k[0] + k[1]*2 + k[2]*2 + k[3])*(dt/6);

  return out;
}

states euler(states y, const double dt, const struct parameters * param)
{
  states out = y + grad(y, param)*dt;
  return out;
}


/* the model gradient (differential equations) */
states grad(const states y, const struct parameters * param)
{
  states out;  // initialise states object to store gradient
  out.hiv_idx = y.hiv_idx;
  out.art_idx = y.art_idx;

  // initialise
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < y.hiv_idx; m++)
        for(size_t u = 0; u < y.art_idx; u++)
          out.X[g][a][m][u] = 0;

  // ageing
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG-1; a++)
      for(size_t m = 0; m < y.hiv_idx; m++)
        for(size_t u = 0; u < y.art_idx; u++){
          out.X[g][a][m][u] -= (1.0/AG_SPAN) * y.X[g][a][m][u];
          out.X[g][a+1][m][u] += (1.0/AG_SPAN) * y.X[g][a][m][u];
        }

  // natural mortality
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < y.hiv_idx; m++)
        for(size_t u = 0; u < y.art_idx; u++)
          out.X[g][a][m][u] -= mx[param->year_idx][g][a] * y.X[g][a][m][u];

  // fertility
  double births_by_age[AG_FERT], hivp_by_age[AG_FERT];
  for(size_t a = 0; a < AG_FERT; a++){
    births_by_age[a] = 0;
    hivp_by_age[a] = 0;
    for(size_t m = 0; m < y.hiv_idx; m++)
      for(size_t u = 0; u < y.art_idx; u++){
        births_by_age[a] += asfr[param->year_idx][a] * y.X[FEMALE][IDX_FERT + a][m][u];
        if(m >= 1)
          hivp_by_age[a] += y.X[FEMALE][IDX_FERT + a][m][u];
      }
  }

  double frac_hivp_births[AG_FERT];
  for(size_t a = 0; a < AG_FERT; a++)
    frac_hivp_births[a] = vert_trans*(1.0 - y.X[FEMALE][IDX_FERT + a][0][0]/(fert_rat[a]*hivp_by_age[a] + y.X[FEMALE][IDX_FERT + a][0][0]));

  for(size_t a = 0; a < AG_FERT; a++){
    out.X[MALE][0][0][0] += births_by_age[a] * (1.0 - frac_hivp_births[a]) * srb[param->year_idx] / (srb[param->year_idx] + 1.0);
    out.X[FEMALE][0][0][0] += births_by_age[a] * (1.0 - frac_hivp_births[a]) * 1.0 / (srb[param->year_idx] + 1.0);
    for(size_t m = 1; m < y.hiv_idx; m++){
      out.X[MALE][0][m][0] += births_by_age[a] * frac_hivp_births[a] * cd4_initdist[MALE][0][m-1] * srb[param->year_idx] / (srb[param->year_idx] + 1.0);
      out.X[FEMALE][0][m][0] += births_by_age[a] * frac_hivp_births[a] * cd4_initdist[FEMALE][0][m-1] * 1.0 / (srb[param->year_idx] + 1.0);
    }
  }


  if(y.hiv_idx > 1){

    // incidence
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

    double inc_rate_15to49 = param->r * ((Xhivp_noart + relinfect_art * Xart)/Xtot + param->iota);

    double inc_rr[NG][AG]; // incidence rate ratio by age and sex
    for(size_t a = 0; a < AG; a++){
      inc_rr[MALE][a] = inc_agerat[param->year_idx][MALE][a];
      inc_rr[FEMALE][a] = inc_sexrat[param->year_idx] * inc_agerat[param->year_idx][FEMALE][a];
    }

    double Xhivn_incrr = 0;
    for(size_t g = 0; g < NG; g++)
      for(size_t a = IDX_15TO49; a < IDX_15TO49+AG_15TO49; a++){
        Xhivn_incrr += inc_rr[g][a] * y.X[g][a][0][0];
      }

    double age_inc[NG][AG];
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++){
        age_inc[g][a] = inc_rate_15to49 * inc_rr[g][a] * Xhivn / Xhivn_incrr;
        out.X[g][a][0][0] -= age_inc[g][a] * y.X[g][a][0][0];
        for(size_t m = 1; m < DS; m++)
          out.X[g][a][m][0] += age_inc[g][a] * cd4_initdist[g][a][m-1] * y.X[g][a][0][0];
      }

    // disease progression and mortality

    // CD4 stage progression
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
        for(size_t m = 1; m < DS-1; m++){
          out.X[g][a][m][0] -= cd4_prog[g][a][m-1] * y.X[g][a][m][0];
          out.X[g][a][m+1][0] += cd4_prog[g][a][m-1] * y.X[g][a][m][0];
        }

    // HIV and ART mortality
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
        for(size_t m = 1; m < y.hiv_idx; m++)
          for(size_t u = 0; u < y.art_idx; u++)
            out.X[g][a][m][u] -= cd4_art_mort[g][a][m-1][u] * y.X[g][a][m][u];

    if(y.art_idx > 1){

      // ART duration stage progression
      for(size_t g = 0; g < NG; g++)
        for(size_t a = 0; a < AG; a++)
          for(size_t m = 1; m < y.hiv_idx; m++)
            for(size_t u = 1; u < TS - 1; u++){
              out.X[g][a][m][u] -= art_prog[u-1] * y.X[g][a][m][u];
              out.X[g][a][m][u+1] += art_prog[u-1] * y.X[g][a][m][u];
            }

      // ART initiation
      double Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig_15plus = 0.0, grad_art_cx = 0.0;
      for(size_t g = 0; g < NG; g++)
        for(size_t a = IDX_15PLUS; a < AG; a++)
          for(size_t m = 1; m < y.hiv_idx; m++){
            if(m >= artelig_idx[param->year_idx]){
              Xartelig_15plus += y.X[g][a][m][0];
              expect_mort_artelig_15plus += cd4_art_mort[g][a][m-1][0] * y.X[g][a][m][0];
            }
            for(size_t u = 1; u < TS; u++){
              Xart_15plus += y.X[g][a][m][u];
              grad_art_cx += out.X[g][a][m][u];
            }
          }

      double art_15plus_anninit = (param->art_ts  - Xart_15plus) / dt - grad_art_cx; // desired number to initiate per yr (elig * rate)

      for(size_t g = 0; g < NG; g++)
        for(size_t a = IDX_15PLUS; a < AG; a++)
          for(size_t m = artelig_idx[param->year_idx]; m < DS; m++){
            double art_initrate = art_15plus_anninit * 0.5 * (1.0/Xartelig_15plus + cd4_art_mort[g][a][m-1][0] / expect_mort_artelig_15plus);
            if(art_initrate > 1.0/dt)
              art_initrate = 1.0/dt;
            out.X[g][a][m][0] -= art_initrate * y.X[g][a][m][0];
            out.X[g][a][m][1] += art_initrate * y.X[g][a][m][0];
          }

    } // if(y.art_idx > 1)

  } // if(y.hiv_idx > 1)

  return out;
}


void RecordFullOutput(states * y, const size_t outIdx, const size_t numOutDates, double * Xout)
{
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++){
      for(size_t m = 0; m < y->hiv_idx; m++){
        for(size_t u = 0; u < y->art_idx; u++)
          Xout[outIdx + (g + NG*(a + AG*(m + DS*u)))*numOutDates] = y->X[g][a][m][u];
        for(size_t u = y->art_idx; u < TS; u++)
          Xout[outIdx + (g + NG*(a + AG*(m + DS*u)))*numOutDates] = 0.0;
      }
      for(size_t m = y->hiv_idx; m < DS; m++)
        for(size_t u = 0; u < TS; u++)
          Xout[outIdx + (g + NG*(a + AG*(m + DS*u)))*numOutDates] = 0.0;
    }

  return;
}

void PrepareHIV(states * y)
{
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = y->hiv_idx; m < DS; m++)
        y->X[g][a][m][0] = 0.0;
  y->hiv_idx = DS;
}

void PrepareART(states * y)
{
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < DS; m++)
        for(size_t u = y->art_idx; u < TS; u++)
          y->X[g][a][m][u] = 0.0;
  y->art_idx = TS;
}


/////////////////////////////////////////
////  Models for force of infection  ////
/////////////////////////////////////////

#define SPLINE_ORDER 4
#define NUM_SPLINES

void fnBSpline(const double * u, const size_t numSplines, const size_t numSplineSteps, double * rVec)
{

  const size_t nbreaks = numSplines - SPLINE_ORDER + 2;
  double dk = ((double) (numSplineSteps - 1))/(nbreaks - 1); // distance between knots

  // declare workspace
  gsl_bspline_workspace *w = gsl_bspline_alloc(SPLINE_ORDER, nbreaks);
  gsl_vector *x = gsl_vector_alloc(numSplines);
  gsl_matrix *designMat = gsl_matrix_alloc(numSplineSteps, numSplines);

  // set the knot locations
  for(size_t i = 0; i < numSplines + SPLINE_ORDER; i++)
    gsl_vector_set(w->knots, i, 1.0 + dk * ((double) i - (SPLINE_ORDER-1)));

  // evaluate design matrix at each step
  for(size_t i = 0; i < numSplineSteps; i++){
    gsl_bspline_eval(1+i, x, w); // 1+i to match 1-based indexing in Dan's R code
    gsl_matrix_set_row(designMat, i, x);
   }
  

  // determine spline coefficients
  gsl_vector *coefs = gsl_vector_alloc(numSplines);
  gsl_vector_set(coefs, 0, u[0]);
  gsl_vector_set(coefs, 1, u[1]);
  for(size_t i = 2; i < numSplines; i++)
    gsl_vector_set(coefs, i, 2*gsl_vector_get(coefs, i-1) - gsl_vector_get(coefs, i-2) + u[i]);


  gsl_vector_view vw_rVec = gsl_vector_view_array(rVec, numSplineSteps);
  gsl_blas_dgemv(CblasNoTrans, 1.0, designMat, coefs, 0.0, &(vw_rVec.vector));

  return;
}

double * fnGenRVec(const double * u, const size_t numSplines)
{
  double * rVec = (double *) calloc(PROJ_STEPS, sizeof(double));
  fnBSpline(u, numSplines, PROJ_STEPS - ts0, &rVec[ts0]);
  return rVec;
}

void fnFreeRVec(double * rVec)
{
  free(rVec);
  return;
}
