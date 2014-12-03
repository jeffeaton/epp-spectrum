#include <stddef.h>
#include <math.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>

#include "states.h"
#include "model.h"
#include "incidence.h"
#include "parameters.h"

#include <stdio.h>

/////////////////////////////
////  Declare functions  ////
/////////////////////////////

void InitialisePopulation(states * y, const parameters * param);
void SetParamTS(const size_t ts, parameters * param);
void PrepareHIV(states * y);
void PrepareART(states * y);

states rk4(states y, double dt, const struct parameters *);
states euler(states y, const double dt, const struct parameters *);
states grad(states y, const struct parameters *);

// Output functions
void RecordFullOutput(states * y, const size_t outIdx, const size_t numOutDates, double * Xout);
double fn15to49prev(states * y, parameters * param);
double fnANCprev(states * y, parameters * param);



// HARD CODED: ART stage progression rate
const double ART_STAGE_PROG = 2.0;


////////////////////////////
////  Define functions  ////
////////////////////////////

void fnSpectrum(struct parameters * param, const size_t numOutDates, double * Xout)
{
  // set parameters
  // struct parameters param;

  double dt = param->dt;

  // initialise population
  states current;
  InitialisePopulation(&current, param);

  // simulate the model
  for(size_t ts = 0; ts < param->n_proj_steps; ts++){

    // record the outputs (midyear)
    if(fmod(param->proj_steps[ts], 1.0) == dt * ((size_t) (1.0/(2*dt)))){
      size_t out_idx = (size_t) ts * dt;
      RecordFullOutput(&current, out_idx, numOutDates, Xout);
    }

    SetParamTS(ts, param);

    if(ts == param->ts_epi_start)
      PrepareHIV(&current);

    if(param->artnum_15plus[param->ts][FEMALE] > 0 && current.art_idx != TS)
      PrepareART(&current);

    current = euler(current, dt, param);

  }

  return;
}

void fnSpectrumPrev(struct parameters * param, struct modprev * out)
{
  // set parameters
  // struct parameters param;

  double dt = param->dt;

  // initialise population
  states current;
  InitialisePopulation(&current, param);

  // simulate the model
  for(size_t ts = 0; ts < param->n_proj_steps; ts++){

    SetParamTS(ts, param);

    if(ts == param->ts_epi_start)
      PrepareHIV(&current);

    if(param->artnum_15plus[param->ts][FEMALE] > 0 && current.art_idx != TS)
      PrepareART(&current);

    current = euler(current, dt, param);

    // record the outputs (midyear)
    if(ts % ((size_t) (1.0/dt)) + 1 == (1.0/(2*dt))){
      size_t out_idx = (size_t) ts * dt;
      out->ANCprev[out_idx] = fnANCprev(&current, param);
      out->a15to49prev[out_idx] = fn15to49prev(&current, param);
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


///////////////////////////////////////////////////////
////  the model gradient (differential equations)  ////
///////////////////////////////////////////////////////

void calcFertility(const states * y, const struct parameters * param, double * num_births, double * frac_hivp_moth)
{

  *num_births = 0.0;
  *frac_hivp_moth = 0.0;

  for(size_t a = 0; a < AG_FERT; a++){
    double births_ag = 0.0, weight_hivp_ag = 0.0;
    for(size_t m = 0; m < y->hiv_idx; m++)
      for(size_t u = 0; u < y->art_idx; u++){
        births_ag += param->asfr[param->ts][a] * y->X[FEMALE][IDX_FERT + a][m][u];
        if(m >= 1){
	  if(u < (TS-1))  // HARD CODED: If on ART < 1 year, same fertility as not on ART, constant subfertility for ART > 1 year
	    weight_hivp_ag += y->X[FEMALE][IDX_FERT + a][m][u] * param->age_fertrat[a] * param->stage_fertrat[m-1];
	  else
	    weight_hivp_ag += y->X[FEMALE][IDX_FERT + a][m][u] * param->age_fertrat[a] * param->art_fertrat;
	}
      }
    double frac_hivp_moth_ag = 1.0 - y->X[FEMALE][IDX_FERT + a][0][0]/(weight_hivp_ag + y->X[FEMALE][IDX_FERT + a][0][0]);
    *frac_hivp_moth += births_ag*frac_hivp_moth_ag;
    *num_births += births_ag;
  }

  *frac_hivp_moth /= *num_births;

  return;
}

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
          out.X[g][a][m][u] -= param->mx[param->ts][g][a] * y.X[g][a][m][u];

  // fertility
  double num_births, frac_hivp_moth;
  calcFertility(&y, param, &num_births, &frac_hivp_moth);

  out.X[MALE][0][0][0] += num_births * (1.0 - param->vert_trans*frac_hivp_moth) * param->srb[param->ts] / (param->srb[param->ts] + 1.0);
  out.X[FEMALE][0][0][0] += num_births * (1.0 - param->vert_trans*frac_hivp_moth) * 1.0 / (param->srb[param->ts] + 1.0);
  for(size_t m = 1; m < y.hiv_idx; m++){
    out.X[MALE][0][m][0] += num_births * param->vert_trans * frac_hivp_moth * param->cd4_initdist[MALE][0][m-1] * param->srb[param->ts] / (param->srb[param->ts] + 1.0);
    out.X[FEMALE][0][m][0] += num_births * param->vert_trans * frac_hivp_moth * param->cd4_initdist[FEMALE][0][m-1] * 1.0 / (param->srb[param->ts] + 1.0);
  }

  if(y.hiv_idx > 1){

    // incidence
    double age_inc[NG][AG];
    fnAgeInc(y, param, age_inc);

    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++){
        out.X[g][a][0][0] -= age_inc[g][a] * y.X[g][a][0][0];
        for(size_t m = 1; m < DS; m++)
          out.X[g][a][m][0] += age_inc[g][a] * param->cd4_initdist[g][a][m-1] * y.X[g][a][0][0];
      }

    // disease progression and mortality

    // CD4 stage progression
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
        for(size_t m = 1; m < DS-1; m++){
          out.X[g][a][m][0] -= param->cd4_prog[g][a][m-1] * y.X[g][a][m][0];
          out.X[g][a][m+1][0] += param->cd4_prog[g][a][m-1] * y.X[g][a][m][0];
        }

    // HIV and ART mortality
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
        for(size_t m = 1; m < y.hiv_idx; m++)
          for(size_t u = 0; u < y.art_idx; u++)
            out.X[g][a][m][u] -= param->cd4_art_mort[g][a][m-1][u] * y.X[g][a][m][u];

    if(y.art_idx > 1){

      // ART duration stage progression
      for(size_t g = 0; g < NG; g++)
        for(size_t a = 0; a < AG; a++)
          for(size_t m = 1; m < y.hiv_idx; m++)
            for(size_t u = 1; u < TS - 1; u++){
              out.X[g][a][m][u] -= ART_STAGE_PROG * y.X[g][a][m][u];
              out.X[g][a][m][u+1] += ART_STAGE_PROG * y.X[g][a][m][u];
            }

      // ART initiation
      // (calculated separately for men and women)

      for(size_t g = 0; g < NG; g++){
	double Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig_15plus = 0.0, grad_art_cx = 0.0;
        for(size_t a = IDX_15PLUS; a < AG; a++)
          for(size_t m = 1; m < y.hiv_idx; m++){
            if(m >= param->artelig_idx[param->ts]){
              Xartelig_15plus += y.X[g][a][m][0];
              expect_mort_artelig_15plus += param->cd4_art_mort[g][a][m-1][0] * y.X[g][a][m][0];
            }
            for(size_t u = 1; u < TS; u++){
              Xart_15plus += y.X[g][a][m][u];
              grad_art_cx += out.X[g][a][m][u];
            }
          }
	
	double art_15plus_anninit = (param->artnum_15plus[param->ts][g]  - Xart_15plus) / param->dt - grad_art_cx; // desired number to initiate per yr (elig * rate)
	
	for(size_t a = IDX_15PLUS; a < AG; a++)
	  for(size_t m = param->artelig_idx[param->ts]; m < DS; m++){
	    double art_initrate = art_15plus_anninit * 0.5 * (1.0/Xartelig_15plus + param->cd4_art_mort[g][a][m-1][0] / expect_mort_artelig_15plus);
	    if(art_initrate > 1.0/param->dt)
	      art_initrate = 1.0/param->dt;
	    out.X[g][a][m][0] -= art_initrate * y.X[g][a][m][0];
	    out.X[g][a][m][1] += art_initrate * y.X[g][a][m][0];
	  }
      }
      
    } // if(y.art_idx > 1)

  } // if(y.hiv_idx > 1)

  return out;
}

void InitialisePopulation(states * y, const parameters * param)
{
  y->hiv_idx = 1;
  y->art_idx = 1;
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < y->hiv_idx; m++)
        for(size_t u = 0; u < y->art_idx; u++)
          y->X[g][a][m][u] = 0.0;
  
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      y->X[g][a][0][0] = param->init_pop[g][a];

  return;
}

void SetParamTS(const size_t ts, parameters * param)
{
  param->ts = ts;
    
  if(ts == param->ts_epi_start)
      param->iota_ts = param->iota;
  else 
    param->iota_ts = 0;

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


/////////////////////////
////  Model outputs  ////
/////////////////////////

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

 double fn15to49prev(states * y, parameters * param)
 {

    double Xhivn = 0.0, Xhivp_noart = 0.0, Xart = 0.0;
    for(size_t g = 0; g < NG; g++)
      for(size_t a = IDX_15TO49; a < IDX_15TO49+AG_15TO49; a++){
        Xhivn += y->X[g][a][0][0];
        for(size_t m = 1; m < y->hiv_idx; m++){
          Xhivp_noart += y->X[g][a][m][0];
          for(size_t u = 1; u < y->art_idx; u++){
            Xart += y->X[g][a][m][u];
          }
        }
      }
    double Xtot = Xhivn + Xhivp_noart + Xart;

    return 1.0 - Xhivn / Xtot;
 }

 double fnANCprev(states * y, parameters * param)
 {

  double num_births, frac_hivp_moth;
  calcFertility(y, param, &num_births, &frac_hivp_moth);

  return frac_hivp_moth;
 }


/////////////////////////////////////////
////  Models for force of infection  ////
/////////////////////////////////////////

#define SPLINE_ORDER 4

void fnBSpline(const double * u, const size_t numSplines, const size_t numSplineSteps, double * rvec)
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
  for(size_t i = 2; i < numSplines; i++){
    gsl_vector_set(coefs, i, 2*gsl_vector_get(coefs, i-1) - gsl_vector_get(coefs, i-2) + u[i]);
  }
  
  gsl_vector_view vw_rvec = gsl_vector_view_array(rvec, numSplineSteps);
  gsl_blas_dgemv(CblasNoTrans, 1.0, designMat, coefs, 0.0, &vw_rvec.vector);
  
  gsl_vector_free(coefs);
  gsl_bspline_free(w);
  gsl_vector_free(x);
  gsl_matrix_free(designMat);
    
  return;
}

/* 
   fnGenRVec allocates memory for r-vector, calls fnBSpline to generate spline
   (appropriately adusted for t0).

   Returns pointer to rvec if all values of rvec are between [0, maxR), else frees
   the memory allocated for rvec and returns NULL
*/

/*
double * fnGenRVec(const double * u, const size_t numSplines, const double maxR)
{
  double * rVec = (double *) calloc(PROJ_STEPS, sizeof(double));
  fnBSpline(u, numSplines, PROJ_STEPS - ts0, &rVec[ts0]);
  
  for(size_t i = ts0; i < PROJ_STEPS; i++)
    if(rVec[i] < 0 || rVec[i] >= maxR){
      free(rVec);
      return NULL;
    }
  
  return rVec;
}

void fnFreeRVec(double * rVec)
{
  free(rVec);
  return;
}
*/
