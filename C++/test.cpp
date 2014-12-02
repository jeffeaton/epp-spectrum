#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include "states.h"
#include "model.h"
#include "parameters.h"
#include "likelihood.h"

#include "test.h" // declare test parameter values
#include "test_parameters.h"

double ll(struct modprev * out);  // forward declare, defined in likelihood.cpp

// !!!! TODO: temporary, delete PROJ_STEPS

void initParams(struct parameters * param, const double iota, const double * rvec, const double inc_sexrat[], const double inc_agerat[][NG][AG], const double rmat[AG][AG])
{
    
  param->dt = input_dt;
  param->n_proj_steps = PROJ_STEPS;

  for(size_t j = 0; j < NG; j++)
    for(size_t k = 0; k < AG; k++)
      param->init_pop[j][k] = input_init_pop[j][k];


  for(size_t i = 0; i < PROJ_STEPS; i++){
    size_t year_idx = (size_t) i*param->dt;

    param->srb[i] = input_srb[year_idx];
    param->artelig_idx[i] = input_artelig_idx[year_idx];
    param->rvec[i] = rvec[i];
      
    for(size_t k = 0; k < AG_FERT; k++)
      param->asfr[i][k] = input_asfr[year_idx][k];

    for(size_t j = 0; j < NG; j++)
      for(size_t k = 0; k < AG; k++)
	param->mx[i][j][k] = input_mx[year_idx][j][k];

    // determine desired number on ART
    double frac_yr = fmod(i*param->dt, 1.0) + param->dt;
    param->artnum_15plus[i] = frac_yr * input_artnum_15plus[year_idx] + (1.0 - frac_yr) * ((year_idx==0)?0.0:input_artnum_15plus[year_idx-1]);
    
  }

  for(size_t k = 0; k < AG_FERT; k++)
    param->fert_rat[k] = input_fert_rat[k];


  // HIV incidence model parameters
  param->iota = iota;
  param->ts_epi_start = input_ts0;
  param->relinfect_art = input_relinfect_art;
  param->vert_trans = input_vert_trans;
  
  if(inc_sexrat != NULL && inc_agerat != NULL){
    param->incmod = INC_INCRR;

    for(size_t i = 0; i < PROJ_STEPS; i++)
      for(size_t k = 0; k < AG; k++){
	param->agesex_incrr[i][MALE][k] = inc_agerat[0][MALE][k];
	param->agesex_incrr[i][FEMALE][k] = inc_sexrat[0] * inc_agerat[0][FEMALE][k];
      }
  }

  if(rmat){
    param->incmod = INC_RMAT;

    for(size_t i = 0; i < AG; i++)
      for(size_t j = 0; j < AG; j++)
        param->rmat[i][j] = rmat[i][j];
  }


  // HIV progression model parameters

  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < (DS-2); m++)
	param->cd4_prog[g][a][m] = input_cd4_prog[g][a][m];

  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < (DS-1); m++)
	for(size_t u = 0; u < (DS-1); u++)
	  param->cd4_art_mort[g][a][m][u] = input_cd4_art_mort[g][a][m][u];


  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < (DS-1); m++)
	param->cd4_initdist[g][a][m] = input_cd4_initdist[g][a][m];
  

  return;
}

int main()
{

  size_t numOutDates = NUM_OUTDATES;
  
  double * Xout = (double *) malloc(numOutDates * NG * AG * DS *TS * sizeof(double));

  struct parameters *param_incrr = alloc_parameters(PROJ_STEPS);
  initParams(param_incrr, iota_test, rVec_incrr_test, inc_sexrat_test, inc_agerat_test, NULL);

  fnSpectrum(param_incrr, numOutDates, Xout);
  
  printf("\n");
  for(size_t a = 0; a < AG; a++)
    printf("%10f %10f\n", Xout[42 + (0 + NG*(a + AG*(0 + DS*0)))*numOutDates],
           Xout[42 + (1 + NG*(a + AG*(0 + DS*0)))*numOutDates]);
  
  free(Xout);
  
  struct modprev out;
  fnSpectrumPrev(param_incrr, &out);
  
  printf("\nANC prev  15 to 49 prev\n");
  for(size_t i = 0; i < numOutDates; i++)
    printf("%8.4f  %13.4f\n", out.ANCprev[i], out.a15to49prev[i]);
  
  // printf("\nll = %f\n", ll(&out));

  free_parameters(param_incrr);
  
  // test age-specific force of infection model
  printf("\n*** Age-specific force of infection model ***\n");

  struct parameters *param_rmat = alloc_parameters(PROJ_STEPS);
  initParams(param_rmat, iota_test, rVec_rmat_test, NULL, NULL, rmat_test);

  fnSpectrumPrev(param_rmat, &out);

  printf("\nANC prev  15 to 49 prev\n");
  for(size_t i = 0; i < numOutDates; i++)
    printf("%8.4f  %13.4f\n", out.ANCprev[i], out.a15to49prev[i]);
  
  // printf("\nll = %f\n", ll(&out));

  free_parameters(param_rmat);

  /*
  const double theta[] = {-27.2776051, 0.2286067, 1.5936094, 2.0196487, -0.4601538, -1.6417664, 1.4846658, 0.1725820, 0.3623600};
  // const double theta[] = {-2.427689764762438e+01, -8.060906058729440e-02, 1.962444061070041e+00, 7.898717801500630e-01, 7.133904109969543e-01, 1.734195483892624e-01, 2.934792020244043e-01, 1.134129669092168e+00, -2.662140188992846e-01};

  gsl_vector_const_view vw_theta = gsl_vector_const_view_array(theta, 9);


  printf("\n::Example with theta vector::\n");
  printf("log(prior): %f\n", log(prior(&(vw_theta.vector))));
  printf("log(likelihood): %f\n", log(likelihood(&(vw_theta.vector))));
  */
    return 0;
}
