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

double ll(struct modprev * out);  // forward declare, defined in likelihood.cpp

void initParams(struct parameters * param, const double iota, const double * rVec, const double inc_sexrat[], const double inc_agerat[][NG][AG], const double rmat[AG][AG])
{
  param->iota = iota;

  for(size_t i = 0; i < PROJ_STEPS; i++){
    param->rVec[i] = rVec[i];
  }

  if(inc_sexrat != NULL && inc_agerat != NULL){
    param->incmod = INC_INCRR;

    for(size_t i = 0; i < PROJ_YEARS; i++)
      param->inc_sexrat[i] = inc_sexrat[i];

    for(size_t i = 0; i < PROJ_YEARS; i++)
      for(size_t j = 0; j < NG; j++)
        for(size_t k = 0; k < AG; k++)
          param->inc_agerat[i][j][k] = inc_agerat[i][j][k];
  }

  if(rmat){
    param->incmod = INC_RMAT;

    for(size_t i = 0; i < AG; i++)
      for(size_t j = 0; j < AG; j++)
        param->rmat[i][j] = rmat[i][j];
  }

  return;
}

int main()
{

  size_t numOutDates = NUM_OUTDATES;
  
  double * Xout = (double *) malloc(numOutDates * NG * AG * DS *TS * sizeof(double));

  struct parameters param_incrr;
  param_incrr.rVec = (double *) malloc(sizeof(double) * PROJ_STEPS);
  initParams(&param_incrr, iota_test, rVec_incrr_test, inc_sexrat_test, inc_agerat_test, NULL);

  fnSpectrum(&param_incrr, numOutDates, Xout);
  
  printf("\n");
  for(size_t a = 0; a < AG; a++)
    printf("%10f %10f\n", Xout[42 + (0 + NG*(a + AG*(0 + DS*0)))*numOutDates],
           Xout[42 + (1 + NG*(a + AG*(0 + DS*0)))*numOutDates]);
  
  free(Xout);
  
  struct modprev out;
  fnSpectrumPrev(&param_incrr, &out);
  
  printf("\nANC prev  15 to 49 prev\n");
  for(size_t i = 0; i < numOutDates; i++)
    printf("%8.4f  %13.4f\n", out.ANCprev[i], out.a15to49prev[i]);
  
  // printf("\nll = %f\n", ll(&out));

  free(param_incrr.rVec);
  
  // test age-specific force of infection model
  printf("\n*** Age-specific force of infection model ***\n");

  struct parameters param_rmat;
  param_rmat.rVec = (double *) malloc(sizeof(double) * PROJ_STEPS);
  initParams(&param_rmat, iota_test, rVec_rmat_test, NULL, NULL, rmat_test);

  fnSpectrumPrev(&param_rmat, &out);

  printf("\nANC prev  15 to 49 prev\n");
  for(size_t i = 0; i < numOutDates; i++)
    printf("%8.4f  %13.4f\n", out.ANCprev[i], out.a15to49prev[i]);
  
  // printf("\nll = %f\n", ll(&out));

  free(param_rmat.rVec);

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
