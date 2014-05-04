#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include "likelihood.h"

// suggested call: $ ./a.out 100 3168947

int main(int argc, char *argv[])
{

  size_t NumParam = 9;
  size_t InitSamples = (unsigned int) atof(argv[1]);
  unsigned int rng_seed = (unsigned int) atof(argv[2]);

  // Declare and configure GSL RNG
  gsl_rng * rng;
  const gsl_rng_type * T;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  gsl_rng_set(rng, rng_seed);

  gsl_matrix * Xmat = gsl_matrix_alloc(InitSamples, NumParam);

  sample_prior(rng, InitSamples, Xmat);

  // for(size_t i = 0; i < InitSamples; i++){
  //   for(size_t j = 0; j < NumParam; j++)
  //     printf("%10.5f ", gsl_matrix_get(Xmat, i, j));
  //   printf("\n");
  // }

  // printf("\nlog prior:\n");
  // for(size_t i = 0; i < InitSamples; i++){
  //   gsl_vector_view v = gsl_matrix_row(Xmat, i);
  //   printf("%8.4f ", log(prior(&v.vector)));
  //   if((i+1) % 10 == 0) printf("\n");
  // }

  // printf("\nlog likelihood:\n");
  // for(size_t i = 0; i < InitSamples; i++){
  //   gsl_vector_view v = gsl_matrix_row(Xmat, i);
  //   printf("%8.4f ", log(likelihood(&v.vector)));
  //   if((i+1) % 10 == 0) printf("\n");
  // }


  #pragma omp parallel for
  for(size_t i = 0; i < InitSamples; i++){
    gsl_vector_const_view v = gsl_matrix_const_row(Xmat, i);
    likelihood(&v.vector);
    if(likelihood(&v.vector) > 0)
      printf("%7zu: %f %f\n", i, log(prior(&v.vector)), log(likelihood(&v.vector)));
  }

  gsl_matrix_free(Xmat);
  gsl_rng_free(rng);

  return 0;
}
