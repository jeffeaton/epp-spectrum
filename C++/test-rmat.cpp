#include "stdio.h"

#include "states.h"
#include "incidence.h"

int main(){

  // double m_mean = 35.0, m_sd = 20.0, f_mean = 25.0, f_sd = 15.0, corr = 0.6;
  double m_mean = 14.936400024; 
  double m_sd = 0.001768254;
  double f_mean = 14.858318354;
  double f_sd = 32.424206654;
  double corr = 0.614077443;

  double rmat[AG][AG];

  createRmat(m_mean, m_sd, f_mean, f_sd, corr, rmat);


  for(size_t i = 0; i < AG; i++){
    for(size_t j = 0; j < AG; j++)
      printf("%7.4f ", rmat[i][j]);
    printf("\n");
  }


  return 0;

}
