#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "imis.h"

int main(int argc, char *argv[])
{

  // argv: <initial samples> <step samples> <final resamples> <max iterations> <random seed> <simulation name>

  size_t B0 = (int) atof(argv[1]);
  size_t B  = (int) atof(argv[2]);
  size_t B_re  = (int) atof(argv[3]);
  size_t number_k = (int) atof(argv[4]);
  unsigned int rng_seed = (unsigned int) atof(argv[5]);
  size_t numParam = 9;

  fnIMIS(B0, B, B_re, number_k, numParam, rng_seed, argv[6]);

  return 0;
}
