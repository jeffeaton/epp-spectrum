#include <stdio.h>

#include "model.h"

int main()
{

  const size_t numSplines = 7;
  const double u[numSplines] = {1.5936094, 2.0196487, -0.4601538, -1.6417664, 1.4846658, 0.1725820, 0.3623600};

  double * rVec = fnGenRVec(u, numSplines);
  for(size_t i = 0; i < 426; i++){
    printf("%f ", rVec[i]);
    if((i+1) % 10 == 0)
      printf("\n");
  }

  fnFreeRVec(rVec);

  return 0;
 }
