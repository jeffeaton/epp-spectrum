#include <R.h>
#include <Rinternals.h>

#include "model.h"


extern "C" {

  void fnSpectrumR(double * iota, double * rVec, int * numOutDates, double * Xout)
  {

    fnSpectrum(*iota, rVec, (unsigned int) *numOutDates, Xout);
    return;
  }

} // extern "C"
