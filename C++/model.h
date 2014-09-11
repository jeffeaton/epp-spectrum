/* void fnSpectrum(const double iota, const double * rVec, const double rmat[AG][AG], const size_t numOutDates, double * Xout); */
/* void fnSpectrumPrev(const double iota, const double * rVec, const double rmat[AG][AG], struct modprev * out); */

void fnSpectrum(struct parameters * param, const size_t numOutDates, double * Xout);
void fnSpectrumPrev(struct parameters * param, struct modprev * out);

double * fnGenRVec(const double * u, const size_t numSplines, const double maxR);
void fnFreeRVec(double * rVec);

/* structure for outputs to be compared in the likelihood */
#define NUM_OUTDATES 43

struct modprev{
  double ANCprev[NUM_OUTDATES];
  double a15to49prev[NUM_OUTDATES];
};
