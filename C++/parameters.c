#include <stdlib.h>
#include "states.h"
#include "parameters.h"

// void alloc_parameters(struct parameters * p, const size_t nps){

struct parameters * alloc_parameters(const size_t nps){

  struct parameters * p = malloc(sizeof(struct parameters));
  if(!p)
    return NULL;

  p->proj_steps    = malloc(nps * sizeof(*p->proj_steps   ));
  p->mx            = malloc(nps * sizeof(*p->mx           ));
  p->asfr          = malloc(nps * sizeof(*p->asfr         ));
  p->srb           = malloc(nps * sizeof(*p->srb          ));
  p->artnum_15plus = malloc(nps * sizeof(*p->artnum_15plus));
  p->artelig_idx   = malloc(nps * sizeof(*p->artelig_idx  ));
  p->rvec          = malloc(nps * sizeof(*p->rvec         ));
  p->agesex_incrr  = malloc(nps * sizeof(*p->agesex_incrr ));

  return p;
}

void free_parameters(struct parameters * param){

  free(param->proj_steps);
  free(param->mx);
  free(param->asfr);
  free(param->srb);
  free(param->artnum_15plus);
  free(param->artelig_idx);
  free(param->rvec);
  free(param->agesex_incrr);

  free(param);
  return;
}
