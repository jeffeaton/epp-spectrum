#define AG_SPAN    5.0   // number of years spanned by age groups
#define PROJ_YEARS 50    // number of projection years, for time varying parametesr

#define IDX_FERT  3      // first fertility index
#define AG_FERT   7      // number of fertility age groups
#define IDX_15TO49 3
#define AG_15TO49  7
#define IDX_15PLUS 3
#define AG_15PLUS  14

#define IDX_500ELIG 2
#define IDX_350ELIG 3
#define IDX_200ELIG 5

#define MALE      0
#define FEMALE    1


// Incidence model flags
#define INC_INCRR  0
#define INC_RMAT   1


///////////////////////////
////  Declare methods  ////
///////////////////////////

// void alloc_parameters(struct parameters * param, const size_t proj_steps);
struct parameters * alloc_parameters(const size_t nps);
void free_parameters(struct parameters * param);


////////////////////////////
////  Define structure  ////
////////////////////////////
  
struct parameters {

  // fixed parameters over duration of simulation

  double dt;
  size_t n_proj_steps;
  size_t ts_epi_start;       // time step for HIV seeding
  double init_pop[NG][AG];
  
  
  // parameters that vary by ts
  double *proj_steps;
  double (*mx)[NG][AG];
  double (*asfr)[AG_FERT];
  double *srb;

  double (*artnum_15plus)[NG];
  size_t *artelig_idx;

  double (*agesex_incrr)[NG][AG];
  double *rvec;

  double age_fertrat[AG_FERT];
  double stage_fertrat[DS-1];
  double art_fertrat;

  // progression parameters
  double cd4_initdist[NG][AG][DS-1];
  double cd4_art_mort[NG][AG][DS-1][TS];
  double cd4_prog[NG][AG][DS-2];


  // incidence model parameters
  char incmod;  // = {INC_*}
  double iota;
  double rmat[AG][AG];
  double relinfect_art;      // relative infectousness of persons on ART
  double vert_trans;         // vertical transmission probability
  

  // declare varying parameters
  size_t ts;
  double iota_ts;

};


