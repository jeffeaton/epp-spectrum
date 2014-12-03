#include <R.h>
#include <Rinternals.h>

#include "states.h"
#include "parameters.h"
#include "model.h"
#include "incidence.h"

void setParameters(SEXP s_param, struct parameters * param);
void setFixedParameters(SEXP s_fixed_param, struct parameters * param);
SEXP getListElement(SEXP list, const char *str);


extern "C" {

  SEXP fnSpectrumR(SEXP s_param, SEXP s_fixpar, SEXP s_numOutDates)
  {

    // set parameter values
    struct parameters * param = alloc_parameters(length(getListElement(s_fixpar, "proj.steps")));
    setFixedParameters(s_fixpar, param);
    setParameters(s_param, param);

    // prepare output array
    PROTECT(s_numOutDates = coerceVector(s_numOutDates, INTSXP));
    SEXP s_Xout, s_Xout_dim;
    PROTECT(s_Xout = allocVector(REALSXP, *INTEGER(s_numOutDates) * NG * AG * DS * TS));
    PROTECT(s_Xout_dim = allocVector(INTSXP, 5));
    INTEGER(s_Xout_dim)[0] = *INTEGER(s_numOutDates);
    INTEGER(s_Xout_dim)[1] = NG;
    INTEGER(s_Xout_dim)[2] = AG;
    INTEGER(s_Xout_dim)[3] = DS;
    INTEGER(s_Xout_dim)[4] = TS;
    setAttrib(s_Xout, R_DimSymbol, s_Xout_dim);
    double * Xout = REAL(s_Xout);

    fnSpectrum(param, *INTEGER(s_numOutDates), Xout);


    UNPROTECT(3);
    free_parameters(param);

    return s_Xout;
  }

  SEXP createRmatR(SEXP s_m_mean, SEXP s_m_sd, SEXP s_f_mean, SEXP s_f_sd, SEXP s_corr)
  {
    double rmat[AG][AG];
    createRmat(*REAL(s_m_mean), *REAL(s_m_sd), *REAL(s_f_mean), *REAL(s_f_sd), *REAL(s_corr), rmat);
    
    // copy to output
    SEXP s_rmat;
    PROTECT(s_rmat = allocMatrix(REALSXP, AG, AG));
    double * ptr_rmat = REAL(s_rmat);
    for(size_t i = 0; i < AG; i++)
      for(size_t j = 0; j < AG; j++)
        ptr_rmat[i + AG*j] = rmat[i][j];
    
    UNPROTECT(1);
    return(s_rmat);
  }
  
} // extern "C"


SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for ( i = 0; i < length(list); i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  /*
    if ( elmt == R_NilValue )
    error("%s missing from list", str);
  */
  return elmt;
}

void setParameters(SEXP s_param, struct parameters * param)
{

  SEXP s_iota, s_rvec, s_inc_sexrat, s_inc_agerat, s_rmat;

  // iota
  s_iota = getListElement(s_param, "iota");
  if(s_iota == R_NilValue)
    error("iota not found");
  else
    param->iota = *REAL(s_iota);

   // rvec
  s_rvec = getListElement(s_param, "rvec");
  if(s_rvec == R_NilValue)
    error("rvec not found");
  else{
    double * ptr_rvec = REAL(s_rvec);
    for(size_t i = 0; i < param->n_proj_steps; i++)
      param->rvec[i] = ptr_rvec[i];
  }

  param->incmod = INC_INCRR;

  // incidence model
  /* s_inc_agerat = getListElement(s_param, "inc.agerat");
  s_inc_sexrat = getListElement(s_param, "inc.sexrat");
  s_rmat = getListElement(s_param, "Rmat");

  if(s_inc_agerat != R_NilValue &&
     s_inc_sexrat != R_NilValue &&
     s_rmat == R_NilValue){
    param->incmod = INC_INCRR;

    double * ptr_inc_sexrat = REAL(s_inc_sexrat);
    double * ptr_inc_agerat = REAL(s_inc_agerat);

    for(size_t i = 0; i < param->n_proj_steps; i++)
      for(size_t a = 0; a < AG; a++){
        size_t year_idx = (size_t) i*param->dt;
        param->agesex_incrr[i][MALE][a] = ptr_inc_agerat[MALE + a*NG + year_idx*NG*AG];
        param->agesex_incrr[i][FEMALE][a] = ptr_inc_sexrat[year_idx] * ptr_inc_agerat[FEMALE + a*NG + year_idx*NG*AG];
      }


  } else if(s_inc_agerat == R_NilValue &&
            s_inc_sexrat == R_NilValue &&
            s_rmat != R_NilValue){

    param->incmod = INC_RMAT;

    double * ptr_rmat = REAL(s_rmat);
    for(size_t i = 0; i < AG; i++)
      for(size_t j = 0; j < AG; j++)
        param->rmat[i][j] = ptr_rmat[i+AG*j];

  } else
    error("incidence model not found");
  */

  return;
}

void setFixedParameters(SEXP s_fixed_param, struct parameters * param)
{

  // TODO: Two potential efficiency improvemnets would be to access the
  //       list elements by index rather than name (no searching), and
  //       to reformulate parameters to assign them as pointers to
  //       contiguous blocks of memory rather than copying values into
  //       new arrays. I doubt either of these will have massive effects
  //       though, so not high priority.


  double *ptr_item;

  size_t nsteps = length(getListElement(s_fixed_param, "proj.steps"));
  param->n_proj_steps = nsteps;

  param->dt = *REAL(getListElement(s_fixed_param, "dt"));

  ptr_item = REAL(getListElement(s_fixed_param, "proj.steps"));
  for(size_t i = 0; i < nsteps; i++)
    param->proj_steps[i] = ptr_item[i];

  param->ts_epi_start = *INTEGER(getListElement(s_fixed_param, "ts.epi.start")) - 1; // -1 for 0-based indexing

  ptr_item = REAL(getListElement(s_fixed_param, "init.pop"));  // initial population size
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      param->init_pop[g][a] = ptr_item[g + a*NG];

  ptr_item = REAL(getListElement(s_fixed_param, "mx.ts"));
  for(size_t i = 0; i < nsteps; i++)
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
        param->mx[i][g][a] = ptr_item[g + a*NG + i*NG*AG];

  ptr_item = REAL(getListElement(s_fixed_param, "asfr.ts"));
  for(size_t i = 0; i < nsteps; i++)
    for(size_t a = 0; a < AG_FERT; a++)
      param->asfr[i][a] = ptr_item[a + i*AG_FERT];

  ptr_item = REAL(getListElement(s_fixed_param, "srb.ts"));
  for(size_t i = 0; i < nsteps; i++)
    param->srb[i] = ptr_item[i];

  ptr_item = REAL(getListElement(s_fixed_param, "agesex.incrr.ts"));
  for(size_t i = 0; i < nsteps; i++)
    for(size_t g = 0; g < NG; g++)
      for(size_t a = 0; a < AG; a++)
        param->agesex_incrr[i][g][a] = ptr_item[g + a*NG + i*NG*AG];

  ptr_item = REAL(getListElement(s_fixed_param, "age.fertrat"));
  for(size_t a = 0; a < AG_FERT; a++)
    param->age_fertrat[a] = ptr_item[a];

  ptr_item = REAL(getListElement(s_fixed_param, "stage.fertrat"));
  for(size_t m = 0; m < (DS-1); m++)
    param->stage_fertrat[m] = ptr_item[m];
  
  param->art_fertrat = *REAL(getListElement(s_fixed_param, "art.fertrat"));
    
  param->vert_trans = *REAL(getListElement(s_fixed_param, "vert.trans"));

  ptr_item = REAL(getListElement(s_fixed_param, "cd4.prog"));
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < (DS-2); m++)
        param->cd4_prog[g][a][m] = ptr_item[g + a*NG + m*NG*AG];

  ptr_item = REAL(getListElement(s_fixed_param, "cd4.initdist"));
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < (DS-1); m++)
        param->cd4_initdist[g][a][m] = ptr_item[g + a*NG + m*NG*AG];

  ptr_item = REAL(getListElement(s_fixed_param, "cd4.art.mort"));
  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < (DS-1); m++)
        for(size_t u = 0; u < TS; u++)
          param->cd4_art_mort[g][a][m][u] = ptr_item[g + a*NG + m*NG*AG + u*NG*AG*(DS-1)];

  param->relinfect_art = *REAL(getListElement(s_fixed_param, "relinfectART"));

  int *ptr_int = INTEGER(getListElement(s_fixed_param, "artelig.idx.ts"));
  for(size_t i = 0; i < nsteps; i++)
    param->artelig_idx[i] = ptr_int[i] - 1; // -1 for 0-based indexing

  ptr_item = REAL(getListElement(s_fixed_param, "artnum.15plus.ts"));
  for(size_t i = 0; i < nsteps; i++)
    for(size_t g = 0; g < NG; g++)
      param->artnum_15plus[i][g] = ptr_item[g + i*NG];

  return;
}
