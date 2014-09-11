#include <R.h>
#include <Rinternals.h>

#include "states.h"
#include "parameters.h"
#include "model.h"
#include "incidence.h"

void setParameters(SEXP s_param, struct parameters * param);
SEXP getListElement(SEXP list, const char *str);


extern "C" {

  SEXP fnSpectrumR(SEXP s_param, SEXP s_numOutDates)
  {

    // set parameter values
    struct parameters param;

    setParameters(s_param, &param);

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

    fnSpectrum(&param, *INTEGER(s_numOutDates), Xout);


    // fnSpectrum(*iota, rVec, (unsigned int) *numOutDates, Xout);

    UNPROTECT(3);
    // free(param.rVec);

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

  SEXP s_iota, s_rVec, s_inc_sexrat, s_inc_agerat_male, s_inc_agerat_female, s_rmat;

  // iota
  s_iota = getListElement(s_param, "iota");
  if(s_iota == R_NilValue)
    error("iota not found");
  else
    param->iota = *REAL(s_iota);

  // rVec
  s_rVec = getListElement(s_param, "rVec");
  if(s_rVec == R_NilValue)
    error("rVec not found");
  else{
    param->rVec = (double *) malloc(length(s_rVec) * sizeof(double));
    for(size_t i = 0; i < length(s_rVec); i++){
      double * ptr_rVec = REAL(s_rVec);
      param->rVec[i] = ptr_rVec[i];
    }
  }

  // incidence model
  s_inc_agerat_male = getListElement(s_param, "inc.agerat.male");
  s_inc_agerat_female = getListElement(s_param, "inc.agerat.female");
  s_inc_sexrat = getListElement(s_param, "inc.sexrat");
  s_rmat = getListElement(s_param, "Rmat");

  if(s_inc_agerat_male != R_NilValue &&
     s_inc_agerat_female != R_NilValue &&
     s_inc_sexrat != R_NilValue &&
     s_rmat == R_NilValue){
    param->incmod = INC_INCRR;

    param->inc_sexrat = *REAL(s_inc_sexrat);

    double * ptr_inc_agerat_male = REAL(s_inc_agerat_male);
    double * ptr_inc_agerat_female = REAL(s_inc_agerat_female);
    for(size_t i = 0; i < AG; i++){
      param->inc_agerat[MALE][i] = ptr_inc_agerat_male[i];
      param->inc_agerat[FEMALE][i] = ptr_inc_agerat_female[i];
    }

  } else if(s_inc_agerat_male == R_NilValue &&
	    s_inc_agerat_female == R_NilValue &&
            s_inc_sexrat == R_NilValue &&
            s_rmat != R_NilValue){

    param->incmod = INC_RMAT;

    double * ptr_rmat = REAL(s_rmat);
    for(size_t i = 0; i < AG; i++)
      for(size_t j = 0; j < AG; j++)
        param->rmat[i][j] = ptr_rmat[i+AG*j];

  } else
    error("incidence model not found");

  return;
}
