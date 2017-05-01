#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP anMC_mvrnormArma(SEXP, SEXP, SEXP, SEXP);
extern SEXP anMC_trmvrnorm_rej_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"anMC_mvrnormArma",       (DL_FUNC) &anMC_mvrnormArma,       4},
  {"anMC_trmvrnorm_rej_cpp", (DL_FUNC) &anMC_trmvrnorm_rej_cpp, 6},
  {NULL, NULL, 0}
};

void R_init_anMC(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
