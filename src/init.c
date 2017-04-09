/*Automatically generated from the routine:
 * package_native_routine_registration_skeleton("NPflow")
*/
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP NPflow_Fmeasure_costC(SEXP);
extern SEXP NPflow_FmeasureC(SEXP, SEXP);
extern SEXP NPflow_FmeasureC_no0(SEXP, SEXP);
extern SEXP NPflow_mmNiWpdfC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mmsNiWpdfC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mmvnpdfC(SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mmvsnpdfC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mmvstpdfC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mmvtpdfC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mvnlikC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mvnpdfC(SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mvsnlikC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_mvstlikC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NPflow_NuMatParC(SEXP, SEXP);
extern SEXP NPflow_sampleClassC(SEXP, SEXP);
extern SEXP NPflow_similarityMat_nocostC(SEXP);
extern SEXP NPflow_similarityMatC(SEXP);
extern SEXP NPflow_traceEpsC(SEXP, SEXP);
extern SEXP NPflow_vclust2mcoclustC(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"NPflow_Fmeasure_costC",        (DL_FUNC) &NPflow_Fmeasure_costC,        1},
  {"NPflow_FmeasureC",             (DL_FUNC) &NPflow_FmeasureC,             2},
  {"NPflow_FmeasureC_no0",         (DL_FUNC) &NPflow_FmeasureC_no0,         2},
  {"NPflow_mmNiWpdfC",             (DL_FUNC) &NPflow_mmNiWpdfC,             7},
  {"NPflow_mmsNiWpdfC",            (DL_FUNC) &NPflow_mmsNiWpdfC,            9},
  {"NPflow_mmvnpdfC",              (DL_FUNC) &NPflow_mmvnpdfC,              4},
  {"NPflow_mmvsnpdfC",             (DL_FUNC) &NPflow_mmvsnpdfC,             5},
  {"NPflow_mmvstpdfC",             (DL_FUNC) &NPflow_mmvstpdfC,             6},
  {"NPflow_mmvtpdfC",              (DL_FUNC) &NPflow_mmvtpdfC,              5},
  {"NPflow_mvnlikC",               (DL_FUNC) &NPflow_mvnlikC,               6},
  {"NPflow_mvnpdfC",               (DL_FUNC) &NPflow_mvnpdfC,               4},
  {"NPflow_mvsnlikC",              (DL_FUNC) &NPflow_mvsnlikC,              7},
  {"NPflow_mvstlikC",              (DL_FUNC) &NPflow_mvstlikC,              8},
  {"NPflow_NuMatParC",             (DL_FUNC) &NPflow_NuMatParC,             2},
  {"NPflow_sampleClassC",          (DL_FUNC) &NPflow_sampleClassC,          2},
  {"NPflow_similarityMat_nocostC", (DL_FUNC) &NPflow_similarityMat_nocostC, 1},
  {"NPflow_similarityMatC",        (DL_FUNC) &NPflow_similarityMatC,        1},
  {"NPflow_traceEpsC",             (DL_FUNC) &NPflow_traceEpsC,             2},
  {"NPflow_vclust2mcoclustC",      (DL_FUNC) &NPflow_vclust2mcoclustC,      1},
  {NULL, NULL, 0}
};

void R_init_NPflow(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}