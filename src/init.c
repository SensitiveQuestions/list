#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void CombineEndorseListProbit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dpoisbinom(void *, void *, void *, void *);
extern void ictregBinom(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMixed(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMulti(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMulti2Level(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMulti3Level(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMulti4Level(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMultiMixed(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMultiProbit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ictregBinomMultiRobit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void R2Robit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RpoisbinomEff(void *, void *, void *, void *);
extern void RpoisbinomEffMatrix(void *, void *, void *, void *, void *, void *);
extern void RpoisbinomReturn(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"CombineEndorseListProbit", (DL_FUNC) &CombineEndorseListProbit, 55},
  {"dpoisbinom",               (DL_FUNC) &dpoisbinom,                4},
  {"ictregBinom",              (DL_FUNC) &ictregBinom,              22},
  {"ictregBinomMixed",         (DL_FUNC) &ictregBinomMixed,         35},
  {"ictregBinomMulti",         (DL_FUNC) &ictregBinomMulti,         23},
  {"ictregBinomMulti2Level",   (DL_FUNC) &ictregBinomMulti2Level,   29},
  {"ictregBinomMulti3Level",   (DL_FUNC) &ictregBinomMulti3Level,   39},
  {"ictregBinomMulti4Level",   (DL_FUNC) &ictregBinomMulti4Level,   49},
  {"ictregBinomMultiMixed",    (DL_FUNC) &ictregBinomMultiMixed,    36},
  {"ictregBinomMultiProbit",   (DL_FUNC) &ictregBinomMultiProbit,   21},
  {"ictregBinomMultiRobit",    (DL_FUNC) &ictregBinomMultiRobit,    22},
  {"R2Robit",                  (DL_FUNC) &R2Robit,                  10},
  {"RpoisbinomEff",            (DL_FUNC) &RpoisbinomEff,             4},
  {"RpoisbinomEffMatrix",      (DL_FUNC) &RpoisbinomEffMatrix,       6},
  {"RpoisbinomReturn",         (DL_FUNC) &RpoisbinomReturn,          4},
  {NULL, NULL, 0}
};

void R_init_list(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
