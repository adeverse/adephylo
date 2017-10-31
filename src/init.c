#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void distalltips(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gearymoran(void *, void *, void *, void *, void *, void *, void *);
extern void MVarianceDecompInOrthoBasis(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void VarianceDecompInOrthoBasis(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"distalltips",                 (DL_FUNC) &distalltips,                  8},
    {"gearymoran",                  (DL_FUNC) &gearymoran,                   7},
    {"MVarianceDecompInOrthoBasis", (DL_FUNC) &MVarianceDecompInOrthoBasis, 14},
    {"VarianceDecompInOrthoBasis",  (DL_FUNC) &VarianceDecompInOrthoBasis,  12},
    {NULL, NULL, 0}
};

void R_init_adephylo(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}