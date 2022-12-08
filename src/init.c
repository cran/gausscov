#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(decode)(int * j, int * k, int * set);
extern void F77_NAME(fstepwise)(double * y, double * x, int * n, int * k, double * x2, double * res, int * ia, double * alpha, int * kmx, double * pp, int * kex, int * intercept, double * minss, double * ss01, double * qq, int * kmn);
extern void F77_NAME(genint)(double * x, double * xx, int * n, int *k, int * kk, int * kex, int * ord, int * ind, int * ji);
extern void F77_NAME(graphst)(double * xxx, double * x, int * n, int * k, double * y, double * x2, double * res, int * ia, double * alpha, int * kmx, double * pp, double * grph, int * ne, int * kexc, double * xinr, double * minss, int * nedge, double * ss01, int * kmn, int * lin, int * iind, double * grphp);
extern void F77_NAME(graphstst)(double * xxx, double * x, int * n, int * k, double * y, double * x2, double * res, int * ia, double * alpha, int * kmx, double * pp, double * grph, int * ne, int * kexc, double * xinr, double * minss, int * nedge, double * ss01, double * rgrph, int * kmn, int * lin, int * iind);
extern void F77_NAME(lagg)(double * x, int * n, int * k, int * lag, double * xl, double * y, int * jj);
extern void F77_NAME(lmmdch)(double * y, double * x, int * n, int * k, double * xx, double * xxx, double * y1, double * y2, double * d, double * r, double * beta, double * xinv, int * ia, int * intercept, double * ss, int * nv, double * ssr, double * alpha, int * q);
extern void F77_NAME(triggen)(int * n, int * m, double tr);

static const R_FortranMethodDef FortranEntries[] = {
    {"decode",      (DL_FUNC) &F77_NAME(decode),       3},
    {"fstepwise",   (DL_FUNC) &F77_NAME(fstepwise),   16},
    {"genint",      (DL_FUNC) &F77_NAME(genint),       9},
    {"graphst",     (DL_FUNC) &F77_NAME(graphst),     22},
    {"graphstst",   (DL_FUNC) &F77_NAME(graphstst),   22},
    {"lagg",        (DL_FUNC) &F77_NAME(lagg),         7},
    {"lmmdch",      (DL_FUNC) &F77_NAME(lmmdch),      19},
    {"triggen",     (DL_FUNC) &F77_NAME(triggen),      3},
    {NULL, NULL, 0}
};

void R_init_gausscov(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
