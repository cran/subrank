#include "lienR.h"
#include "copc.h"
#include "Tir.h"

SEXP InterCopulation(
  SEXP Rrechant,
  SEXP Rimaxechant, SEXP Rimaxssech, SEXP Rimaxdim,
  SEXP Rimixties,
  SEXP Riu, SEXP Rimaxtir)
{
  Rrechant=coerceVector(Rrechant,REALSXP);
  double *rechant=REAL(Rrechant);
  int *imaxechant=INTEGER(Rimaxechant), *imaxssech=INTEGER(Rimaxssech),
    *imaxdim=INTEGER(Rimaxdim), *imaxtir=INTEGER(Rimaxtir), 
    *imixties=INTEGER(Rimixties), *iu=INTEGER(Riu);
  const int imaxcop=floor(.5+pow((int)*imaxssech,(int)*imaxdim));
  SEXP Ricop;
  int *icop;
  PROTECT(Ricop = allocVector(INTSXP, imaxcop+2));
  icop=INTEGER(Ricop);
  Copulation(rechant, imaxechant, imaxssech, imaxdim, imixties, iu, imaxtir, icop);
  UNPROTECT(1);
  return(Ricop);
}

SEXP InterPredFly(
  SEXP Rnbcomp, SEXP Rnbexps, SEXP Rnbinc, SEXP Rnbpreds,
  SEXP Rsubsampsize, SEXP Rmixties, SEXP Rmaxtirs,
  SEXP Rcompleteobs, SEXP Rincompleteobs)
{
  PROTECT(Rcompleteobs=coerceVector(Rcompleteobs,REALSXP));
  PROTECT(Rincompleteobs=coerceVector(Rincompleteobs,REALSXP));
  double *completeobs=REAL(Rcompleteobs), *incompleteobs=REAL(Rincompleteobs);
  int *nbcomp=INTEGER(Rnbcomp), *nbexps=INTEGER(Rnbexps),
    *nbinc=INTEGER(Rnbinc), *nbpreds=INTEGER(Rnbpreds),
    *subsampsize=INTEGER(Rsubsampsize),
    *mixties=INTEGER(Rmixties), *maxtirs=INTEGER(Rmaxtirs) ;
  SEXP Rcompletion;
  int *completion;
  PROTECT(Rcompletion = allocVector(INTSXP, *nbpreds * *nbinc));
  completion=INTEGER(Rcompletion);
  PredFly( nbcomp, nbexps, nbinc, nbpreds,
         subsampsize, mixties, maxtirs,
         completeobs, incompleteobs, completion );
  UNPROTECT(3);
  return(Rcompletion);
}

SEXP InterTir(
  SEXP Rnbobsconnues, SEXP Rnbdimconnues, SEXP Rnbdiminc, SEXP Rtailsousech, SEXP Runif,
  SEXP Rcop, SEXP Rrangconnues, SEXP Rdimconnues, SEXP Rdimincs
  )
{
  PROTECT(Runif=coerceVector(Runif,REALSXP));
  PROTECT(Rcop=coerceVector(Rcop,REALSXP));
  PROTECT(Rrangconnues=coerceVector(Rrangconnues,INTSXP));
  PROTECT(Rdimconnues=coerceVector(Rdimconnues,INTSXP));
  PROTECT(Rdimincs=coerceVector(Rdimincs,INTSXP));
  int *nbobsconnues=INTEGER(Rnbobsconnues), *nbdimconnues=INTEGER(Rnbdimconnues),
    *nbdiminc=INTEGER(Rnbdiminc),
    *tailsousech=INTEGER(Rtailsousech), *rangconnues=INTEGER(Rrangconnues),
    *dimconnues=INTEGER(Rdimconnues), *dimincs=INTEGER(Rdimincs) ;
  double *unif=REAL(Runif), *cop=REAL(Rcop);
  SEXP Rrangprevues;
  int *rangprevues;
  PROTECT(Rrangprevues = allocVector(INTSXP, *nbdiminc**nbobsconnues));
  rangprevues=INTEGER(Rrangprevues);
  TirMultCop(nbobsconnues,nbdimconnues,nbdiminc,tailsousech,unif,
  cop, rangconnues, dimconnues, dimincs,
  rangprevues);
  UNPROTECT(6);
  return(Rrangprevues);
}

#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
    CALLDEF(InterCopulation, 7),
    CALLDEF(InterPredFly, 9),
    CALLDEF(InterTir, 9),
    {NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_subrank(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
