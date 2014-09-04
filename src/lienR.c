#include "lienR.h"
#include "copc.h"
#include "Tir.h"

SEXP InterCopulation(
  SEXP Rrechant,
  SEXP Rimaxechant, SEXP Rimaxssech, SEXP Rimaxdim,
  SEXP Riu, SEXP Rimaxtir)
{
  Rrechant=coerceVector(Rrechant,REALSXP);
  double *rechant=REAL(Rrechant);
  int *imaxechant=INTEGER(Rimaxechant), *imaxssech=INTEGER(Rimaxssech),
    *imaxdim=INTEGER(Rimaxdim), *imaxtir=INTEGER(Rimaxtir), *iu=INTEGER(Riu);
  const int imaxcop=floor(.5+pow((int)*imaxssech,(int)*imaxdim));
  SEXP Ricop;
  int *icop;
  PROTECT(Ricop = allocVector(INTSXP, imaxcop+2));
  icop=INTEGER(Ricop);
  Copulation(rechant, imaxechant, imaxssech, imaxdim, iu, imaxtir, icop);
  UNPROTECT(1);
  return(Ricop);
}

SEXP InterTir(
  SEXP Rnbobsconnues, SEXP Rnbdimconnues, SEXP Rnbdiminc, SEXP Rtailsousech, SEXP Runif,
  SEXP Rcop, SEXP Rrangconnues, SEXP Rdimconnues, SEXP Rdimincs
  )
{
  Runif=coerceVector(Runif,REALSXP);
  Rcop=coerceVector(Rcop,REALSXP);
  Rrangconnues=coerceVector(Rrangconnues,INTSXP);
  Rdimconnues=coerceVector(Rdimconnues,INTSXP);
  Rdimincs=coerceVector(Rdimincs,INTSXP);
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
  UNPROTECT(1);
  return(Rrangprevues);
}
