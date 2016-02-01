#include <R.h>
#include <Rinternals.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "randomkit.h"

SEXP InterCopulation(
  SEXP Rrechant,
  SEXP Rimaxechant, SEXP Rimaxssech, SEXP Rimaxdim,
  SEXP Rimixties,
  SEXP Riu, SEXP Rimaxtir);

SEXP InterTir(
  SEXP Rnbobsconnues, SEXP Rnbdimconnues, SEXP Rnbdiminc, SEXP Rtailsousech,
  SEXP Runif,
  SEXP Rcop, SEXP Rrangconnues, SEXP Rdimconnues, SEXP Rdimincs);

SEXP InterPredFly(
  SEXP Rnbcomp, SEXP Rnbexps, SEXP Rnbinc, SEXP Rnbpreds,
  SEXP Rsubsampsize, SEXP Rmixties, SEXP Rmaxtirs, SEXP Rstep,
  SEXP Rcompleteobs, SEXP Rincompleteobs);

