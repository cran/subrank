#include <stdlib.h>
#include <math.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif
#include "randomkit.h"

void Tri( double *cle, int *trace, int lon );

int Ajoutscopule( double *rechant, rk_state *rkfil, int *issech, int *iajoutscop,
  const int imaxssech, const int imaxechant, const int imaxdim ,
  const int imixties ) ;

void Permutation( int *permute, rk_state *rkfil, int n );

void TirSech( int *issech, rk_state *rkfil,
  const int imaxechant, const int imaxssech );

double NumComb( int n, int p );

void Combinaison( int *ssech, int numero, int n, int p );

void PredFly( int *nbcomp, int *nbexps, int *nbinc, int *nbpreds,
  int *subsampsize, int *mixties, int *maxtirs,
  double *completeobs, double *incompleteobs, int *completion );

int PredFlyUnic( int nbcomp, int nbexps, int nbinc,
  int subsampsize, int mixties,
  rk_state *rkfil, int *permutinc,
  double *completeobs, double *incompleteobs, int *completion );

void Copulation(double *rechant, 
  int *imaxechant, int *imaxssech, int *imaxdim,
  int *imixties,
  int *iu, int *imaxtir,
  int *icop);

void CopulationStoRed( double *rechant, 
  const int imaxechant, const int imaxssech, const int imaxdim,
  const int imixties,
  const int iu, const int imaxtir,
  int *icop );

void CopulationStoAto( double *rechant, 
  const int imaxechant, const int imaxssech, const int imaxdim,
  const int imixties,
  const int iu, const int imaxtir,
  int *icop );

void CopulationDet( double *rechant,
  const int imaxechant, const int imaxssech, const int imaxdim,
  const int imixties,
  int *icop );
