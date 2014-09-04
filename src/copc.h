#include <stdlib.h>
#include <math.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif
#include "randomkit.h"

void Tri( double *cle, int *trace, int lon );

int Ajoutscopule( double *rechant, int *issech, int *iajoutscop,
  const int imaxssech, const int imaxechant, const int imaxdim );

void TirSech( int *issech, rk_state *rkfil,
  const int imaxechant, const int imaxssech );

double NumComb( int n, int p );

void Combinaison( int *ssech, int numero, int n, int p );

void Copulation(double *rechant, 
  int *imaxechant, int *imaxssech, int *imaxdim,
  int *iu, int *imaxtir,
  int *icop);

void CopulationStoRed( double *rechant, 
  const int imaxechant, const int imaxssech, const int imaxdim,
  const int iu, const int imaxtir,
  int *icop );

void CopulationStoAto( double *rechant, 
  const int imaxechant, const int imaxssech, const int imaxdim,
  const int iu, const int imaxtir,
  int *icop );

void CopulationDet( double *rechant,
  const int imaxechant, const int imaxssech, const int imaxdim,
  int *icop );
