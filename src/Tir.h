#include <stdlib.h>
#include <math.h>

int Suivant(int nbchiffres, int base, int *chiffres);

void TirUnicCop(const int *nbdimconnues, const int *nbdiminc, const int *tailsousech,
  double *unif,
  const double *cop, const int *rangconnues, const int *dimconnues, const int *dimincs,
  int *rangprevues);

void TirMultCop(const int *nbobsconnues, const int *nbdimconnues,
  const int *nbdiminc, const int *tailsousech,
  double *unif,
  const double *cop, const int *rangconnues, const int *dimconnues, const int *dimincs,
  int *inthreads,
  int *rangprevues);
