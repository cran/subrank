#include "Tir.h"

int Suivant(int nbchiffres, int base, int *chiffres)
{
  int position;
  chiffres[0]++;
  for( position=0; position<nbchiffres-1 ; position++ )
  {
    if (chiffres[position]==base)
    {
      chiffres[position]=0;
      chiffres[position+1]++;
    }
  }
  if (chiffres[nbchiffres-1]>=base) { return 0; }
  else return 1;
}

void TirUnicCop(const int *nbdimconnues, const int *nbdiminc, const int *tailsousech,
  double *unif,
  const double *cop, const int *rangconnues, const int *dimconnues, const int *dimincs,
  int *rangprevues)
{
  int position, NumEtat, NumConnues=0 ;
  double ProbaCondition=0;
  for( position=0; position<*nbdimconnues ; position++ )
  {
    NumConnues += rangconnues[position]*
      round(pow((double)*tailsousech,dimconnues[position]));
  }
  for( position=0; position<*nbdiminc ; position++ )
  { rangprevues[position]=0; }
  if (*nbdimconnues>1)
  {
    do
    {
      NumEtat=NumConnues;
      for( position=0; position<*nbdiminc ; position++ )
      {
        NumEtat+=
          rangprevues[position]*round(pow(*tailsousech,dimincs[position]));
      }
      ProbaCondition+=cop[NumEtat];
    }
    while(Suivant(*nbdiminc,*tailsousech,rangprevues));
  }
  else
  {
    if (*nbdimconnues==0) { ProbaCondition=1; }
    if (*nbdimconnues==1) { ProbaCondition=1/(double)*tailsousech; }
  }
  ProbaCondition=(*unif)*ProbaCondition;
  for( position=0; position<*nbdiminc ; position++ )
  { rangprevues[position]=0; }
  do
  {
    NumEtat=NumConnues;
    for( position=0; position<*nbdiminc ; position++ )
    {
      NumEtat+=
        rangprevues[position]*round(pow(*tailsousech,dimincs[position]));
    }
    ProbaCondition-=cop[NumEtat];
  }
  /* utilisation de la short-circuit evaluation, 6.5.13 ds norme C1X */
  while ( ProbaCondition>0 && Suivant(*nbdiminc,*tailsousech,rangprevues) );
}

void TirMultCop(const int *nbobsconnues, const int *nbdimconnues,
  const int *nbdiminc, const int *tailsousech,
  double *unif,
  const double *cop, const int *rangconnues, const int *dimconnues, const int *dimincs,
  int *inthreads,
  int *rangprevues)
{
  #pragma omp parallel num_threads(*inthreads)
  {
    int NumObs ;
    #pragma omp for schedule(dynamic)
    for( NumObs=0; NumObs<*nbobsconnues ; NumObs++ )
    {
      TirUnicCop(nbdimconnues,nbdiminc,tailsousech,
      &unif[NumObs],
      cop, &rangconnues[NumObs**nbdimconnues],dimconnues,dimincs,
      &rangprevues[NumObs**nbdiminc]);
    }
  }
}
