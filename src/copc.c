#include "copc.h"


/*
conventions de nommage :
i pour entiers, r pour reels
parmi les entiers :
- indxxx pour un indice
- imax pour une valeur max d'indice
les noms sont en francais, sauf quand il s'agit d'ex-aequo,
pour lequel on prefere "tie".
*/


/*
www.iti.fh-flensburg.de/lang/algorithmen/sortieren/shell/shellen.htm
discutable de ne pas reutiliser un tri deja fait
mais une justification : itrace qui est une trace de la permutation
c'est la fonction order de R
*/
void Tri( double *rcle, int *itrace, int imax )
{
  const int isaut[]={4592,1968,861,336,112,48,21,7,3,1};
  const int imaxsaut=sizeof(isaut)/sizeof(isaut[0]);
  int k;
  for( k=0; k<imax ; k++ )
  { itrace[k]=k; }
  double vcle, vtrace ;
  for( k=0; k<imaxsaut; k++ )
  {
    int i;
    for( i=isaut[k]; i<imax; i++ )
    {
      vcle=rcle[i];
      vtrace=itrace[i];
      int j=i;
      while( j>=isaut[k] && rcle[j-isaut[k]]>vcle )
      {
        rcle[j]=rcle[j-isaut[k]];
        itrace[j]=itrace[j-isaut[k]];
        j=j-isaut[k];
      }
      rcle[j]=vcle;
      itrace[j]=vtrace;
    }
  }
}


/*
programmation de C_n^p, recursive et pas tres efficace,
mais pas important ici
*/
double NumComb( int n, int p )
{
  if( 2*p>n ) { return NumComb(n,n-p); }
  if( p>0 ) { return (NumComb(n,p-1)*(n-p+1))/p; }
  else { return 1; }
}


/*
OpenMP ne peut pas paralleliser une boucle while, seulement une boucle for.
Or, on connait le nombre de sous-ensembles de p elements parmi n.
Il faut donc enumerer dans l'ordre les sous-ensemble de p éléments parmi n.
On utilise une propriété connue: choisir p elements parmi n, c'est choisir
le nieme, puis en choisir p-1 parmi n-1$, ou bien exclure le nieme,
puis choisir p parmi n-1.
Il faut donc comparer m à C^{p-1}_{n-1}, et ensuite recurrer.
Programmation recursive pas tres efficace, mais pas important ici
car profondeur de la recursion toujours faible, sinon pas exhaustif
et deteministe
*/
void Combinaison( int *ssech, int numero, int n, int p )
{
  if ( p>0 )
  {
    int i, cnmoins1pmoins1 ;
    cnmoins1pmoins1=(int)NumComb(n-1,p-1);
    if ( p==n )
    { for( i=0;i<p;i++ ) { ssech[i]=i; }; }
    else if( numero<=cnmoins1pmoins1 )
      {
        ssech[p-1]=n-1;
        Combinaison(ssech,numero,n-1,p-1);
      }
      else{ Combinaison(ssech,numero-cnmoins1pmoins1,n-1,p); }
  }
}


/*
tirage d'un ss-echantillon
http://stackoverflow.com/questions/2394246/algorithm-to-select-a-single-random-combination-of-values
et avant doc SAS sur Surveyselect
*/
void TirSech( int *issech, rk_state *rkfil,
  const int imaxechant, const int imaxssech )
{
  int indechant, indajout=0 ;
  for( indechant=imaxechant-imaxssech; indechant<imaxechant; indechant++ )
  {
    int indverif, iajout, iegalite=0 ;
    iajout =rk_interval(indechant, rkfil);
    for( indverif=0 ; indverif<indajout ; indverif++ )
    { if (issech[indverif]==iajout) {iegalite=1;}; }
    if (iegalite==0) {issech[indajout]=iajout;}
    else  {issech[indajout]=indechant;}
    indajout++;
  }
}


/*
dans la suite, chaque point de la copule estimee est identifie par son
developpement en base imaxssech, ou la premiere coordonnee est celle
de plus fort poids.
Ajoutscopule prend un echantillon de dimension d, les indices
definissant un ss-echantillon, des tailles, et met dans un tableau
le nombre identifiant le point de la copule estimee ou il faut
ajouter 1.
Methode :
On voudrait mettre en i la valeur du rang r(i), mais on n'a pas calcule
directement le rang. La trace donnee par la fonction Tri est la permutation
reciproque du rang, qu'on peut noter o. Donc, au lieu de mettre
r(i) en i, on met i en o(i).
Pour le developpement en base imaxssech, multiplier tout le vecteur des
positions par imaxssech revient evidemment a decaler tous les chiffres vers la
gauche
*/
int Ajoutscopule( double *rechant, int *issech, int *iajoutscop,
  const int imaxssech, const int imaxechant, const int imaxdim )
{
  int indssech, inddim;
  double rechtemp[imaxdim][imaxssech];
  for( indssech=0; indssech<imaxssech; indssech++ )
  { 
    iajoutscop[indssech]=0;
    for( inddim=0; inddim<imaxdim; inddim++ )
    { 
      rechtemp[inddim][indssech]=
        rechant[inddim*imaxechant+issech[indssech]];
    }
  }
  for( inddim=0; inddim<imaxdim; inddim++ )
  {
    int permutation[imaxssech];
    Tri(&rechtemp[inddim][0],permutation,imaxssech);
    {
      for( indssech=1; indssech<imaxssech ; indssech++ )
      {
        if( rechtemp[inddim][indssech]==rechtemp[inddim][indssech-1] )
        { return 1; }
      }
    }
    for( indssech=0; indssech<imaxssech; indssech++ )
    {
      iajoutscop[permutation[indssech]]=
        iajoutscop[permutation[indssech]]*imaxssech+indssech;
    }
  }
  return 0;
}


/*
On doit sommer toutes les fois ou on atteint un point donne de la copule
discrete. Quelques remarques :
- Reduction sur array impossible en C, meme avec openMP 3.0.
- La copule discrete peut etre enorme.
- Remplissage de la copule discrete dans une region critical 
cause enorme perte de temps.
- Si utilisation de atomic, de moins en moins de collisions quand la copule
discrete grandit.
- Si echantillon et sous-echnatillon petits, calcul exhaustif possible.
Donc :
Si utilisateur acccepte de faire plus de calculs que ce que demanderait
l'exhaustivite, on fait exhaustif. Et alors, copies de copules.
Sinon, si copule assez petite (8000, par exemple), on copie, sinon, on ecrit
tous ensemble dans le meme tableau, en esperant qu'il est assez grand pour
qu'il y ait peu de collisions.

*/
void Copulation(
  double *rechant,
  int *imaxechant, int *imaxssech, int *imaxdim,
  int *iu, int *imaxtir,
  int *icop )
{
  /* initialisations */
  const int CopMaxRed = 8000 ;
  const int imaxcop=floor(.5+pow((int)*imaxssech,(int)*imaxdim)) ;
  {
    int indcop;
    for( indcop=0; indcop<imaxcop+1 ; indcop++ )
    {icop[indcop]=0;}
  }
  /* choix entre deterministe et stochastique */
  if ( NumComb(*imaxechant,*imaxssech)>*imaxtir )
  {
    if ( imaxcop>CopMaxRed )
    { 
      CopulationStoAto
        (rechant,*imaxechant,*imaxssech,*imaxdim,*iu,*imaxtir,icop);
    }
    else
    { 
      CopulationStoRed
        (rechant,*imaxechant,*imaxssech,*imaxdim,*iu,*imaxtir,icop);
    }
    icop[imaxcop+1]=*imaxtir;
  }
  else 
  {
    CopulationDet(rechant,*imaxechant,*imaxssech,*imaxdim,icop);
    icop[imaxcop+1]=NumComb(*imaxechant,*imaxssech);
  }
}


/*
calcul de copule par simulation,
avec utilisation de "atomic" pour faire une reduction
si nombre de sous-echantillons trop grand
*/
void CopulationStoAto( double *rechant,
  const int imaxechant, const int imaxssech, const int imaxdim,
  int iu, int imaxtir,
  int *icop )
{
  /* initialisations */
  int indtir, itie=0 ;
  #pragma omp parallel reduction (+:itie)
  {
    rk_state rkfil ;
    rk_randomseed(&rkfil);
    int iajoutscop[imaxssech], issech[imaxssech] ;
    #pragma omp for schedule(dynamic)
    for( indtir=0; indtir<imaxtir; indtir++ )
    {
      TirSech( issech, &rkfil, imaxechant, imaxssech );
      if( Ajoutscopule(rechant, issech, iajoutscop, imaxssech,
        imaxechant, imaxdim) == 0 )
      {
        int indpoint ;
        for( indpoint=0 ; indpoint<imaxssech ; indpoint++ )
        {
          #pragma omp atomic
          icop[iajoutscop[indpoint]]++;
        }
      }
      else { itie++; }
    }
  }
  icop[(int)floor(.5+pow(imaxssech,imaxdim))]=itie;
}


/*
calcul de copule par simulation,
avec une copie de la copule par thread
si nombre de sous-echantillons trop grand
*/
void CopulationStoRed( double *rechant,
  const int imaxechant, const int imaxssech, const int imaxdim,
  int iu, int imaxtir,
  int *icop )
{
  const int imaxcop=floor(.5+pow(imaxssech,imaxdim)) ;
  int indtir, itie=0 ; /* ici ou apres parallel ??? */
  /* enumeration des sous-echantillons */
  #pragma omp parallel reduction (+:itie)
  {
    rk_state rkfil ;
    rk_randomseed(&rkfil);
    int *icoppart ;
    icoppart=(int *)malloc (imaxcop * sizeof (int));
    {
      int indcop;
      for( indcop=0; indcop<imaxcop ; indcop++ )
        { icoppart[indcop]=0; }
    }
    #pragma omp for schedule(dynamic)
    for( indtir=0; indtir<imaxtir; indtir++ )
   {
      int indssech, issech[imaxssech], iajoutscop[imaxssech] ;
      TirSech( issech, &rkfil, imaxechant, imaxssech );
      if (Ajoutscopule(rechant, issech, iajoutscop,
        imaxssech, imaxechant, imaxdim)==0)
        for( indssech=0; indssech<imaxssech; indssech++ )
          {icoppart[ iajoutscop[indssech] ] ++ ;}
      else { itie++; }
    }
    #pragma omp critical
    {
      int indcop;
      for( indcop=0; indcop<imaxcop ; indcop++ )
        { icop[indcop] += icoppart[indcop]; }
    }
    free(icoppart);
  }
  icop[imaxcop]=itie;
}


/*
calcul de copule deterministe et exhaustif,
si nombre de sous-echantillons assez petit
avec une copie de la copule par thread
*/
void CopulationDet(double *rechant,
  const int imaxechant, const int imaxssech, const int imaxdim,
  int *icop)
{
  const int imaxcop=floor(.5+pow(imaxssech,imaxdim)) ;
  /* initialisations */
  int indcomb, itie=0 ;
  /* enumeration des sous-echantillons */
  #pragma omp parallel reduction (+:itie)
  {
    int *icoppart ;
    icoppart=(int *)malloc (imaxcop * sizeof (int));
    {
      int indcop;
      for( indcop=0; indcop<imaxcop ; indcop++ )
        { icoppart[indcop]=0; }
    }
    #pragma omp for schedule(static)
    for( indcomb=(int)NumComb(imaxechant,imaxssech);
      indcomb>0; indcomb-- )
    {
      int indssech, issech[imaxssech], iajoutscop[imaxssech] ;
      Combinaison(issech,indcomb,imaxechant,imaxssech);
      if (Ajoutscopule(rechant, issech, iajoutscop,
        imaxssech, imaxechant, imaxdim)==0)
        for( indssech=0; indssech<imaxssech; indssech++ )
          {icoppart[ iajoutscop[indssech] ] ++ ;}
      else { itie++; }
    }
    #pragma omp critical
    {
      int indcop;
      for( indcop=0; indcop<imaxcop ; indcop++ )
        { icop[indcop] += icoppart[indcop]; }
    }
    free(icoppart);
  }
  icop[imaxcop]=itie;
}
