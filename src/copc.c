#include "copc.h"
#include <stdio.h>


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

/* http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle */
void Permutation( int *permute, rk_state *rkfil, int n )
{
  int i, j, stockage ;
  for( i=n-1;i>0;i-- )
  {
    j = rk_interval(i, rkfil);
    stockage=permute[i];
    permute[i]=permute[j];
    permute[j]=stockage;
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
int Ajoutscopule( double *rechant, rk_state *rkfil, int *issech, int *iajoutscop,
  const int imaxssech, const int imaxechant, const int imaxdim,
  const int imixties )
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
      if ( imixties==0 )
      {
        for( indssech=1; indssech<imaxssech ; indssech++ )
        {
          if( rechtemp[inddim][indssech]==rechtemp[inddim][indssech-1] )
          { return 1; }
        }
      }
      else
      {
        int egaliteavant=0, egalite, debutmix, finmix ;
        for( indssech=1; indssech<imaxssech ; indssech++ )
        {
          egalite=(int)( rechtemp[inddim][indssech]==rechtemp[inddim][indssech-1] );
          if ( egalite==1 && egaliteavant==0 )
          { debutmix=indssech; }
          if ( (egalite==0 || indssech==imaxssech-1 ) && egaliteavant==1 )
          {
            finmix=indssech;
            Permutation( &permutation[debutmix-1], rkfil, finmix-debutmix+2);
          }
          egaliteavant=egalite;
        }
      }
      for( indssech=0; indssech<imaxssech; indssech++ )
      {
        iajoutscop[permutation[indssech]]=
          iajoutscop[permutation[indssech]]*imaxssech+indssech;
      }
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
  int *imixties,
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
        (rechant,*imaxechant,*imaxssech,*imaxdim,*imixties,*iu,*imaxtir,icop);
    }
    else
    { 
      CopulationStoRed
        (rechant,*imaxechant,*imaxssech,*imaxdim,*imixties,*iu,*imaxtir,icop);
    }
    icop[imaxcop+1]=*imaxtir;
  }
  else 
  {
    CopulationDet(rechant,*imaxechant,*imaxssech,*imaxdim,*imixties,icop);
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
  const int imixties,
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
      if( Ajoutscopule(rechant, &rkfil, issech, iajoutscop, imaxssech,
        imaxechant, imaxdim,imixties) == 0 )
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
  const int imixties,
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
      if (Ajoutscopule(rechant, &rkfil, issech, iajoutscop,
        imaxssech, imaxechant, imaxdim,imixties)==0)
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
  const int imixties,
  int *icop)
{
  const int imaxcop=floor(.5+pow(imaxssech,imaxdim)) ;
  /* initialisations */
  int indcomb, itie=0 ;
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
    #pragma omp for schedule(static)
    for( indcomb=(int)NumComb(imaxechant,imaxssech);
      indcomb>0; indcomb-- )
    {
      int indssech, issech[imaxssech], iajoutscop[imaxssech] ;
      Combinaison(issech,indcomb,imaxechant,imaxssech);
      if (Ajoutscopule(rechant, &rkfil, issech, iajoutscop,
        imaxssech, imaxechant, imaxdim,imixties)==0)
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

void PredFly( int *nbcomp, int *nbexps, int *nbinc, int *nbpreds,
  int *subsampsize, int *mixties, int *maxtirs,
  double *completeobs, double *incompleteobs, int *completion )
{
  const int bloctirs=1000 ;
  int permutinc[*nbexps**nbinc], tabnbpredsfaites[*nbinc],
    inddim, indinc, indtirs, resteapredire=1 ;
  for (indinc=0 ; indinc<*nbinc; indinc++ )
  { tabnbpredsfaites[indinc]=0; }
  for (indinc=*nbinc**nbpreds-1 ; indinc>=0; indinc-- )
  { completion[indinc]=-1; }
  for( inddim=0; inddim<*nbexps; inddim++ )
  { Tri(&incompleteobs[inddim**nbinc],&permutinc[inddim**nbinc],*nbinc); }
  for ( indtirs=0 ; indtirs<*maxtirs && resteapredire==1 ; indtirs+= bloctirs )
  {
    #pragma omp parallel
    {
      rk_state rkfil ;
      rk_randomseed(&rkfil) ;
      int *completemp ;
      completemp=(int *)malloc (bloctirs**nbinc*4 * sizeof (int));
      int nbpredsfaites, indtirs2 ;
      nbpredsfaites=0 ;
      #pragma omp for schedule(static)
      for ( indtirs2=0 ; indtirs2< bloctirs; indtirs2++ )
      {
        nbpredsfaites+=PredFlyUnic( *nbcomp, *nbexps, *nbinc,
          *subsampsize, *mixties,
          &rkfil, permutinc,
          completeobs, incompleteobs, &completemp[2*nbpredsfaites] ) ;
      }
      #pragma omp critical
      {
        int indpredfaite, indincourant ;
        for ( indpredfaite=0; indpredfaite<nbpredsfaites ; indpredfaite++ )
        {
          indincourant=completemp[2*indpredfaite+1];
          if (tabnbpredsfaites[indincourant]<*nbpreds)
          {
            completion[indincourant**nbpreds+tabnbpredsfaites[indincourant]]=
              completemp[2*indpredfaite] ;
            tabnbpredsfaites[indincourant]++; /* ??????????????????? */
          }
        }
        resteapredire=0;
        for (indincourant=0 ; indincourant<*nbinc; indincourant++ )
        { if (tabnbpredsfaites[indincourant]<*nbpreds) {
          resteapredire=1 ; } }
      }
      free(completemp);
    }
  }
}


int PredFlyUnic( int nbcomp, int nbexps, int nbinc,
  int subsampsize, int mixties,
  rk_state *rkfil, int *permutincloc,
  double *completeobs, double *incompleteobs, int *completion )
{
  int issech[subsampsize], permutss[nbexps][subsampsize],
    voisins[nbexps][nbinc][2],
    inddim, nbpredsfaites=0 ;
  double rechtemp[nbexps][subsampsize];
  TirSech( issech, rkfil, nbcomp, subsampsize );
  for( inddim=0; inddim<nbexps; inddim++ )
  {
    int indcomp ;
    for( indcomp=0; indcomp<subsampsize; indcomp++ )
    { 
      rechtemp[inddim][indcomp]=
        completeobs[inddim*nbcomp+issech[indcomp]];
    }
    Tri(&rechtemp[inddim][0],&permutss[inddim][0],subsampsize);
  }
  for( inddim=0; inddim<nbexps; inddim++ )
  {
    int ranginc=0, rangcomp=0,indobsqcq=0, indcompdernier,
      numobsqcq[nbexps][nbinc+subsampsize], obsinc[nbexps][nbinc+subsampsize]  ;
    while( rangcomp<subsampsize && ranginc<nbinc )
    {
      if ( rechtemp[inddim][rangcomp]<incompleteobs[inddim*nbinc+ranginc] )
      {
        numobsqcq[inddim][indobsqcq]=permutss[inddim][rangcomp];
        obsinc[inddim][indobsqcq]=0;
        rangcomp++;
        indobsqcq++;
      }
      else
      {
        numobsqcq[inddim][indobsqcq]=permutincloc[inddim*nbinc+ranginc];
        obsinc[inddim][indobsqcq]=1;
        ranginc++;
        indobsqcq++;
      }
    }
    while (rangcomp<subsampsize)
    {
      numobsqcq[inddim][indobsqcq]=permutss[inddim][rangcomp];
      obsinc[inddim][indobsqcq]=0;
      rangcomp++;
      indobsqcq++;
    }
    while (ranginc<nbinc)
    {
      numobsqcq[inddim][indobsqcq]=permutincloc[inddim*nbinc+ranginc];
      obsinc[inddim][indobsqcq]=1;
      ranginc++;
      indobsqcq++;
    }
    indcompdernier=-1;
    for(indobsqcq=0; indobsqcq<subsampsize+nbinc; indobsqcq++)
    { 
      if (obsinc[inddim][indobsqcq]==1)
      {
        voisins[inddim][numobsqcq[inddim][indobsqcq]][0]=indcompdernier;
      }
      else
      {
        indcompdernier=numobsqcq[inddim][indobsqcq];
      }
    }
    indcompdernier=-1;
    for(indobsqcq=subsampsize+nbinc-1; indobsqcq>=0; indobsqcq--)
    { 
      if (obsinc[inddim][indobsqcq]==1)
      {
        voisins[inddim][numobsqcq[inddim][indobsqcq]][1]=indcompdernier;
      }
      else
      {
        indcompdernier=numobsqcq[inddim][indobsqcq];
      }
    }
  }
  {
    int indinc, apres ;
    for ( indinc=0; indinc<nbinc; indinc++ )
    {
      for (apres=0 ; apres<=1; apres++)
      {
        for( inddim=1; inddim<nbexps; inddim++ )
        {
          if (voisins[0][indinc][apres]!=voisins[inddim][indinc][0] &&
              voisins[0][indinc][apres]!=voisins[inddim][indinc][1])
              break ;
        }
        if (inddim==nbexps && voisins[0][indinc][apres]>=0)
        {
          completion[2*nbpredsfaites]=issech[voisins[0][indinc][apres]];
          completion[2*nbpredsfaites+1]=indinc;
          nbpredsfaites++;
        }
      }      
    }
  }
  return nbpredsfaites ;
}
