#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrice.h"


void annulle_matrice(double *tab, int lignes, int colonnes)
{
  int i, j;

  for (i=0; i < lignes; i++)
  {
    for (j=0; j < colonnes; j++)
    {

      *(tab + i * colonnes + j) = 0.0;

    } /* end for (j=0; j < colonnes; j++) */

  } /* end for (i=0; i < lignes; i++) */

} /* end void annulle_matrice(double *tab, int lignes, int colonnes) */

 /* Détermine la matrice inverse d'une matrice carrée non nulle. La méthode
    exploitée est la décomposition en deux matrices triangulaires (L - R)
	qui sont chacune facilement inversable et dont le produit donne la 
	matrice inverse cherchée.

    Cette méthode souffre du défaut qu'il faut que les mineurs de la diagonale
	ne soient pas null. Pour palier à cet inconvénient, on introduit au départ
	une légère modification de la matrice (ajout d'une petite valeur à tous les
	éléments de la diagonale, ce qui permet de toujours calculer. 

    Cette erreur volontaire est minimisée ensuite en corrigeant la matrice obtenue
	par itérations jusqu'a obtenir une erreur inférieure à un miliardième
	de l'erreur introduite au départ.

    Cette façon de faire n'est pas très scientifique, mais c'était écris pour le
	fun et de plus cela marche.
   
    Les entrées sont:
     - un pointeur sur le premier éléments de la matrice;
     - un entier donnant le nombre de ligne et de colonne de la matrice
       carrée;
     - un pointeur sur le permier élément de la matrice inverse. Cette 
       matrice doit être définie dans le code appelant la fonction.

    La fonction retourne le déterminant approché de la matrice à inverser.
	Si ce déterminant est nul, la matrice inverse déterminée est fausse (ou
	problème d'allocation dynamique de mémoire)
 
 */

double inverse_matrice(double *mat, int dim, double *mat_inv)
{
  
  double *U; /* matrice de travail */
  double *L; /* matrice de travail */
  double *Y; /* matrice de travail */
  double *Z; /* matrice de travail */
  double somme, erreur_acceptee, test;
  int i,j,k,p;

  if((U = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}
  if((L = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}
  if((Y = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}
  if((Z = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}

  /* annullation des matrices de travail */
  
  annulle_matrice((double*)U, dim, dim);  
  annulle_matrice((double*)L, dim, dim);
  annulle_matrice((double*)mat_inv, dim, dim);
  annulle_matrice((double*)Y, dim, dim);
  annulle_matrice((double*)Z, dim, dim);
  
  /* 
     détermination de l'ordre de grandeur de la matrice pour éviter
     les erreurs de division par zéro. Lorsque sur la diagonale, se 
     trouve des zéros ou si un même nombre forme un coin dans la matrice,
     il y a risque de division par zéro. Pour annuller ce risque, la  
     fonction ajoute aux éléments de la diagonale un nombre égal au
     milièmme de la valeur absolue non nulle la plus base. Cette 
     procédure introduit une erreur de calcul de l'ordre du milième 
     mais permet de calculer la matrice inverse de n'importe qu'elle 
     matrice carrée non nulle.
	 
	 le déterminant est une valeur approchée (au millième près)
	 
	 la matrice est ensuite corrigée par un algorhitme récursif
	 réalisé au maximum 10 fois ou jusqu'à ce que l'erreur entre
	 la matrice unitaire et la multiplication de la matrice par son
	 inverse soit plus petit qu'un miliardième de l'erreur acceptée.
  */
  
  erreur_acceptee = 1000.0;
  somme = 0.0;
  for(k=0;k<dim;k++)      /* détermination de la valeurs absolue non nulle */
  {                       /* la plus basse */
    for(i=0;i<dim;i++)
    {
    test = fabs(*(mat + k*dim +i));
    if(somme == 0.0)
      { 
        somme = test;
      } /* end if */
    if(test>0 && test< somme)
      {
        somme = test;
      } /* end if */
    } /* end for(i=0;i<dim;i++) */
  } /* end for(k=0;k<dim;k++) */
  erreur_acceptee = somme/erreur_acceptee;

  for(k=0;k<dim;k++)  /* modifie la diagonale de la matrice d'origine */
  {
    (*(mat+k*dim +k)) += erreur_acceptee; 
  }

  /* détermination des paramètres de L et U */

  for (k=0; k<dim; k++)
  {
    for (j=k; j<dim; j++)
    {
      somme = 0.0;
      
      for(p=0; p<k; p++)
      {
  		  somme += (*(L+k*dim+p)) * (*(U+p*dim+j));
      } /* end for(p=0; p<k; p++) */

  	  *(U+k*dim+j) = (*(mat + k*dim + j)) - somme;
   
      somme = 0.0;
      
      for(p=0; p<k; p++)
      {
  		  somme += (*(L+j*dim+p)) * (*(U+p*dim+k));
      } /* end for(p=0; p<k; p++) */

  	  *(L+j*dim+k) = 1.0;
      
      if(j!=k)
      {
  		  *(L+j*dim+k) = ((*(mat + j*dim + k)) - somme)/(*(U+k*dim+k));
      }

    } /* end for (j=k; j<dim; j++) */
  } /* end for (k=0; k<dim; k++) */

  for(k=0;k<dim;k++)  /* restitue la matrice d'origine */
  {
    (*(mat+k*dim+k)) -= erreur_acceptee;  
  }

  /* inverse les matrices L et U */

  for (j=0; j<dim; j++)
  {
  	  *(Y+j*dim+j)=1.0;

    for (i=j+1; i<dim; i++)
    {
      somme = 0.0;
      
      for(k=j; k<i; k++)
      {
  		  somme += (*(L+i*dim+k)) * (*(Y+k*dim+j));
      } /* end for(k=j; k<i; k++) */

  	  *(Y+i*dim+j) = (0.0 - somme)/(*(L+i*dim+i));

    } /* end for (i=j+1; i<dim; i++) */
  } /* end for (j=0; j<dim; j++) */

  for (j=0; j<dim; j++)
  {
  	  *(Z+j*dim+j)=1.0/(*(U+j*dim+j));

    for (i=j-1; i>-1; i--)
    {
      somme = 0.0;
      
      for(k=i+1; k<j+1; k++)
      {
  		somme += (*(U+i*dim+k)) * (*(Z+k*dim+j));
      } /* end for(k=j; k<i-1; k++) */

  	  *(Z+i*dim+j) = (0.0 - somme)/(*(U+i*dim+i));

    } /* end for (i=j+1; i<dim; i++) */
  } /* end for (j=0; j<dim; j++) */
  

/* calcul de la matrice inverse par multiplication de Z par Y */

  for (i=0; i<dim; i++)
  {  
    for (j=0; j<dim; j++)
      {
        for (k=0; k<dim; k++)  
        {
  			(*(mat_inv + i*dim +j)) += (*(Z+i*dim+k)) * (*(Y+k*dim+j));
        } /* end for (k=0; k<dim; k++) */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */

  /* 
     calcul du déterminant de la matrice par multiplication des
     éléments de la diagonale de U 
  */

  somme = 1.0;
  for(k=0;k<dim;k++)
  {
  	somme *= (*(U+k*dim+k));
  } /* end for(k=0;k<dim;k++) */
  
  /* l'introduction de l'erreur_acceptée ne permet plus de trouver un
     déterminant nul. Le test suivant impose un déterminant nul
     si celui-ci est de l'ordre de grandeur de l'erreur_acceptée */
  
  if(fabs(somme) < erreur_acceptee*10) 
    {
      somme = 0.0;
    } /* end if(somme < erreur_acceptee*10) */
	
/* correction de la matrice inverse si le déterminant n'est pas nul */

if (somme != 0.0)
  {

  annulle_matrice((double*)U, dim, dim);  
  annulle_matrice((double*)L, dim, dim);
  annulle_matrice((double*)Y, dim, dim);
  annulle_matrice((double*)Z, dim, dim);

  /* creation de la matrice unite U */
  for(k=0;k<dim;k++)
    {
  	  *(U+k*dim+k) = 1.0;
    } /* end for(k=0;k<dim;k++) */

    /* L est le resultat de la multiplication de la matrice par son inverse approchee */
    mutiplication_matrice((double*)mat, dim, dim, dim, (double*)mat_inv, (double*)L);
   /* Y est le resultat de la matrice unite (U) moins la matrice L */
    soustraction_matrice((double*)U, dim, dim, (double*)L, (double*)Y);
  
  k = 0;    
  	  test = 0.0;
      for(i=0;i<dim;i++)
       {
	     for(j=0;j<dim;j++)
		  {
			if(test<fabs(*(Y+i*dim+j))) {test = fabs(*(Y+i*dim+j));}
		  } /* end for(j=0;j<dim;j++) */
       } /* end for(i=0;i<dim;i++) */
	   
	   printf("pour k = %d, test = %e\n", k, test);

  do
    {
  
       /* Z est le resultat de la multiplication de la matrice inverse par Y */
      mutiplication_matrice((double*)mat_inv, dim, dim, dim, (double*)Y, (double*)Z);
  
      /* la matrice inverse corrigée est la somme de la matrice inverse et de Z */
      addition_matrice((double*)mat_inv, dim, dim, (double*)Z, (double*)mat_inv);
	  
	  /* test pour la fin du traitement sur base de la difference entre U et mat * mat_inv 
	     et preparation du calcul pour le tour suivant si necessaire */
    
	  /* L est le resultat de la multiplication de la matrice par son inverse approchee */
      mutiplication_matrice((double*)mat, dim, dim, dim, (double*)mat_inv, (double*)L);
  	  
     /* Y est le resultat de la matrice unite (U) moins la matrice L */
      soustraction_matrice((double*)U, dim, dim, (double*)L, (double*)Y);
	  
	  test = 0.0;
      for(i=0;i<dim;i++)
       {
	     for(j=0;j<dim;j++)
		  {
			if(test<fabs(*(Y+i*dim+j))) {test = fabs(*(Y+i*dim+j));}
		  } /* end for(j=0;j<dim;j++) */
       } /* end for(i=0;i<dim;i++) */
	   
	   k++;
 	   
    } while (test>erreur_acceptee * 1e-9 && k<10);
	
  } /* end if (somme != 0.0) */

  free(U);
  free(L);
  free(Y);
  free(Z);

  return(somme); /* si la somme est nulle, la matrice inverse est fausse */

} /* end double inverse_matrice(double *mat, int dim, double *mat_inv) */

/*  
  Multiplication de deux matrices mat1*mat2 de dimentions:
    
    - mat1[ligne1][colonne1]
    - mat2[colonne1][colonne2]

  dans une matrice resultat[ligne1][colonne2]

  Les données sont:

    - pointeur sur le premier élément de la matrice 1
    - entier du nombre de lignes de la matrice 1
    - entier du nombre de colonnes de la matrice 1 = nombre de lignes de la
      matrice 2
    - entier du nombre de colonnes de la matrice 2
    - pointeur sur le premier élément de la matrice 2
    - pointeur sur le premier élément de la matrice de resultat

*/

void mutiplication_matrice(double *mat1, int ligne1, int colonne1, int colonne2, double *mat2, double *resultat)
{
  int i, j, k;  

  annulle_matrice((double*)resultat, ligne1, colonne2); 

  for (i=0; i<ligne1; i++)
  {  
    for (j=0; j<colonne2; j++)
      {
        for (k=0; k<colonne1; k++)  
        {
          /* (*(resultat + i*ligne1 +j)) += (*(mat1 + i*ligne1 + k)) * (*(mat2 + k*colonne1 +j)); !!! faux faux pour matrice non carrée */
			(*(resultat + i*colonne2 +j)) += (*(mat1 + i*colonne1 + k)) * (*(mat2 + k*colonne2 +j)); /* RR2.0 : 21/06/2012 : correction bug */
        } /* end for (k=0; k<dim; k++) */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */


} /* end void multiplication_matrice(double *mat, int dim, double *mat_inv, double *unitaire) */
 
/*  
  addition de deux matrices mat1 et mat2 de dimentions:
    
    - mat1[lignes][colonnes]
    - mat2[lignes][colonnes]

  dans une matrice resultat[lignes][colonnes]

  Les données sont:

    - pointeur sur le premier élément de la matrice 1
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    - pointeur sur le premier élément de la matrice 2
    - pointeur sur le premier élément de la matrice de resultat

*/

void addition_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat)
{

int i, j;

  for (i=0; i<lignes; i++)
  {  
    for (j=0; j<colonnes; j++)
      {
          /* (*(resultat + i*lignes +j)) = (*(mat1 + i*lignes + j)) + (*(mat2 + i*lignes +j)); !!!!! faux pour matrice non carrée */
		  (*(resultat + i*colonnes +j)) = (*(mat1 + i*colonnes + j)) + (*(mat2 + i*colonnes +j)); /* RR2.0 : 21/06/2012 : correction bug */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */



} /* end void addition_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat) */
 
/*  
  soustraction de deux matrices mat1 - mat2 de dimentions:
    
    - mat1[lignes][colonnes]
    - mat2[lignes][colonnes]

  dans une matrice resultat[lignes][colonnes]

  Les données sont:

    - pointeur sur le premier élément de la matrice 1
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    - pointeur sur le premier élément de la matrice 2
    - pointeur sur le premier élément de la matrice de resultat

*/
 
void soustraction_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat)
{

int i, j;

  for (i=0; i<lignes; i++)
  {  
    for (j=0; j<colonnes; j++)
      {
         /* (*(resultat + i*lignes +j)) = (*(mat1 + i*lignes + j)) - (*(mat2 + i*lignes +j)); !!!!! faux pour matrice non carrée */
		  (*(resultat + i*colonnes +j)) = (*(mat1 + i*colonnes + j)) - (*(mat2 + i*colonnes +j)); /* RR2.0 : 21/06/2012 : correction bug */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */


} /* end void soustraction_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat)*/

/*  
  copie d'une matrice dans une autre deux matrices mat_copie - mat de dimentions identiques:
    
    - mat_copie[ligne][colonne]
    - mat2[ligne][colonne]

  
  Les données sont:

    - pointeur sur le premier élément de la matrice recevant la copie
	- pointeur sur le premier élément de la matrice à copier
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    
*/


void copie_matrice(double *mat_copie, double *mat, int ligne, int colonne)
{
	int i, j;

	for (i=0; i<ligne; i++)
	{  
		for (j=0; j<colonne; j++)
		{
			/* (*(mat_copie + i*ligne +j)) = (*(mat + i*ligne + j)); !!!!! faux pour matrice non carrée */
		  (*(mat_copie + i*colonne +j)) = (*(mat + i*colonne + j)); /* RR2.0 : 21/06/2012 : correction bug */
		} /* end for (j=0; j<dim; j++) */
	} /* end for (i=0; i<dim; i++) */
}

/*  
  multiplie une matrice par un réel:
          
  Les données sont:

    - pointeur sur le premier élément de la matrice à multiplier
	- un double = réel multiplicateur
    - entier du nombre de lignes de la matrice
    - entier du nombre de colonnes de la matrice
    
*/

void reelfoismatrice(double *mat, double reel, int ligne, int colonne) /* multiplie une matrice par un réel, le résultat est dans la matrice de départ */
{
	int i, j;

	for (i=0; i<ligne; i++)
	{  
		for (j=0; j<colonne; j++)
		{
			/* (*(mat + i*ligne +j)) *= reel; !!!!! faux pour matrice non carrée */
		   (*(mat + i*colonne +j)) *= reel; /* RR2.0 : 21/06/2012 : correction bug */
		} /* end for (j=0; j<dim; j++) */
	} /* end for (i=0; i<dim; i++) */

}

/*  
  addition de deux matrices mat et mat_add de dimentions identiques:
    
    - mat[lignes][colonnes]
    - mat_add[lignes][colonnes]

  le résultat étant placé dans la matrice mat[lignes][colonnes]

  Les données sont:

    - pointeur sur le premier élément de la matrice 1
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    - pointeur sur le premier élément de la matrice 2
    
*/


void addition_dans_matrice(double *mat, int lignes, int colonnes, double *mat_add) /* addition de mat+mat_add avec résultats dans mat */
{
  int i, j;

  for (i=0; i<lignes; i++)
  {  
    for (j=0; j<colonnes; j++)
      {
          /* (*(mat + i*lignes +j)) += (*(mat_add + i*lignes + j)); !!!!! faux pour matrice non carrée */
		   (*(mat + i*colonnes +j)) += (*(mat_add + i*colonnes + j)); /* RR2.0 : 21/06/2012 : correction bug */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */



}

/* inversion de la matrice à l'aide des coéfficients d'un polynome caractéristique:

  recherche du polynome (p) par la méthode de Leverrier, 

          1
   p  = - - * (s  + p  * s    + .... + p    * s )
    n     n     n    1    n-1           n-1    1

	
  puis, calcul de la matrice inverse par la formule:

   -1    1     n-1         n-2
  A  = - - * (A    + p  * A     + ..... + p   * E)
         p            1                    n-1
          n

  Les données sont:

     - un pointeur sur le premier éléments de la matrice;
     - un entier donnant le nombre de ligne et de colonne de la matrice
       carrée;
     - un pointeur sur le permier élément de la matrice inverse. Cette 
       matrice doit être définie dans le code appelant la fonction.

  retourne 1 si matrice est inversée ou 0 si défaut d'allocation dynamique de mémoire
  ou si un diviseur est nul

 */

int inverse_matrice_poly(double *mat, int dim, double *mat_inv)
{
	double *tr1, *tr2, *s, *p;
	int i, j;

	if((s = (double *) malloc(sizeof(double) * (dim+1)))==NULL) {return(0);} /* intermédiaires pour le calcul de p */
	if((p = (double *) malloc(sizeof(double) * (dim+1)))==NULL) {return(0);} /* coéfficients du polynome caractéristique de mat */
	if((tr1 = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0);} /* matrices de travail */
	if((tr2 = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0);} /* matrices de travail */

	annulle_matrice((double*)tr1, dim, dim);
	annulle_matrice((double*)tr2, dim, dim);

	/* calcul des coéfficients du polynome caractéristiques
	par le méthode de Leverrier(avec les matrices de puissance -> dim */

	/* les coéficients intermédiaires s sont les sommes des éléments des
	diagonales des matrices de puissance 1 à dim-1 */

	/* calcul de s1 */
	*(s) = 0.0;

	for (i=0; i < dim; i++)
	{
		 *(s) += *(mat + i * dim + i);
		  //printf("mat = %f\n", *(mat + i * dim + i));			  	    
	} /* end for (i=0; i < dim; i++) */

	copie_matrice(tr1, mat, dim, dim);

	/* calculs de s2 à s_dim */

	for(i=1; i<dim;i++)
	{
		mutiplication_matrice(tr1, dim, dim, dim, mat, tr2);
		*(s+i) = 0.0;

		for (j=0; j < dim; j++)
		{
			 *(s+i) += *(tr2 + j * dim + j);
		} /* end for (j=0; j < dim; j++) */

		copie_matrice(tr1, tr2, dim, dim);
	}

	/* détermination des coéficients p sur base des coéficients s */

	*(p) = -(*(s));

	for(i=1; i<dim;i++)
	{
		*(p+i) = *(s+i);
		for(j = 0; j<i; j++)
		{
			*(p+i)+= (*(p+j)) * (*(s+i-j-1));
		}

		(*(p+i))/= -(i+1);
	}

	if(*(p+dim-1) == 0)
	{
		/* inversion matrice impossible */
		free(s);
		free(p);
		free(tr1);
		free(tr2);
		return(0);
	}

	annulle_matrice((double*)mat_inv, dim, dim);

	for(i=0; i<dim;i++) /* matrice unité * p_n-1 */
	{
		*(mat_inv + i*dim + i) = *(p+dim-2);
	}

	copie_matrice(tr1, mat, dim, dim); /* copie mat dans tr1 */
	copie_matrice(tr2, mat, dim, dim); /* copie mat dans tr2 */

	for(i=2; i<dim; i++) /* calcul les éléments suivants */
	{
		reelfoismatrice(tr1, *(p+dim-(i+1)), dim, dim);
		addition_dans_matrice(mat_inv, dim, dim, tr1);
		copie_matrice(tr1, tr2, dim, dim); /* remet matrice puissance i dans tr1 avant mult par mat car tr1 modifier par fois reel */
		mutiplication_matrice(tr1, dim, dim, dim, mat, tr2);
		copie_matrice(tr1, tr2, dim, dim);		
	}

	addition_dans_matrice(mat_inv, dim, dim, tr1);
	reelfoismatrice(mat_inv, -1/(*(p+dim-1)), dim, dim);

	free(s);
	free(p);
	free(tr1);
	free(tr2);

	return(1); /* calcul terminé */

}

void print_matrice(double *tab, int lignes, int colonnes)
{
  int i, j;

  for (i=0; i < lignes; i++)
  {
    for (j=0; j < colonnes; j++)
    {

      printf("%f ", *(tab + i * colonnes + j));

    } /* end for (j=0; j < colonnes; j++) */

    printf("\n");

  } /* end for (i=0; i < lignes; i++) */

  printf("\n");

} /* end void print_matrice(double *tab, int lignes, int colonnes) */


/* inversion de la matrice par la méthode du pivot.

  le pivot est choisi dans chaque ligne comme étant le plus grand coéfficient non nul.

  La ligne de la matrice est dividée par le pivot, puis cette ligne est retranchée des
  autres lignes en la multipliant par le coéfficient de la ligne ayant le même indice
  que le pivot choisi. La colonne du pivot devient donc nulle sauf pour la ligne du pivot
  ou l'élément est égal à 1.

  On fait en parallèle la même chose dans le matrice unité de départ qui devient, en
  fin de calcul la matrice inverse dont les lignes sont mélangées. En effet, si pour la
  première ligne de la matrice, c'est l'élément x qui est le plus grand, le traitement
  donne le résultat pour la ligne x et non pas pour la première ligne. Comme toutes les
  lignes vont être traitée de 0 à n-1, les lignes de la matrice inverse seront mélangées.
  En gardant en mémoire l'ordre de traitement des éléments, il est facile en fin de
  calcul de la matrice inverse mélangées de recopier les lignes aux bon endroits dans
  la matrice inverse (dans notre cas, la première ligne doit être copiée dans la ligne
  x de la matrice inverse finale.

  Comme le système détruit la matrice d'entrée, celle-ci est copiée avant traitement.

  Le déterminant est le produit de tous les pivots.

  Si une ligne ne contient plus que des éléments nuls, la matrice n'est pas inversable
  et son déterminant est nul.
  
  Les données sont:

     - un pointeur sur le premier éléments de la matrice;
     - un entier donnant le nombre de ligne et de colonne de la matrice
       carrée;
     - un pointeur sur le permier élément de la matrice inverse. Cette 
       matrice doit être définie dans le code appelant la fonction.

 */

double inverse_matrice_pivot(double *mat, int dim, double *mat_inv) /* détermine la matrice inverse d'une matrice carrée non nulle par la méthode du pivot*/
{
	double *mtr, *melange; /* matrice de départ de travail et matrice inverse mélangée */
	int *memoire;
	int i, j, k, indice;
	double determinant = 1.0;
	double pivot, diviseur;
	double max;

	if((mtr = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);} /* matrice de travail */
	if((melange = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);} /* matrice inverse mélangée */
	if((memoire = (int *) malloc(sizeof(int) * dim))==NULL) {return(0.0);}; /* tableau de correspondance entre la ligne et la colonne traitée */
	copie_matrice(mtr, mat, dim, dim); /* copie mat dans mtr */
	annulle_matrice(melange, dim, dim); /* anulle tous les éléments de la matrice */
	annulle_matrice(mat_inv, dim, dim);
	
	for(i=0; i<dim;i++) /* matrice inverse mélangée commence comme matrice unité */
	{
		*(melange + i*dim + i) = 1.0;
	}
	
	/* inversion de la matrice par recherche dans chaque ligne du pivot le plus grand en absolu */

	for(i=0;i<dim;i++) /* on traite la ligne i */
	{
		indice = 0;
		
		/* recherche de l'indice de l'élément le plus grand en absolu 
		
		!!!! en prenant l'élément le plus grand, on ne traite plus la ligne correspondant
		 à la i, mais la ligne correspondant à l'indice de l'élément choisi. On construit donc la
		 matrice inverse avec des lignes mélangées. Il faut donc garder une trace de la 
		 correspondance entre la ligne traitée et la ligne correspondante dans la matrice
		 inverse finale pour par la suite recopier le tout dans le bon ordre.
		
		*/
		
		max = fabs(*(mtr+i*dim+indice));
		j = 1;
		
		while(j<dim)
		{
			if(fabs(*(mtr+i*dim+j)) > max)
			{
				max = fabs(*(mtr+i*dim+j));
				indice = j;
			}		
			j++;
		}
		
		if(max == 0.0) /* a trouver une ligne dont tous les coéfficients sont nuls */
		{
			free(mtr);
			free(melange);
			free(memoire);
			return(0.0); /* matrice n'est pas inversible car déterminant nul */
		}
		
		pivot = *(mtr+i*dim+indice);

		*(memoire+i) = indice; /* la ligne i correspond à la ligne finale indice */

		/* calcule le déterminant */
		determinant*=pivot;
		
		/* divise la ligne des deux matrices par le pivot */
		for(j = 0; j<dim ; j++)
		{
			*(mtr+i*dim+j)/=pivot;
			*(melange+i*dim+j)/=pivot;
		}
		
		/* soustrait la ligne du pivot des autres lignes * le coéfficient de la colonne du pivot */

		for(j=0; j<i; j++)
		{
			diviseur = *(mtr+j*dim+indice);
		
			for(k=0; k<dim; k++)
			{
				*(mtr+j*dim+k)-= diviseur * (*(mtr+i*dim+k));
				*(melange+j*dim+k)-= diviseur * (*(melange+i*dim+k));
			}			
		}

		for(j=i+1; j<dim; j++)
		{
			diviseur = *(mtr+j*dim+indice);
			
			for(k=0; k<dim; k++)
			{
				*(mtr+j*dim+k)-= diviseur * (*(mtr+i*dim+k));
				*(melange+j*dim+k)-= diviseur * (*(melange+i*dim+k));
			}			
		}		
	}

	/* copie des lignes de la matrice inverse mélangée à la bonne place dans la matrice
	inverse finale : ligne i de mélange = ligne *(memoire+i) de la matrice inverse  */

	for(i=0;i<dim;i++) /* on traite la ligne i */
	{
		for(j=0; j<dim; j++)
		{
			*(mat_inv+(*(memoire+i))*dim+j) = *(melange+i*dim+j);
		}

	}

	free(mtr);
	free(melange);
	free(memoire);

	return(determinant);
}
