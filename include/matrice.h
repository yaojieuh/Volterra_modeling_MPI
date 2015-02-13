/* matrice.h */  

/* 

bibliothèque d'utilisation des matrices

  il y a trois méthodes différentes pour l'inversion d'une matrice carrée.

  La fonction "inverse_matrice" est historiquement la première et avait été écrite
  juste pour voir si cela marchait. Testée valide jusqu'à une matrice de 600*600.
  Le temps d'inversion est le plus long des trois méthodes. A cause de la partie
  itérative pour arriver normalement à un meilleur résultat, il se pourrait que la
  matrice inverse soit fausse (jamais vu, mais pas impossible si itération ne converge
  pas).

  La fonction "inverse_matrice_poly" a été écrite car la méthode d'inversion est 
  jolie. A cause de la construction des matrices puissance, la méthode n'est exacte
  que jusqu'à une matrice de 50*50 (après la précision est trop faible et il y a
  risque de dépassement de capacité des double).

  La fonction "inverse_matrice_pivot" a été écrite pour être réellement utilisée.
  C'est la fonction la plus rapide en calcul et toute matrice inversable sera
  inversée quelque soient la disposition des éléments dans la matrice. Cette fonction
  a été testée valide jusqu'à une matrice de 1000*1000 (inversion en 190 secondes sur
  un PIII à 450 Hz)


*/

/*
  Historique du programme:/
  ************************

  R1.0: 03 septembre 1996 : première édition par BDC 
  R1.01 : 27 janvier 2000 : ajout addition, soustraction et
          modification de l'inversion (calcul corrige) par BDC
  R1.02 : 21 juin 2004 : inversion de matrice à l'aide des 
          coefficients d'un polynome caractéristique (BDC)
  R1.03 : 24 juin 2004 : inversion matrice par la méthode du pivot (BDC)
  R2.0  : 21 juin 2012 : correction bug pour matrices non carrées (merci à Cyprien DESCHEEMAEKER) (BDC)

*/

/* déclaration des structures */


/* déclaration des fonctions */

void annulle_matrice(double *tab, int lignes, int colonnes); /* annulation des éléments de la matrice */
double inverse_matrice(double *mat, int dim, double *mat_inv); /* détermine la matrice inverse d'une matrice carrée non nulle */
void mutiplication_matrice(double *mat1, int ligne1, int colonne1, int colonne2, double *mat2, double *resultat); /* multiplication des deux matrices */
void addition_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat);
void soustraction_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat);
int inverse_matrice_poly(double *mat, int dim, double *mat_inv); /* détermine la matrice inverse d'une matrice carrée non nulle méthode des polynomes caractéristique*/
void copie_matrice(double *mat_copie, double *mat, int ligne, int colonne); /* copie la matrice mat dans la matrice mat_copie */
void reelfoismatrice(double *mat, double reel, int ligne, int colonne); /* multiplie une matrice par un réel, le résultat est dans la matrice de départ */
void addition_dans_matrice(double *mat, int lignes, int colonnes, double *mat_add); /* addition de mat+mat_add avec résultats dans mat */
void print_matrice(double *tab, int lignes, int colonnes);
double inverse_matrice_pivot(double *mat, int dim, double *mat_inv); /* détermine la matrice inverse d'une matrice carrée non nulle par la méthode du pivot*/

	
