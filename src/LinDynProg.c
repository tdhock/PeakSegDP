#include <R.h>
#include <Rinternals.h>
#include "math.h"


// Memory efficient and with constraint Poisson dynamic programing
// The constraint is \mu1 < \mu2 > \mu3 < \mu4 >....
// For the Poisson we don't need to consider the factorial part
// For the Poisson in the stand DP version we only to consider 
// the \sum x_i ln(lambda) [at the max likelihood \sum_i \hat{\lambda_i} = \sum x_i

void cPoisLinProgDyn
(int *sequence, int *weights, 
 int *lgSeq, int *nStep, 
 double *res1, int *res2, double *res3
  ){
  /* Lecture des données */
  /* On stocke les données comme un liste de vecteurs */
  //printf("ligne 11\n");
	
  char c = 13;
  //gsl_matrix_view matResult1 = gsl_matrix_view_array(res1, *nStep, *lgSeq);
  int i, j, k, l;
  i=0;
  while(i < *lgSeq-1){ 
    k=0;
    while(k < *nStep){
      res1[(*lgSeq)*k+i] = INFINITY;
      res3[(*lgSeq)*k+i] = INFINITY;
      k++;
    }
    i++;	
  }
	
  /* 	STRT INITIALISATION 	   */

  double SommeSeq, SommeWei, CostSeg0toi, MeanSeg0toi;

  // update Somme and SommeWei
  SommeSeq=weights[0]*sequence[0];
  SommeWei= weights[0];
	
  // compute cost and mean for next step
  if(SommeSeq != 0){
    // chech what it does when Somme = 0.
    CostSeg0toi = - SommeSeq * (log(SommeSeq) - log(SommeWei)); 
  } else {
    CostSeg0toi = 0;
  }
  MeanSeg0toi = SommeSeq/SommeWei;
	
  i=0;
  while(i < *lgSeq-1){
    //fprintf(stderr, "Boucle ligne 31 : %d \n", i);
    res1[i] = CostSeg0toi;
    res3[i] = MeanSeg0toi;
    res2[i] = 0;
		
    // update Somme and SommeWei.
    SommeSeq = SommeSeq + weights[i+1]*sequence[i+1];
    SommeWei= SommeWei + weights[i+1];
		
    // compute cost and mean for next step.
    if(SommeSeq != 0){
      // chech what it does when Somme = 0.
      CostSeg0toi = - SommeSeq * (log(SommeSeq) - log(SommeWei)); 
    } else {
      CostSeg0toi = 0;
    }
    MeanSeg0toi = SommeSeq/SommeWei;

    i++;
  }
  res1[i] = CostSeg0toi;
  res3[i] = MeanSeg0toi;
  res2[i] = 0;
	
  /* 	END INITIALISATION 	   */

  /* 	CONSIDER ALL SEG 	   */
	
  int minim;
  double coutTraj, CostSegitoj, MeanSegitoj;
  i=1;
  /* i : point de depart, j arrivée, k nombre de pas*/
  while(i < *lgSeq){
    //fprintf(stderr, "%c Node :   %d  / %d  ", c, i, *lgSeq);
    if( i+1 <= *nStep){ minim=i+1; } 
    else{ minim=*nStep; }
    //printf("Sommet %d et Mini %d\n", i, minim);

		
    SommeSeq=weights[i]*sequence[i];
    SommeWei= weights[i];

    if(SommeSeq != 0){
      // chech what it does when Somme = 0
      CostSegitoj = - SommeSeq * (log(SommeSeq) - log(SommeWei)); 
    } else {
      CostSegitoj = 0;
    }
    MeanSegitoj = SommeSeq/SommeWei;

    j=i+1;
    while(j < *lgSeq){
			
      k=1;
      R_CheckUserInterrupt();

      while(k < minim){
	if( ( (k%2 == 1) & 
	      ( MeanSegitoj  - res3[(*lgSeq)*(k-1)+i-1]  > 0)) | 
	    ((k%2 == 0) & 
	     ( MeanSegitoj  - res3[(*lgSeq)*(k-1)+i-1]  < 0)) ){ 
				
	  coutTraj = CostSegitoj + res1[(*lgSeq)*(k-1)+i-1];
	  //printf("I= %d, K=%d, J=%d, Poids: %f, %f, Means=%f - %f \n", i, k, j, CostSegitoj, coutTraj, MeanSegitoj, gsl_matrix_get(&matConstraint.matrix, k-1, i-1));
	  //printf("I= %d, K=%d, J=%d, Poids: %f, %f\n", i, k, j, Poids, coutTraj);
	  if(  coutTraj  < res1[(*lgSeq)*k+j-1]   ){
	    res1[(*lgSeq)*k+j-1] = coutTraj;
	    res3[(*lgSeq)*k+j-1] = MeanSegitoj;
	    res2[(*lgSeq)*k+j-1] = i;
	  }
	  //} else {
	  //printf("Not- I= %d, K=%d, J=%d, Poids: %f, %f, Means=%f - %f\n", i, k, j, CostSegitoj, coutTraj, MeanSegitoj, gsl_matrix_get(&matConstraint.matrix, k-1, i-1));
	}
	k++;
      }
		
      SommeSeq = SommeSeq + weights[j]*sequence[j];
      SommeWei= SommeWei + weights[j];
      if(SommeSeq != 0){
	// chech what it does when Somme = 0
	CostSegitoj = - SommeSeq * (log(SommeSeq) - log(SommeWei)); 
      } else {
	CostSegitoj = 0;
      }
      MeanSegitoj = SommeSeq/SommeWei;

				
      j++;
    }
		
    // last point //
    k=1;
    while(k < minim){
      if( ( (k%2 == 1) & ( MeanSegitoj  - res3[(*lgSeq)*(k-1)+i-1]  > 0)) | ((k%2 == 0) & ( MeanSegitoj  - res3[(*lgSeq)*(k-1)+i-1]  < 0)) ){ 
				
	coutTraj = CostSegitoj + res1[(*lgSeq)*(k-1)+i-1];
	//printf("I= %d, K=%d, J=%d, Poids: %f, %f, Means=%f - %f \n", i, k, j, CostSegitoj, coutTraj, MeanSegitoj, gsl_matrix_get(&matConstraint.matrix, k-1, i-1) );
	//printf("I= %d, K=%d, J=%d, Poids: %f, %f\n", i, k, j, Poids, coutTraj);
	if(  coutTraj  < res1[(*lgSeq)*k+j-1]  ){
	  res1[(*lgSeq)*k+j-1] = coutTraj;
	  res3[(*lgSeq)*k+j-1] = MeanSegitoj;
	  res2[(*lgSeq)*k+j-1] = i;
	}
	//} else {
	//	printf("Not- I= %d, K=%d, J=%d, Poids: %f, %f, Means=%f - %f \n", i, k, j, CostSegitoj, coutTraj, MeanSegitoj, gsl_matrix_get(&matConstraint.matrix, k-1, i-1));
      }
      k++;
    }
		
    i++;
  }
  //printf("\n End of Dynamic Programming\n");
  
}


