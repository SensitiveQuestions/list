/* 
  list: Multivariate Statistical Analysis for the Item Count Technique
  Author: Graeme Blair and Kosuke Imai
  Efficient calculation of poisson binomial density
  using the method developed by Chen, Dempster, and Liu (1994) Biometrika
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"

/* direct use of recursive formula */
double Rpoisbinom(int k, double *p, int l) {

  double dtmp = 0.0, sumT;
  int i, j; 

  if (k == 0) {
    dtmp = 1.0;
  } else if (k > 0) {    
    dtmp = 0.0;
    for (i = 1; i <= k; i++) {
      sumT = 0.0;
      for (j = 0; j < l; j++) {
	sumT += R_pow_di(p[j]/(1-p[j]), i);
      }
      dtmp += R_pow_di(-1.0, i+1) * sumT * Rpoisbinom(k-i, p, l);
    } 
    dtmp /= k;
  } else {
    error("Rpoisbinom: invalid input for k.\n");
  }
  return(dtmp);
}

/* called from within R */
void RpoisbinomReturn(int *k, int *n, double *p, int *l, double *res) {

  int i; 

  for (i = 0; i < *n; i++) {
    res[i] = Rpoisbinom(k[i], p, *l);
  }
}

/* Poisson-binomial density function called from within R */
void dpoisbinom(int *y, int *n, double *p, int *k, double *res) {

  int i;
  double dtmp;

  dtmp = 1.0;
  for (i = 0; i < *k; i++) {
      dtmp *= (1-p[i]);
  }    
  for (i = 0; i < *n; i++) {
    res[i] = dtmp * Rpoisbinom(y[i], p, *k);
  }
}


/* alternative methods; tend to be faster when k is large */

void RpoisbinomEff(int *k, double *p, int *l, double *Rs) {

  double *sumT;
  int i, j; 
  sumT = doubleArray(*k);

  Rs[0] = 1.0;
  if (*k > 0) {
    for (i = 1; i <= *k; i++) {
      Rs[i] = 0.0;
      sumT[i-1] = 0.0;
      for (j = 0; j < *l; j++) {
	sumT[i-1] += R_pow_di(p[j]/(1-p[j]), i);
      }
      for (j = 1; j <= i; j++) {
	Rs[i] += R_pow_di(-1.0, j+1) * sumT[j-1] * Rs[i-j];
      }
      Rs[i] /= i;
    }
  }

  free(sumT);
}

/* this one genrealizes the previous code and accepts a vector of k (length m)
   and matrix of p (m copies of l dimensional vector) */
void RpoisbinomEffMatrix(int *k, int *maxk, double *p, int *l, int *m, double *Rs) {
  
  double ptmp, *dtmp, *sumT;
  int h, i, j;
  dtmp = doubleArray(*maxk+1); 
  sumT = doubleArray(*maxk);

  for (h = 0; h < *m; h++) {
    dtmp[0] = 1.0; 
    if (k[h] > 0) {
      for (i = 1; i <= k[h]; i++) {
	dtmp[i] = 0.0;
	sumT[i-1] = 0.0;
	for (j = 0; j < *l; j++) {
	  ptmp = p[h*l[0]+j];
	  sumT[i-1] += R_pow_di(ptmp/(1-ptmp), i);
	}
	for (j = 1; j <= i; j++) {
	  dtmp[i] += R_pow_di(-1.0, j+1) * sumT[j-1] * dtmp[i-j];
	}
	dtmp[i] /= i;
      }
    }
    Rs[h] = dtmp[k[h]];
  }
  
  free(dtmp);
  free(sumT);
}
