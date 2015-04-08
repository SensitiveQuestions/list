#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

/** 
  Item Count Technique Binomial Regression for the Standard Design
**/

void ictregBinom(int *Y,             /* outcome vector */
		 int *J,             /* # of control items */
		 int *n_samp,        /* sample size */
		 int *n_draws,       /* # of MCMC draws */
		 int *treat,         /* treatment indicator vector: 0 or 1 */
		 double *Xall,       /* covariates in a vector form */
		 double *delta,      /* coefs for sensitive item */
		 double *psi,        /* coefs for control items */ 
		 int *n_cov,         /* # of covariates */
		 double *delta0,     /* prior mean for delta */
		 double *psi0,       /* prior mean for psi */
		 double *A0deltaAll, /* prior precision for delta */
		 double *A0psiAll,   /* prior precision for psi */
		 double *deltaVar,   /* proposal variance for delta */
		 double *psiVar,     /* proposal variance for psi */
		 int *unconst,       /* is this unconstrained model? */
		 int *ceiling,       /* ceiling effects */
		 int *floor,         /* floor effects */
		 int *burnin,        /* number of burnins */
		 int *keep,          /* keep every *th draw */
		 int *verbose,       /* want to print progress? */
		 double *allresults  /* storage for results */
		 ) {

  int i, j, n_covp, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int itempZstar0, itempZstar1;
  int *deltaCounter = intArray(1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);    /* acceptance ratio for psi */
  int *psi1Counter = intArray(1);   /* acceptance ratio for psi1 */
  double dtemp, dtemp1, dtemp2, dtemp3;

  /* intercept only model */
  if (*unconst == 2) {
    n_covp = *n_cov + 1;
  } else {
    n_covp = *n_cov;
  }
  if (((*ceiling == 1) || (*floor == 1)) && (*unconst > 0)) {
    error("Only constrained models are allowed with ceiling and floor effects\n");
  }

  /** get random seed **/
  GetRNGstate();
  /* additional parameters */
  double *psi1 = doubleArray(n_covp);
  double *psi01 = doubleArray(n_covp);
  
  if (*unconst == 1) {
    for (j = 0; j < n_covp; j++) {
      psi1[j] = psi[n_covp + j];
      psi01[j] = psi0[n_covp + j];;
    }
  }

  /** Data **/
  int *Zstar = intArray(*n_samp);
  int *Y0 = intArray(*n_samp);
  int *Y0temp0 = intArray(*n_samp);
  int *Y0temp1 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_covp); 
  double **Xtemp0 = doubleMatrix(*n_samp, n_covp); 
  double **Xtemp1 = doubleMatrix(*n_samp, n_covp); 
  double *Xdelta = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double *Xpsi1 = doubleArray(*n_samp);
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  /** Prior **/
  double **A0delta = doubleMatrix(*n_cov, *n_cov);
  double **A0psi = doubleMatrix(n_covp, n_covp);
  double **A0psi1 = doubleMatrix(n_covp, n_covp);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0delta[i][j] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < n_covp; j++)
    for (i = 0; i < n_covp; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  if (*unconst == 1) {
    itemp = 0;
    for (j = 0; j < n_covp; j++)
      for (i = 0; i < n_covp; i++)
	A0psi1[i][j] = A0psiAll[itemp++];
  }

  /** Proposal precision **/
  double **deltaPro = doubleMatrix(*n_cov, *n_cov);
  double **psiPro = doubleMatrix(n_covp, n_covp);
  double **psi1Pro = doubleMatrix(n_covp, n_covp);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      deltaPro[i][j] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < n_covp; j++)
    for (i = 0; i < n_covp; i++)
      psiPro[i][j] = psiVar[itemp++];

  if (*unconst == 1) {
    itemp = 0;
    for (j = 0; j < n_covp; j++)
      for (i = 0; i < n_covp; i++)
	psi1Pro[i][j] = psiVar[itemp++];
  }

  /** known Z star **/
  for (i = 0; i < *n_samp; i++) {
    if ((treat[i] == 1) && (Y[i] == (*J+1))) {
      Zstar[i] = 1;
      if (*ceiling == 1) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");
      }
    } else if ((treat[i] == 1) && (Y[i] == 0)) {
      Zstar[i] = 0;
      if (*floor == 1) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    } else { /* random draw if not known */
      Zstar[i] = (unif_rand() < 0.5);
    }
    if (Zstar[i] == 1) {
      Y0[i] = Y[i] - treat[i];
    } else {
      Y0[i] = Y[i];
    }
    if (*unconst == 2) 
      X[i][n_covp-1] = Zstar[i];
  }
  
  /** MCMC **/
  itemp = 0; deltaCounter[0] = 0; psiCounter[0] = 0; psi1Counter[0] = 0;
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... *** \n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    itempZstar0 = 0; itempZstar1 = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar */
      if ((treat[i] == 0) || ((treat[i] == 1) && (Y[i] < (*J + 1)) && (Y[i] > 0))) {
	Xdelta[i] = 0;  Xpsi[i] = 0;  Xpsi1[i] = 0;
	for (j = 0; j < *n_cov; j++) { 
	  Xdelta[i] += X[i][j]*delta[j];
	  Xpsi[i] += X[i][j]*psi[j];
	  if (*unconst == 1) {
	    Xpsi1[i] += X[i][j]*psi1[j];
	  } 
	}
	if (*unconst == 2) {
	  Xpsi1[i] = Xpsi[i] + psi[n_covp-1];
	}
	if ((*ceiling == 1) && (treat[i] == 1) && (Y[i] == *J)) { /* ceiling effects */
	  dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J-1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	  dtemp2 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  dtemp = unif_rand();
	  if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 1;
	    Y0[i] = *J-1;
	  } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 1;
	    Y0[i] = *J;
	  } else {
	    Zstar[i] = 0;
	    Y0[i] = *J;
	  }
	} else if ((*floor == 1) && (treat[i] == 1) && (Y[i] == 1)) { /* floor effects */
	  dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1));	  
	  dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	  dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	  dtemp = unif_rand();
	  if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 1;
	    Y0[i] = 0;
	  } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 0;
	    Y0[i] = 1;
	  } else {
	    Zstar[i] = 0;
	    Y0[i] = 0;
	  }
	} else { /* no ceiling and floor effects */
	  if (*unconst > 0) {
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i]-treat[i], *J, 1 / (1 + exp(-Xpsi1[i])), 1));
	  } else {
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i]-treat[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  }
	  dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	    Zstar[i] = 1;
	    Y0[i] = Y[i] - treat[i];
	  } else { 
	    Zstar[i] = 0;
	    Y0[i] = Y[i];
	  }
	}
	if (*unconst == 1) {
	  if (Zstar[i] == 1) {
	    Y0temp0[itempZstar0] = Y0[i];
	    for (j = 0; j < n_covp; j++) {
	      Xtemp0[itempZstar0][j] = X[i][j];
	    }
	    itempZstar0++;
	  } else {
	    Y0temp1[itempZstar1] = Y0[i];
	    for (j = 0; j < n_covp; j++) {
	      Xtemp1[itempZstar1][j] = X[i][j];
	    }
	    itempZstar1++;
	  }	    
	}
	if (*unconst == 2) { 
	  X[i][n_covp-1] = Zstar[i];
	}
      }
    }
    
    /* Sample delta */
    BinomLogit(Zstar, X, delta, *n_samp, 1, *n_cov, delta0, A0delta, deltaPro, 1, deltaCounter);
    
    /* Sample psi */
    if (*unconst == 1) {
      BinomLogit(Y0temp0, Xtemp0, psi, itempZstar0, *J, n_covp, psi0, A0psi, psiPro, 1, psiCounter);
      BinomLogit(Y0temp1, Xtemp1, psi1, itempZstar1, *J, n_covp, psi01, A0psi1, psi1Pro, 1, psi1Counter);
    } else {
      BinomLogit(Y0, X, psi, *n_samp, *J, n_covp, psi0, A0psi, psiPro, 1, psiCounter);
    }

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = delta[j];
	}      
	for (j = 0; j < n_covp; j++) {
	  allresults[itemp++] = psi[j];
	}
	if (*unconst == 1) {
	  for (j = 0; j < n_covp; j++) {
	    allresults[itemp++] = psi1[j];
	  }
	}
	allresults[itemp++] = ((double) *deltaCounter / (double) (main_loop + 1));
	allresults[itemp++] = ((double) *psiCounter / (double) (main_loop + 1));
	if (*unconst == 1) {
	  allresults[itemp++] = ((double) *psi1Counter / (double) (main_loop + 1));
	}
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("\n%3d percent done.\n    Metropolis acceptance ratios\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("      Sensitive item model: %3g\n", 
		fprec((double) *deltaCounter / (double) (main_loop + 1), 3));
	if (*unconst == 1) {
	  Rprintf("      Control items model: %3g (psi)  %3g (psi1)\n", 
		  fprec((double) *psiCounter / (double) (main_loop + 1), 3),
		  fprec((double) *psi1Counter / (double) (main_loop + 1), 3));
	} else {
	  Rprintf("      Control items model: %3g\n", 
		  fprec((double) *psiCounter / (double) (main_loop + 1), 3));
	}
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(Zstar);
  free(Y0);
  free(Y0temp0);
  free(Y0temp1);
  free(Xdelta);
  free(Xpsi);
  free(Xpsi1);
  free(psi1);
  free(psi01);
  free(deltaCounter);
  free(psiCounter);
  free(psi1Counter);
  FreeMatrix(X, *n_samp);
  FreeMatrix(Xtemp0, *n_samp);
  FreeMatrix(Xtemp1, *n_samp);
  FreeMatrix(A0delta, *n_cov);
  FreeMatrix(A0psi, n_covp);
  FreeMatrix(A0psi1, n_covp);
  FreeMatrix(deltaPro, *n_cov);
  FreeMatrix(psiPro, n_covp);
  FreeMatrix(psi1Pro, n_covp);
}


/** 
   Item Count Technique Binomial Regression for the Multiple Sensitive Item Design
**/


void ictregBinomMulti(int *Y,             /* outcome vector */
		      int *J,             /* # of control items */
		      int *n_samp,        /* sample size */
		      int *n_draws,       /* # of MCMC draws */
		      int *treat,         /* treatment indicator vector: 0, ..., tmax */
		      int *tmax,          /* number of sensitive items */
		      double *Xall,       /* covariates in a vector form */
		      double *delta,      /* coefs for sensitive item */
		      double *psi,        /* coefs for control items */ 
		      int *n_cov,         /* # of covariates */
		      double *delta0,     /* prior mean for delta */
		      double *psi0,       /* prior mean for psi */
		      double *A0deltaAll, /* prior precision for delta */
		      double *A0psiAll,   /* prior precision for psi */
		      double *deltaVar,   /* proposal variance for delta */
		      double *psiVar,     /* proposal variance for psi */
		      int *unconst,       /* is this unconstrained model? */
		      int *ceiling,       /* ceiling effects */
		      int *floor,         /* floor effects */
		      int *burnin,        /* number of burnins */
		      int *keep,          /* keep every *th draw */
		      int *verbose,       /* want to print progress? */
		      double *allresults  /* storage for results */
		      ) {

  int i, j, k, main_loop, itemp, itempS, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int **deltaCounter = intMatrix(*tmax, 1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);        /* acceptance ratio for psi */
  int *treatSum = intArray(*tmax);
  double dtemp, dtemp1, dtemp2, dtemp3;

  /** get random seed **/
  GetRNGstate();

  /** Parameters for sensitive items **/
  int n_dim = *n_cov;
  if (*unconst) { /* dimension of delta */
    n_dim = *n_cov + 1;
  }

  if (*unconst > 0) {
    itemp = 0;
    for (i = 0; i < *tmax; i++) {
      itemp += ceiling[i];
      itemp += floor[i];
    }
    Rprintf("%5d\n", itemp);
    if (itemp > 0)
      error("Only constrained models are allowed with ceiling and floor effects\n");
  }

  double **deltaMatrix = doubleMatrix(*tmax, n_dim);

  itemp = 0; 
  for (i = 0; i < *tmax; i++) 
    for (j = 0; j < n_dim; j++)
      deltaMatrix[i][j] = delta[itemp++];

  /** Data **/
  int **Zstar = intMatrix(*tmax, *n_samp);
  int *Y0 = intArray(*n_samp);
  double *Xdelta = doubleArray(*n_samp);
  double *Xdelta1 = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_dim); 
  double ***Xtemp = doubleMatrix3D(*tmax, *n_samp, n_dim); 
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  for (i = 0; i < *tmax; i++)
    treatSum[i] = 0;

  for (i = 0; i < *n_samp; i++) {
    Y0[i] = Y[i];
    if (treat[i] > 0) {
      for (j = 0; j < *n_cov; j++) 
	Xtemp[treat[i]-1][treatSum[treat[i]-1]][j] = X[i][j];
      treatSum[treat[i]-1]++;
      if ((ceiling[treat[i]-1] == 1) && (Y[i] == (*J+1))) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");	
      }
      if ((floor[treat[i]-1] == 1) && (Y[i] == 0)) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    }
  }

  /** Prior **/
  double ***A0delta = doubleMatrix3D(*tmax, n_dim, n_dim);
  double **A0psi = doubleMatrix(*n_cov, *n_cov);
  double **m0delta = doubleMatrix(*tmax, n_dim);

  itemp = 0;
  for (i = 0; i < *tmax; i++)
    for (k = 0; k < n_dim; k++)
      for (j = 0; j < n_dim; j++)
	A0delta[i][j][k] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  itemp = 0;
  for (i = 0; i < *tmax; i++)
    for (j = 0; j < n_dim; j++)
      m0delta[i][j] = delta0[itemp++];

  /** Proposal precisoin **/
  double ***deltaPro = doubleMatrix3D(*tmax, n_dim, n_dim);
  double **psiPro = doubleMatrix(*n_cov, *n_cov);

  itemp = 0;
  for (i = 0; i < *tmax; i++)
    for (k = 0; k < n_dim; k++)
      for (j = 0; j < n_dim; j++)
	deltaPro[i][j][k] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      psiPro[i][j] = psiVar[itemp++];

  /** MCMC **/
  itempS = 0; psiCounter[0] = 0; 
  for (i = 0; i < *tmax; i++) {
    deltaCounter[i][0] = 0; 
  }
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... *** \n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *tmax; j++)
      treatSum[j] = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar for treated units */
      if (treat[i] > 0) {
	if (Y[i] == (*J+1)) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0;  Xpsi[i] = 0;  
	  for (j = 0; j < *n_cov; j++) {
	    Xdelta[i] += X[i][j] * deltaMatrix[treat[i]-1][j];
	    Xpsi[i] += X[i][j] * psi[j];
	  }
	  if ((ceiling[treat[i]-1] == 1) && (Y[i] == *J)) { /* ceiling effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J-1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp2 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J-1;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = *J;
	    }
	  } else if ((floor[treat[i]-1] == 1) && (Y[i] == 1)) { /* floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = 0;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 1;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 0;
	    }
	  } else { /* no ceiling and floor effects */
	    if (*unconst) {
	      Xdelta1[i] = Xdelta[i] + (Y[i] - 1) * deltaMatrix[treat[i]-1][*n_cov];
	      Xdelta[i] += Y[i] * deltaMatrix[treat[i]-1][*n_cov];
	      dtemp1 = exp(Xdelta1[i] - log1p(exp(Xdelta1[i])) + 
			   dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    } else {
	      dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			   dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    } 
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = Y[i] - 1;
	    } else { 
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = Y[i];
	    }
	  }
	}
	if (*unconst) {
	  X[i][*n_cov] = Y0[i];
	  Xtemp[treat[i]-1][treatSum[treat[i]-1]][*n_cov] = Y0[i];
	}
	treatSum[treat[i]-1]++;
      }
    }

    /* Sample delta */
    for (k = 0; k < *tmax; k++) {
      BinomLogit(Zstar[k], Xtemp[k], deltaMatrix[k], treatSum[k], 1, n_dim, 
		 m0delta[k], A0delta[k], deltaPro[k], 1, deltaCounter[k]);
    }

    /* Sample psi */
    BinomLogit(Y0, X, psi, *n_samp, *J, *n_cov, psi0, A0psi, psiPro, 1, psiCounter);

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	for (k = 0; k < *tmax; k++) {
	  for (j = 0; j < n_dim; j++) {
	    allresults[itempS++] = deltaMatrix[k][j];
	  }      
	}
	for (j = 0; j < *n_cov; j++) {
	  allresults[itempS++] = psi[j];
	}
	for (k = 0; k < *tmax; k++) {
	  allresults[itempS++] = ((double) deltaCounter[k][0] / (double) (main_loop + 1));
	}
	allresults[itempS++] = ((double) *psiCounter / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n    Metropolis acceptance ratios\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("      Sensitive item model: %3g", 
		fprec((double) deltaCounter[0][0] / (double) (main_loop + 1), 3));
	for (k = 1; k < *tmax; k++) {
	  Rprintf(" %3g", fprec((double) deltaCounter[k][0] / (double) (main_loop + 1), 3));
	}
	Rprintf("\n    Control items model: %3g\n", 
		fprec((double) *psiCounter / (double) (main_loop + 1), 3));
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(psiCounter);
  free(treatSum);
  free(Y0);
  free(Xdelta);
  free(Xdelta1);
  free(Xpsi);
  FreeMatrix(m0delta, *tmax);
  FreeintMatrix(deltaCounter, *tmax);
  FreeMatrix(deltaMatrix, *tmax);
  FreeintMatrix(Zstar, *tmax);
  FreeMatrix(X, *n_samp);
  Free3DMatrix(Xtemp, *tmax, *n_samp);
  Free3DMatrix(A0delta, *tmax, n_dim);
  FreeMatrix(A0psi, *n_cov);
  Free3DMatrix(deltaPro, *tmax, n_dim);
  FreeMatrix(psiPro, *n_cov);
}


/** 
  Item Count Technique Binomial Mixed Effects Regression for the Standard Design
  
  unconst = 0 or 2 is allowed

**/

void ictregBinomMixed(int *Y,             /* outcome vector */
		      int *J,             /* # of control items */
		      int *n_samp,        /* sample size */
		      int *n_draws,       /* # of MCMC draws */
		      int *treat,         /* treatment indicator vector: 0 or 1 */
		      double *Xall,       /* fixed effects covariates in a vector form */
		      double *delta,      /* fixed effects coefs for sensitive item */
		      double *psi,        /* fixed effects coefs for control items */ 
		      int *n_cov,         /* # of fixed effects covariates */
		      double *delta0,     /* prior mean for delta */
		      double *psi0,       /* prior mean for psi */
		      double *A0deltaAll, /* prior precision for delta */
		      double *A0psiAll,   /* prior precision for psi */
		      double *deltaVar,   /* proposal variance for delta (fixed effects) */
		      double *psiVar,     /* proposal variance for psi (fixed effects) */
		      int *unconst,       /* is this unconstrained model? */
		      int *ceiling,       /* ceiling effects? */
		      int *floor,         /* floor effects? */
		      int *grp,           /* group indicator, 0, 1, ..., G-1 */
		      int *n_grp,         /* number of groups, G */
		      int *max_samp_grp,  /* max # of obs within each group */
		      double *Xall_rand,  /* random effects covariates */
		      int *n_rand,        /* # of random effects covariates */
		      double *tune_gamma, /* tuning constants for random effects (sensitive) */
		      double *tune_zeta,  /* tuning constants for random effects (control) */
		      double *SigmaAll,   /* covariance for random effects (sensitive) */
		      double *PhiAll,     /* covariance for random effects (control) */
		      int *s0,            /* prior df for Sigma (sensitive) */
		      double *S0All,      /* prior scale matrix for Sigma (sensitive) */
		      int *t0,            /* prior df for Phi (control) */
		      double *T0All,      /* prior scale matrix for Phi (control) */
		      int *burnin,        /* number of burnins */
		      int *keep,          /* keep every *th draw */
		      int *verbose,       /* want to print progress? */
		      double *allresults  /* storage for all results */
		      ) {

  int i, j, k, n_covp, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int *deltaCounter = intArray(1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);    /* acceptance ratio for psi */
  int *gammaCounter = intArray(*n_grp);
  int *zetaCounter = intArray(*n_grp);
  int *vitemp = intArray(*n_grp);
  double dtemp, dtemp1, dtemp2, dtemp3;

  /* intercept only model */
  if (*unconst == 2) {
    n_covp = *n_cov + 1;
  } else {
    n_covp = *n_cov;
  }
  if (((*ceiling == 1) || (*floor == 1)) && (*unconst > 0)) {
    error("Only constrained models are allowed with ceiling and floor effects\n");
  }

  /** get random seed **/
  GetRNGstate();

  /** parameters for random effects **/
  double **gamma = doubleMatrix(*n_grp, *n_rand); /* sensitive */
  double *gamma0 = doubleArray(*n_rand);
  double **Sigma = doubleMatrix(*n_rand, *n_rand);
  double **SigInv = doubleMatrix(*n_rand, *n_rand);
  double **zeta = doubleMatrix(*n_grp, *n_rand);  /* control */
  double *zeta0 = doubleArray(*n_rand);
  double **Phi = doubleMatrix(*n_rand, *n_rand);
  double **PhiInv = doubleMatrix(*n_rand, *n_rand);

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      Sigma[i][j] = SigmaAll[itemp++];
  dinv(Sigma, *n_rand, SigInv);
  
  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      Phi[i][j] = PhiAll[itemp++];
  dinv(Phi, *n_rand, PhiInv);

  /* starting values for random effects from prior */
  for (j = 0; j < *n_rand; j++) {
    gamma0[j] = 0; 
    zeta0[j] = 0;
  }
  for (j = 0; j < *n_grp; j++) { 
    rMVN(gamma[j], gamma0, Sigma, *n_rand);
    rMVN(zeta[j], zeta0, Phi, *n_rand);
  }

  /** Data **/
  int *Zstar = intArray(*n_samp);
  int *Y0 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_covp); 
  double *Xdelta = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  double *Xpsi1 = doubleArray(*n_samp);
  double ***Vgrp = doubleMatrix3D(*n_grp, *max_samp_grp, *n_rand);
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++)
      X[i][j] = Xall[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++) {
    vitemp[j] = 0;
  }
  for (i = 0; i < *n_samp; i++) {
    for (j = 0; j < *n_rand; j++)
      Vgrp[grp[i]][vitemp[grp[i]]][j] = Xall_rand[itemp++];
    vitemp[grp[i]]++;
  }
 
  /* PdoubleMatrix3D(Vgrp, *n_grp, *max_samp_grp, *n_rand); */

  /** Prior **/
  double **A0delta = doubleMatrix(*n_cov, *n_cov);
  double **A0psi = doubleMatrix(n_covp, n_covp);
  double **S0 = doubleMatrix(*n_rand, *n_rand);
  double **T0 = doubleMatrix(*n_rand, *n_rand);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0delta[i][j] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < n_covp; j++)
    for (i = 0; i < n_covp; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      S0[i][j] = S0All[itemp++];

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      T0[i][j] = T0All[itemp++];

  /** Proposal precisoin **/
  double **deltaPro = doubleMatrix(*n_cov, *n_cov);
  double **psiPro = doubleMatrix(n_covp, n_covp);

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      deltaPro[i][j] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < n_covp; j++)
    for (i = 0; i < n_covp; i++)
      psiPro[i][j] = psiVar[itemp++];

  /** known Z star **/
  for (i = 0; i < *n_samp; i++) {
    if ((treat[i] == 1) && (Y[i] == (*J+1))) {
      Zstar[i] = 1;
      Y0[i] = *J;
      if (*ceiling == 1) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");
      }
    } else if ((treat[i] == 1) && (Y[i] == 0)) {
      Zstar[i] = 0;
      Y0[i] = 0;
      if (*floor == 1) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    } else { /* random draw if not known */
      Zstar[i] = (unif_rand() < 0.5);
    }
    if (Zstar[i] == 1)
      Y0[i] = Y[i] - treat[i];
    else 
      Y0[i] = Y[i];
    if (*unconst == 2) 
      X[i][n_covp-1] = Zstar[i];
  }
  
  /** MCMC **/
  itemp = 0; deltaCounter[0] = 0; psiCounter[0] = 0; 
  for (j = 0; j < *n_grp; j++) {
    gammaCounter[j] = 0;
    zetaCounter[j] = 0;
  }
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... ***\n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar */
      if ((treat[i] == 0) || ((treat[i] == 1) && (Y[i] < (*J + 1)) && (Y[i] > 0))) {
	Xdelta[i] = 0;  Xpsi[i] = 0;  
	for (j = 0; j < *n_cov; j++) { /* fixed effects */
	  Xdelta[i] += X[i][j]*delta[j];
	  Xpsi[i] += X[i][j]*psi[j];
	}
	for (j = 0; j < *n_rand; j++) { /* random effects */
	  Xdelta[i] += Vgrp[grp[i]][vitemp[grp[i]]][j] * gamma[grp[i]][j];
	  Xpsi[i] += Vgrp[grp[i]][vitemp[grp[i]]][j] * zeta[grp[i]][j];
	}
	if (*unconst == 2) 
	  Xpsi1[i] = Xpsi[i] + psi[n_covp-1]; 
	if ((*ceiling == 1) && (treat[i] == 1) && (Y[i] == *J)) { /* ceiling effects */
	  dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J-1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	  dtemp2 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  dtemp = unif_rand();
	  if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 1;
	    Y0[i] = *J-1;
	  } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 1;
	    Y0[i] = *J;
	  } else {
	    Zstar[i] = 0;
	    Y0[i] = *J;
	  }
	} else if ((*floor == 1) && (treat[i] == 1) && (Y[i] == 1)) { /* floor effects */
	  dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1));	  
	  dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	  dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	  dtemp = unif_rand();
	  if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 1;
	    Y0[i] = 0;
	  } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	    Zstar[i] = 0;
	    Y0[i] = 1;
	  } else {
	    Zstar[i] = 0;
	    Y0[i] = 0;
	  }
	} else { /* no ceiling and floor effects */
	  if (*unconst == 2) {
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i]-treat[i], *J, 1 / (1 + exp(-Xpsi1[i])), 1));
	  } else {
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i]-treat[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  }
	  dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	  if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	    Zstar[i] = 1;
	    Y0[i] = Y[i] - treat[i];
	  } else { 
	    Zstar[i] = 0;
	    Y0[i] = Y[i];
	  }
	}
	if (*unconst == 2) {
	  X[i][n_covp-1] = Zstar[i];
	}
      }
      vitemp[grp[i]]++;
    }

    /* Sample delta */
    BinomLogitMixed(Zstar, X, Vgrp, grp, delta, gamma, SigInv, *n_samp, 1, *n_cov, *n_rand, *n_grp, 
		    delta0, A0delta, *s0, S0, deltaPro, tune_gamma, 1, deltaCounter, gammaCounter);
    
    /* Sample psi */
    BinomLogitMixed(Y0, X, Vgrp, grp, psi, zeta, PhiInv, *n_samp, *J, n_covp, *n_rand, *n_grp, 
		    psi0, A0psi, *t0, T0, psiPro, tune_zeta, 1, psiCounter, zetaCounter);

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	/* fixed effects */
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = delta[j];
	}      
	for (j = 0; j < n_covp; j++) {
	  allresults[itemp++] = psi[j];
	}
	/* Var-Cov matrices */
	dinv(SigInv, *n_rand, Sigma);
	for (j = 0; j < *n_rand; j++) {
	  for (k = j; k < *n_rand; k++) {
	    allresults[itemp++] = Sigma[j][k];
	  }      
	}
	dinv(PhiInv, *n_rand, Phi);
	for (j = 0; j < *n_rand; j++) {
	  for (k = j; k < *n_rand; k++) {
	    allresults[itemp++] = Phi[j][k];
	  }      
	}
	/* random effects */
	for (j = 0; j < *n_grp; j++) {
	  for (k = 0; k < *n_rand; k++) {
	    allresults[itemp++] = gamma[j][k];
	  }
	}
	for (j = 0; j < *n_grp; j++) {
	  for (k = 0; k < *n_rand; k++) {
	    allresults[itemp++] = zeta[j][k];
	  }
	}
	/* acceptance ratios */
	allresults[itemp++] = ((double) *deltaCounter / (double) (main_loop + 1));
	allresults[itemp++] = ((double) *psiCounter / (double) (main_loop + 1));
	for (j = 0; j < *n_grp; j++)
	  allresults[itemp++] = ((double) gammaCounter[j] / (double) (main_loop + 1));
	for (j = 0; j < *n_grp; j++)
	  allresults[itemp++] = ((double) zetaCounter[j] / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("\n%3d percent done.\n    Metropolis acceptance ratios", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("\n    Fixed effects: %3g (sensitive) %3g (control)", 
		fprec((double) *deltaCounter / (double) (main_loop + 1), 3),
		fprec((double) *psiCounter / (double) (main_loop + 1), 3));
	Rprintf("\n    Sensitive item model random effects (for each group):");
	for (j = 0; j < *n_grp; j++) 
	  Rprintf(" %3g", fprec((double) gammaCounter[j] / (double) (main_loop + 1), 3));
	Rprintf("\n    Control items model random effects (for each group):");
	for (j = 0; j < *n_grp; j++) 
	  Rprintf(" %3g", fprec((double) zetaCounter[j] / (double) (main_loop + 1), 3));
	Rprintf("\n");
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(Zstar);
  free(Y0);
  free(Xdelta);
  free(Xpsi);
  free(Xpsi1);
  free(deltaCounter);
  free(psiCounter);
  free(gammaCounter);
  free(zetaCounter);
  free(vitemp);
  free(gamma0);
  free(zeta0);
  FreeMatrix(Sigma, *n_rand);
  FreeMatrix(SigInv, *n_rand);
  FreeMatrix(Phi, *n_rand);
  FreeMatrix(PhiInv, *n_rand);
  FreeMatrix(X, *n_samp);
  FreeMatrix(A0delta, *n_cov);
  FreeMatrix(A0psi, *n_cov);
  FreeMatrix(S0, *n_rand);
  FreeMatrix(T0, *n_rand);
  FreeMatrix(deltaPro, *n_cov);
  FreeMatrix(psiPro, *n_cov);
  Free3DMatrix(Vgrp, *n_grp, *max_samp_grp);
}

/** 
  Item Count Technique Binomial Mixed Effects Regression for the Multiple Sensitive Item Design
  
  for now, constrained model only

**/

void ictregBinomMultiMixed(int *Y,             /* outcome vector */
			   int *J,             /* # of control items */
			   int *n_samp,        /* sample size */
			   int *n_draws,       /* # of MCMC draws */
			   int *treat,         /* treatment indicator vector: 0 or 1 */
			   int *tmax,          /* number of sensitive items */
			   double *Xall,       /* fixed effects covariates in a vector form */
			   double *delta,      /* fixed effects coefs for sensitive item */
			   double *psi,        /* fixed effects coefs for control items */ 
			   int *n_cov,         /* # of fixed effects covariates */
			   double *delta0,     /* prior mean for delta */
			   double *psi0,       /* prior mean for psi */
			   double *A0deltaAll, /* prior precision for delta */
			   double *A0psiAll,   /* prior precision for psi */
			   double *deltaVar,   /* proposal variance for delta (fixed effects) */
			   double *psiVar,     /* proposal variance for psi (fixed effects) */
			   int *unconst,       /* unconst = 1 includes Y(0) */
			   int *ceiling,       /* ceiling effects */
			   int *floor,         /* floor effects */
			   int *grp,           /* group indicator, 0, 1, ..., G-1 */
			   int *n_grp,         /* number of groups, G */
			   int *max_samp_grp,  /* max # of obs within each group */
			   double *Xall_rand,  /* random effects covariates */
			   int *n_rand,        /* # of random effects covariates */
			   double *tune_gamma, /* tuning constants for random effects (sensitive) */
			   double *tune_zeta,  /* tuning constants for random effects (control) */
			   double *SigmaAll,   /* covariance for random effects (sensitive) */
			   double *PhiAll,     /* covariance for random effects (control) */
			   int *s0,            /* prior df for Sigma (sensitive) */
			   double *S0All,      /* prior scale matrix for Sigma (sensitive) */
			   int *t0,            /* prior df for Phi (control) */
			   double *T0All,      /* prior scale matrix for Phi (control) */
			   int *burnin,        /* number of burnins */
			   int *keep,          /* keep every *th draw */
			   int *verbose,       /* want to print progress? */
			   double *allresults  /* storage for all results */
			   ) {
  
  int i, j, k, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int **deltaCounter = intMatrix(*tmax, 1);  /* acceptance ratio for delta */
  int *psiCounter = intArray(1);    /* acceptance ratio for psi */
  int **gammaCounter = intMatrix(*tmax, *n_grp);
  int *zetaCounter = intArray(*n_grp);
  int *vitemp = intArray(*n_grp);
  int *treatSum = intArray(*tmax);
  double dtemp, dtemp1, dtemp2, dtemp3;

  /** get random seed **/
  GetRNGstate();

  /** Parameters for sensitive items **/
  int n_dim = *n_cov;
  if (*unconst) { /* dimension of delta */
    n_dim = *n_cov + 1;
  }

  if (*unconst > 0) {
    itemp = 0;
    for (i = 0; i < *tmax; i++) {
      itemp += ceiling[i];
      itemp += floor[i];
    }
    if (itemp > 0)
      error("Only constrained models are allowed with ceiling and floor effects\n");
  }

  /** parameters for fixed effects **/
  double **deltaMatrix = doubleMatrix(*tmax, n_dim);
  double **delta0Matrix = doubleMatrix(*tmax, n_dim);

  itemp = 0; 
  for (i = 0; i < *tmax; i++) 
    for (j = 0; j < n_dim; j++)
      deltaMatrix[i][j] = delta[itemp++];

  itemp = 0; 
  for (i = 0; i < *tmax; i++) 
    for (j = 0; j < n_dim; j++)
      delta0Matrix[i][j] = delta0[itemp++];

  /** parameters for random effects **/
  double ***gamma = doubleMatrix3D(*tmax, *n_grp, *n_rand); /* sensitive */
  double *gamma0 = doubleArray(*n_rand);
  double ***Sigma = doubleMatrix3D(*tmax, *n_rand, *n_rand);
  double ***SigInv = doubleMatrix3D(*tmax, *n_rand, *n_rand);
  double **zeta = doubleMatrix(*n_grp, *n_rand);  /* control */
  double *zeta0 = doubleArray(*n_rand);
  double **Phi = doubleMatrix(*n_rand, *n_rand);
  double **PhiInv = doubleMatrix(*n_rand, *n_rand);

  itemp = 0;
  for (k = 0; k < *tmax; k++)
    for (j = 0; j < *n_rand; j++)
      for (i = 0; i < *n_rand; i++)
	Sigma[k][i][j] = SigmaAll[itemp++];
  for (k = 0; k < *tmax; k++)
    dinv(Sigma[k], *n_rand, SigInv[k]);
  
  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      Phi[i][j] = PhiAll[itemp++];
  dinv(Phi, *n_rand, PhiInv);

  /* starting values for random effects from prior */
  for (j = 0; j < *n_rand; j++) {
    gamma0[j] = 0; 
    zeta0[j] = 0;
  }
  for (i = 0; i < *tmax; i++) {
    for (j = 0; j < *n_grp; j++) 
      rMVN(gamma[i][j], gamma0, Sigma[i], *n_rand);
  }
  for (j = 0; j < *n_grp; j++) {
    rMVN(zeta[j], zeta0, Phi, *n_rand);
  }

  /** Data **/
  int **Zstar = intMatrix(*tmax, *n_samp);
  int *Y0 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_dim);
  double ***Xtemp = doubleMatrix3D(*tmax, *n_samp, n_dim);
  double *Xdelta = doubleArray(*n_samp);
  double *Xdelta1 = doubleArray(*n_samp);
  double *Xpsi = doubleArray(*n_samp);
  int **gtmp = intMatrix(*tmax, *n_samp);
  double ****Vgtmp = doubleMatrix4D(*tmax, *n_grp, *max_samp_grp, *n_rand); 
  int **Mitemp = intMatrix(*tmax, *n_grp); 
  double ***Vgrp = doubleMatrix3D(*n_grp, *max_samp_grp, *n_rand);
  
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = Xall[itemp++];
  
  itemp = 0;
  for (j = 0; j < *n_grp; j++) 
    vitemp[j] = 0; 
  for (i = 0; i < *tmax; i++) {
    treatSum[i] = 0;
    for (j = 0; j < *n_grp; j++) 
      Mitemp[i][j] = 0; 
  }  
  for (i = 0; i < *n_samp; i++) {
    Y0[i] = Y[i];
    for (j = 0; j < *n_rand; j++)
      Vgrp[grp[i]][vitemp[grp[i]]][j] = Xall_rand[itemp++];
    if (treat[i] > 0) {
      gtmp[treat[i]-1][treatSum[treat[i]-1]] = grp[i];
      for (j = 0; j < *n_cov; j++) 
	Xtemp[treat[i]-1][treatSum[treat[i]-1]][j] = X[i][j];
      for (j = 0; j < *n_rand; j++)
	Vgtmp[treat[i]-1][grp[i]][Mitemp[treat[i]-1][grp[i]]][j] = Vgrp[grp[i]][vitemp[grp[i]]][j];
      treatSum[treat[i]-1]++;
      Mitemp[treat[i]-1][grp[i]]++;
      if ((ceiling[treat[i]-1] == 1) && (Y[i] == (*J+1))) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");	
      }
      if ((floor[treat[i]-1] == 1) && (Y[i] == 0)) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    }
    vitemp[grp[i]]++;
  }

  /** Prior **/
  double ***A0delta = doubleMatrix3D(*tmax, n_dim, n_dim);
  double **A0psi = doubleMatrix(*n_cov, *n_cov);
  double ***S0 = doubleMatrix3D(*tmax, *n_rand, *n_rand);
  double **T0 = doubleMatrix(*n_rand, *n_rand);

  itemp = 0;
  for (k = 0; k < *tmax; k++)
    for (j = 0; j < n_dim; j++)
      for (i = 0; i < n_dim; i++)
	A0delta[k][i][j] = A0deltaAll[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      A0psi[i][j] = A0psiAll[itemp++];

  itemp = 0;
  for (k = 0; k < *tmax; k++)
    for (j = 0; j < *n_rand; j++)
      for (i = 0; i < *n_rand; i++)
	S0[k][i][j] = S0All[itemp++];

  itemp = 0;
  for (j = 0; j < *n_rand; j++)
    for (i = 0; i < *n_rand; i++)
      T0[i][j] = T0All[itemp++];

  /** Proposal precisoin **/
  double ***deltaPro = doubleMatrix3D(*tmax, n_dim, n_dim);
  double **psiPro = doubleMatrix(*n_cov, *n_cov);

  itemp = 0;
  for (k = 0; k < *tmax; k++)
    for (j = 0; j < n_dim; j++)
      for (i = 0; i < n_dim; i++)
	deltaPro[k][i][j] = deltaVar[itemp++];

  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_cov; i++)
      psiPro[i][j] = psiVar[itemp++];

  /** Initializing the counters **/
  itemp = 0; psiCounter[0] = 0; 
  for (i = 0; i < *tmax; i++) {
    deltaCounter[i][0] = 0; 
    for (j = 0; j < *n_grp; j++) 
      gammaCounter[i][j] = 0;
  }
  for (j = 0; j < *n_grp; j++) 
    zetaCounter[j] = 0;

  /** MCMC **/
  if (*verbose) {
    Rprintf("\n *** Starting posterior sampling... ***\n");
  }
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *n_grp; j++)
      vitemp[j] = 0;
    for (j = 0; j < *tmax; j++)
      treatSum[j] = 0;
    for (i = 0; i < *n_samp; i++) {
      /* Sample Zstar */
      if (treat[i] > 0) { 
	if (Y[i] == (*J+1)) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0;  Xpsi[i] = 0;  
	  for (j = 0; j < *n_cov; j++) { /* fixed effects */
	    Xdelta[i] += X[i][j]*deltaMatrix[treat[i]-1][j];
	    Xpsi[i] += X[i][j]*psi[j];
	  }
	  for (j = 0; j < *n_rand; j++) { /* random effects */
	    Xdelta[i] += Vgrp[grp[i]][vitemp[grp[i]]][j] * gamma[treat[i]-1][grp[i]][j];
	    Xpsi[i] += Vgrp[grp[i]][vitemp[grp[i]]][j] * zeta[grp[i]][j];
	  }
	  if ((ceiling[treat[i]-1] == 1) && (Y[i] == *J)) { /* ceiling effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J-1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp2 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J-1;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = *J;
	    }
	  } else if ((floor[treat[i]-1] == 1) && (Y[i] == 1)) { /* floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(1, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xpsi[i])), 1)); 
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = 0;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 1;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 0;
	    }
	  } else { /* no ceiling and floor effects */
	    if (*unconst) {
	      Xdelta1[i] = Xdelta[i] + (Y[i] - 1) * deltaMatrix[treat[i]-1][*n_cov];
	      Xdelta[i] += Y[i] * deltaMatrix[treat[i]-1][*n_cov];
	      dtemp1 = exp(Xdelta1[i] - log1p(exp(Xdelta1[i])) + 
			   dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    } else {
	      dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			   dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    }
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xpsi[i])), 1));
	    if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = Y[i] - 1;
	    } else { 
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = Y[i];
	    }
	  }
	}
	if (*unconst) {
	  X[i][*n_cov] = Y0[i];
	  Xtemp[treat[i]-1][treatSum[treat[i]-1]][*n_cov] = Y0[i];
	}
	treatSum[treat[i]-1]++;
      }
      vitemp[grp[i]]++;
    }

    /* Sample delta */
    for (k = 0; k < *tmax; k++) {
      BinomLogitMixed(Zstar[k], Xtemp[k], Vgtmp[k], gtmp[k], deltaMatrix[k], gamma[k], SigInv[k], treatSum[k], 1,
		      n_dim, *n_rand, *n_grp, delta0Matrix[k], A0delta[k], s0[k], S0[k], deltaPro[k], 
		      tune_gamma, 1, deltaCounter[k], gammaCounter[k]);
    } 

    /* Sample psi */
    BinomLogitMixed(Y0, X, Vgrp, grp, psi, zeta, PhiInv, *n_samp, *J, *n_cov, *n_rand, *n_grp, 
		    psi0, A0psi, *t0, T0, psiPro, tune_zeta, 1, psiCounter, zetaCounter);
    
    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	/* fixed effects */
	for (i = 0; i < *tmax; i++) {
	  for (j = 0; j < n_dim; j++) {
	    allresults[itemp++] = deltaMatrix[i][j];
	  }      
	}
	for (j = 0; j < *n_cov; j++) {
	  allresults[itemp++] = psi[j];
	}
	/* Var-Cov matrices */
	for (i = 0; i < *tmax; i++) {
	  dinv(SigInv[i], *n_rand, Sigma[i]);
	  for (j = 0; j < *n_rand; j++) {
	    for (k = j; k < *n_rand; k++) {
	      allresults[itemp++] = Sigma[i][j][k];
	    }      
	  }
	}
	dinv(PhiInv, *n_rand, Phi);
	for (j = 0; j < *n_rand; j++) {
	  for (k = j; k < *n_rand; k++) {
	    allresults[itemp++] = Phi[j][k];
	  }      
	}
	/* random effects */
	for (i = 0; i < *tmax; i++) {
	  for (j = 0; j < *n_grp; j++) {
	    for (k = 0; k < *n_rand; k++) {
	      allresults[itemp++] = gamma[i][j][k];
	    }
	  }
	}
	for (j = 0; j < *n_grp; j++) {
	  for (k = 0; k < *n_rand; k++) {
	    allresults[itemp++] = zeta[j][k];
	  }
	}
	/* acceptance ratios */
	for (i = 0; i < *tmax; i++) 
	  allresults[itemp++] = ((double) deltaCounter[i][0] / (double) (main_loop + 1));
	allresults[itemp++] = ((double) *psiCounter / (double) (main_loop + 1));
	for (i = 0; i < *tmax; i++) 
	  for (j = 0; j < *n_grp; j++)
	    allresults[itemp++] = ((double) gammaCounter[i][j] / (double) (main_loop + 1));
	for (j = 0; j < *n_grp; j++)
	  allresults[itemp++] = ((double) zetaCounter[j] / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("\n%3d percent done.\n    Metropolis acceptance ratios\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("      Sensitive item model (fixed effects):");
	for (i = 0; i < *tmax; i++) {
	  Rprintf(" %3g", fprec((double) deltaCounter[i][0] / (double) (main_loop + 1), 3));
	}
	Rprintf("\n      Control items model (fixed effects): %3g",
		fprec((double) *psiCounter / (double) (main_loop + 1), 3));
	for (i = 0; i < *tmax; i++) {
	  Rprintf("\n      Sensitive item model (random effects for each group): %2d):", 
		  i + 1);
	  for (j = 0; j < *n_grp; j++) 
	    Rprintf(" %3g", fprec((double) gammaCounter[i][j] / (double) (main_loop + 1), 3));
	}
	Rprintf("\n      Control item model (random effects for each group):");
	for (j = 0; j < *n_grp; j++) 
	  Rprintf(" %3g", fprec((double) zetaCounter[j] / (double) (main_loop + 1), 3));
	Rprintf("\n");
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(Y0);
  free(Xdelta);
  free(Xdelta1);
  free(Xpsi);
  free(psiCounter);
  free(zetaCounter);
  free(vitemp);
  free(treatSum);
  free(gamma0);
  free(zeta0);
  FreeintMatrix(gtmp, *tmax);
  FreeintMatrix(Mitemp, *tmax); 
  FreeintMatrix(Zstar, *tmax);
  FreeMatrix(zeta, *n_grp);
  Free3DMatrix(gamma, *tmax, *n_grp);
  Free4DMatrix(Vgtmp, *tmax, *n_grp, *max_samp_grp); 
  FreeintMatrix(deltaCounter, *tmax);
  FreeintMatrix(gammaCounter, *tmax);
  FreeMatrix(deltaMatrix, *tmax);
  FreeMatrix(delta0Matrix, *tmax);
  Free3DMatrix(Sigma, *tmax, *n_rand);
  Free3DMatrix(SigInv, *tmax, *n_rand);
  FreeMatrix(Phi, *n_rand);
  FreeMatrix(PhiInv, *n_rand);
  FreeMatrix(X, *n_samp);
  Free3DMatrix(Xtemp, *tmax, *n_samp);
  Free3DMatrix(A0delta, *tmax, *n_cov);
  FreeMatrix(A0psi, *n_cov);
  Free3DMatrix(S0, *tmax, *n_rand);
  FreeMatrix(T0, *n_rand);
  Free3DMatrix(deltaPro, *tmax, *n_cov);
  FreeMatrix(psiPro, *n_cov);
  Free3DMatrix(Vgrp, *n_grp, *max_samp_grp);
}



/** 
  Item Count Technique Binomial Mixed Effects Regression for the Multiple Sensitive Item Design
  
  This is the 2-level multilevel model 

**/

void ictregBinomMulti2Level(int *Y,             /* outcome vector */
			    int *J,             /* # of control items */
			    int *n_samp,        /* sample size */
			    int *n_draws,       /* # of MCMC draws */
			    int *treat,         /* treatment indicator vector: 0, 1, or 2 */
			    int *tmax,          /* number of sensitive items */
			    double *Xall,       /* individual-level covariates in a vector form */
			    double *Vall,       /* village-level covariates in a vector form */
			    double *beta,       /* individual-level coefs for control and sensitive items */
			    double *beta_village, /* village-level coefs for control and sensitive items */
			    int *n_cov,         /* # of individual-level covariates */
			    int *n_covV,        /* # of village-level covariates */
			    double *beta0,         /* prior mean for beta */
			    double *beta0_village, /* prior mean for beta_village */
			    double *A0betaAll,  /* prior precision for beta */
			    double *A0beta_villageAll,  /* prior precision for beta_village */
			    double *betaPro,    /* proposal precision for beta */
			    int *ceiling,       /* ceiling effects */
			    int *floor,         /* floor effects */
			    int *village,       /* village indicator, 0, 1, ..., G-1 */
			    int *n_village,     /* number of villages, G */
			    double *alphaPro,   /* proposal precision for beta */
			    double *sigma2All,  /* variances for random effects  */
			    int *nu0,           /* prior df for sigma */
			    double *s0,         /* prior scale for sigma */
			    int *burnin,        /* number of burnins */
			    int *keep,          /* keep every *th draw */
			    int *verbose,       /* want to print progress? */
			    double *allresults  /* storage for all results */
			    ) {
  
  int i, j, k, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  int n_grp = *tmax + 1;
  int *vitemp = intArray(*n_village);
  int *treatSum = intArray(*tmax);
  double dtemp, dtemp1, dtemp2, dtemp3;

  /** get random seed **/
  GetRNGstate();

  /** Individual level Data **/
  int n_dim = *n_cov + *n_village;
  int **Zstar = intMatrix(*tmax, *n_samp);
  int *Y0 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_dim); /* village dummies plus individual covariates */
  double ***Xtemp = doubleMatrix3D(*tmax, *n_samp, n_dim);
  double *Xdelta = doubleArray(*n_samp);
  double *Xdelta0 = doubleArray(*n_samp);
  int **gtmp = intMatrix(*tmax, *n_samp);
  
  itemp = 0;
  for (j = *n_village; j < n_dim; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = Xall[itemp++];
  
  for (j = 0; j < *n_village; j++) 
    vitemp[j] = 0; 
  for (i = 0; i < *tmax; i++)
    treatSum[i] = 0;

  itemp = 0;
  for (i = 0; i < *n_samp; i++) {
    Y0[i] = Y[i];
    /* village dummies */
    for (j = 0; j < *n_village; j++)
      X[i][j] = 0;
    X[i][village[i]] = 1;
    if (treat[i] > 0) {
      /* village indicator for each treatment group */
      gtmp[treat[i]-1][treatSum[treat[i]-1]] = village[i];
      /* X matrix for each treatment group */
      for (j = 0; j < n_dim; j++) 
	Xtemp[treat[i]-1][treatSum[treat[i]-1]][j] = X[i][j];
      treatSum[treat[i]-1]++;
      if ((ceiling[treat[i]-1] == 1) && (Y[i] == (*J+1))) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");	
      }
      if ((floor[treat[i]-1] == 1) && (Y[i] == 0)) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    }
    vitemp[village[i]]++;
  }

  /** Village level data **/
  double **V = doubleMatrix(*n_village, *n_covV);

  itemp = 0;
  for (j = 0; j < *n_covV; j++)
    for (i = 0; i < *n_village; i++) 
      V[i][j] = Vall[itemp++];

  /** parameters for fixed and random effects with starting values **/
  double **delta = doubleMatrix(n_grp, n_dim);
  double **delta_village = doubleMatrix(n_grp, *n_covV);
  double **sigma2_village = doubleMatrix(n_grp, 1);

  itemp = 0;
  for (i = 0; i < n_grp; i++)
    sigma2_village[i][0] = sigma2All[itemp++];
  /* PdoubleMatrix(sigma2_village, n_grp, 1); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = 0; j < *n_covV; j++)
      delta_village[i][j] = beta_village[itemp++];
  }
  /* PdoubleMatrix(delta_village, n_grp, *n_covV); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = *n_village; j < n_dim; j++)
      delta[i][j] = beta[itemp++];
    for (j = 0; j < *n_village; j++) {
      delta[i][j] = norm_rand() * sqrt(sigma2_village[i][0]);
      for (k = 0; k < *n_covV; k++)
	delta[i][j] += delta_village[i][k] * V[j][k];
    }
  }
  /* PdoubleMatrix(delta, n_grp, n_dim); */

  /** Prior **/
  double **delta0 = doubleMatrix(n_grp, n_dim);
  double ***A0delta = doubleMatrix3D(n_grp, n_dim, n_dim);
  double **delta0_village = doubleMatrix(n_grp, *n_covV);
  double ***A0delta_village = doubleMatrix3D(n_grp, *n_covV, *n_covV);

  itemp = 0; 
  for (i = 0; i < n_grp; i++) { 
    for (j = *n_village; j < n_dim; j++)
      delta0[i][j] = beta0[itemp++];
    for (j = 0; j < *n_village; j++) {
      delta0[i][j] = 0;
      for (k = 0; k < *n_covV; k++)
	delta0[i][j] += delta_village[i][k] * V[j][k];
    }  
  }
  /* PdoubleMatrix(delta0, n_grp, n_dim); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {     
    for (j = 0; j < *n_covV; j++)
      delta0_village[i][j] = beta0_village[itemp++];
  }
  /* PdoubleMatrix(delta0_village, n_grp, *n_covV); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = *n_village; j < n_dim; j++) {
      for (i = *n_village; i < n_dim; i++)
	A0delta[k][i][j] = A0betaAll[itemp++];
      for (i = 0; i < *n_village; i++)
	A0delta[k][i][j] = 0;
    }
    for (j = 0; j < *n_village; j++) {
      for (i = 0; i < n_dim; i++) {
	if (i == j) 
	  A0delta[k][i][j] = sigma2_village[k][0];
	else 
	  A0delta[k][i][j] = 0;
      } 
    }
  }
  /* PdoubleMatrix3D(A0delta, n_grp, n_dim, n_dim); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = 0; j < *n_covV; j++) {
      for (i = 0; i < *n_covV; i++)
	A0delta_village[k][i][j] = A0beta_villageAll[itemp++];
    }
  }
  /* PdoubleMatrix3D(A0delta_village, n_grp, *n_covV, *n_covV); */

  /** Proposal precisoin **/
  double ***deltaPro = doubleMatrix3D(n_grp, n_dim, n_dim);

  itemp = 0;
  for (k = 0; k < n_grp; k++)
    for (j = 0; j < *n_village; j++) {
      for (i = 0; i < *n_village; i++)
	if (i == j) 
	  deltaPro[k][i][j] = alphaPro[itemp++];
	else 
	  deltaPro[k][i][j] = 0;
      for (i = *n_village; i < n_dim; i++)
	deltaPro[k][i][j] = 0;
    }

  itemp = 0;
  for (k = 0; k < n_grp; k++)
    for (j = *n_village; j < n_dim; j++) {
      for (i = *n_village; i < n_dim; i++) 
	if (i == j) 
	  deltaPro[k][i][j] = betaPro[itemp++];
	else
	  deltaPro[k][i][j] = 0;
      for (i = 0; i < *n_village; i++) 
	deltaPro[k][i][j] = 0;
    }
  /* PdoubleMatrix3D(deltaPro, n_grp, n_dim, n_dim); */

  /** Counter to calculate acceptance ratio for delta */
  int **deltaCounter = intMatrix(n_grp, 1);  

  itemp = 0; 
  for (i = 0; i < n_grp; i++) 
    deltaCounter[i][0] = 0; 

  /** MCMC **/
  if (*verbose) 
    Rprintf("\n *** Starting posterior sampling... ***\n");
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *n_village; j++)
      vitemp[j] = 0;
    for (j = 0; j < *tmax; j++)
      treatSum[j] = 0;
    /* Sample Zstar */
    for (i = 0; i < *n_samp; i++) {
      if (treat[i] > 0) { 
	if (Y[i] == (*J+1)) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0; Xdelta0[i] = 0;
	  for (j = 0; j < n_dim; j++) { 
	    Xdelta[i] += X[i][j]*delta[treat[i]][j];
	    Xdelta0[i] += X[i][j]*delta[0][j];
	  }
	  if ((ceiling[treat[i]-1] == 1) && (Y[i] == *J)) { /* ceiling effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J-1, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp2 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J-1;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = *J;
	    }
	  } else if ((floor[treat[i]-1] == 1) && (Y[i] == 1)) { /* floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(1, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = 0;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 1;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 0;
	    }
	  } else { /* no ceiling and floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = Y[i] - 1;
	    } else { 
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = Y[i];
	    }
	  }
	}
	treatSum[treat[i]-1]++;
      }
      vitemp[village[i]]++;
    }

    /* Sample delta for sensitive items */
    for (k = 1; k < n_grp; k++) {
      BinomLogit(Zstar[k-1], Xtemp[k-1], delta[k], treatSum[k-1], 1, n_dim, delta0[k], A0delta[k], deltaPro[k], 1, 
		 deltaCounter[k]);
    } 

    /* Sample delta for control items */
    BinomLogit(Y0, X, delta[0], *n_samp, *J, n_dim, delta0[0], A0delta[0], deltaPro[0], 1, deltaCounter[0]);
    
    /* Update the village level model */
    for (k = 0; k < n_grp; k++) {
      bNormalReg(delta[k], V, delta_village[k], sigma2_village[k], *n_village, *n_covV, 1, 1, delta0_village[k], A0delta_village[k], 
		 1, s0[k], nu0[k], 0, 0);
    }

    /* Update the prior information */
    for (i = 0; i < n_grp; i++) { 
      for (j = 0; j < *n_village; j++) {
	delta0[i][j] = 0;
	for (k = 0; k < *n_covV; k++)
	  delta0[i][j] += delta_village[i][k] * V[j][k];
      }  
    }

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	/* fixed and random effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = *n_village; j < n_dim; j++) {
	    allresults[itemp++] = delta[i][j];
	  }      
	}
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_village; j++) {
	    allresults[itemp++] = delta[i][j];
	  }      
	}
	/* village level coefficients */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_covV; j++) {
	    allresults[itemp++] = delta_village[i][j];
	  }      
	}
	/* sigmas */
	for (i = 0; i < n_grp; i++) {
	  allresults[itemp++] = sigma2_village[i][0];
	}
	/* acceptance ratios */
	for (i = 0; i < n_grp; i++) 
	  allresults[itemp++] = ((double) deltaCounter[i][0] / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("\n%3d percent done.\n    Metropolis acceptance ratios\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("      Control item model:");
	Rprintf(" %3g", fprec((double) deltaCounter[0][0] / (double) (main_loop + 1), 3));
	Rprintf("\n");
	for (i = 1; i < n_grp; i++) {
	  Rprintf("      Sensitive item model %3d:", i);
	  Rprintf(" %3g", fprec((double) deltaCounter[i][0] / (double) (main_loop + 1), 3));
	  Rprintf("\n");
	}
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(vitemp);
  free(treatSum);
  FreeintMatrix(Zstar, *tmax);
  free(Y0);
  FreeMatrix(X, *n_samp);
  Free3DMatrix(Xtemp, *tmax, *n_samp);
  free(Xdelta);
  free(Xdelta0);
  FreeintMatrix(gtmp, *tmax);
  FreeMatrix(V, *n_village);
  FreeMatrix(delta, n_grp);
  FreeMatrix(delta_village, n_grp);
  FreeMatrix(sigma2_village, n_grp);
  FreeMatrix(delta0, n_grp);
  Free3DMatrix(A0delta, n_grp, n_dim);
  FreeMatrix(delta0_village, n_grp);
  Free3DMatrix(A0delta_village, n_grp, *n_covV);
  Free3DMatrix(deltaPro, n_grp, n_dim);
  FreeintMatrix(deltaCounter, n_grp);

}



/** 
  Item Count Technique Binomial Mixed Effects Regression for the Multiple Sensitive Item Design
  
  This is the 3-level multilevel model 

**/

void ictregBinomMulti3Level(int *Y,             /* outcome vector */
			    int *J,             /* # of control items */
			    int *n_samp,        /* sample size */
			    int *n_draws,       /* # of MCMC draws */
			    int *treat,         /* treatment indicator vector: 0, 1, or 2 */
			    int *tmax,          /* number of sensitive items */
			    double *Xall,       /* individual-level covariates in a vector form */
			    double *Vall,       /* village-level covariates in a vector form */
			    double *Wall,       /* district-level covariates in a vector form */
			    double *beta,       /* individual-level coefs for control and sensitive items */
			    double *beta_village,  /* village-level coefs for control and sensitive items */
			    double *beta_district, /* district-level coefs for control and sensitive items */
			    int *n_cov,         /* # of individual-level covariates */
			    int *n_covV,        /* # of village-level covariates */
			    int *n_covW,        /* # of district-level covariates */
			    double *beta0,         /* prior mean for beta */
			    double *beta0_village, /* prior mean for beta_village */
			    double *beta0_district, /* prior mean for beta_district */
			    double *A0betaAll,  /* prior precision for beta */
			    double *A0beta_villageAll,  /* prior precision for beta_village */
			    double *A0beta_districtAll,  /* prior precision for beta_district */
			    double *betaPro,    /* proposal precision for beta */
			    int *ceiling,       /* ceiling effects */
			    int *floor,         /* floor effects */
			    int *village,       /* village indicator, 0, 1, ..., G-1 */
			    int *n_village,     /* number of villages, G */
			    int *district,      /* district indicator, 0, 1, ..., D-1 */
			    int *n_district,    /* number of districts, D */
			    double *alphaPro,   /* proposal precision for beta */
			    double *sigma2_villageAll,  /* variances for village level random effects  */
			    double *sigma2_districtAll, /* variances for district level random effects  */
			    int *nu0_village,           /* prior df for sigma2_village */
			    double *s0_village,         /* prior scale for sigma2_vilalge */
			    int *nu0_district,          /* prior df for sigma2_district */
			    double *s0_district,        /* prior scale for sigma2_district */
			    int *burnin,        /* number of burnins */
			    int *keep,          /* keep every *th draw */
			    int *verbose,       /* want to print progress? */
			    double *allresults  /* storage for all results */
			    ) {
  
  int i, j, k, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  double dtemp, dtemp1, dtemp2, dtemp3;

  /** get random seed **/
  GetRNGstate();

  /** Individual level Data **/
  int n_grp = *tmax + 1;
  int n_dim = *n_cov + *n_village;
  int **Zstar = intMatrix(*tmax, *n_samp);
  int *Y0 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_dim); /* village dummies plus individual covariates */
  double ***Xtemp = doubleMatrix3D(*tmax, *n_samp, n_dim);
  double *Xdelta = doubleArray(*n_samp);
  double *Xdelta0 = doubleArray(*n_samp);
  int **gtmp = intMatrix(*tmax, *n_samp);
  int *vitemp = intArray(*n_village);
  int *treatSum = intArray(*tmax);
  
  itemp = 0;
  for (j = *n_village; j < n_dim; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = Xall[itemp++];
  
  for (j = 0; j < *n_village; j++) 
    vitemp[j] = 0; 
  for (i = 0; i < *tmax; i++)
    treatSum[i] = 0;

  itemp = 0;
  for (i = 0; i < *n_samp; i++) {
    Y0[i] = Y[i];
    /* village dummies */
    for (j = 0; j < *n_village; j++)
      X[i][j] = 0;
    X[i][village[i]] = 1;
    if (treat[i] > 0) {
      /* village indicator for each treatment group */
      gtmp[treat[i]-1][treatSum[treat[i]-1]] = village[i];
      /* X matrix for each treatment group */
      for (j = 0; j < n_dim; j++) 
	Xtemp[treat[i]-1][treatSum[treat[i]-1]][j] = X[i][j];
      treatSum[treat[i]-1]++;
      if ((ceiling[treat[i]-1] == 1) && (Y[i] == (*J+1))) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");	
      }
      if ((floor[treat[i]-1] == 1) && (Y[i] == 0)) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    }
    vitemp[village[i]]++;
  }
  /* PdoubleMatrix(X, *n_samp, n_dim); */

  /** Village level data **/
  int n_dimV = *n_covV + *n_district;
  int *vitemp1 = intArray(*n_district);
  double **V = doubleMatrix(*n_village, n_dimV); /* district dummies plus village level covariates */

  for (i = 0; i < *n_district; i++) 
    vitemp1[i] = 0;

  itemp = 0;
  for (j = *n_district; j < n_dimV; j++)
    for (i = 0; i < *n_village; i++) 
      V[i][j] = Vall[itemp++];

  for (i = 0; i < *n_village; i++) {
    for (j = 0; j < *n_district; j++)
      V[i][j] = 0;
    V[i][district[i]] = 1;
    vitemp1[district[i]]++;
  }
  /* PdoubleMatrix(V, *n_village, n_dimV); */

  /** District level data **/
  double **W = doubleMatrix(*n_district, *n_covW);

  itemp = 0;
  for (j = 0; j < *n_covW; j++)
    for (i = 0; i < *n_district; i++) 
      W[i][j] = Wall[itemp++];
  /* PdoubleMatrix(W, *n_district, *n_covW); */

  /** parameters for fixed and random effects with starting values **/
  double **delta = doubleMatrix(n_grp, n_dim);
  double **delta_village = doubleMatrix(n_grp, n_dimV);
  double **delta_district = doubleMatrix(n_grp, *n_covW);
  double **sigma2_village = doubleMatrix(n_grp, 1);
  double **sigma2_district = doubleMatrix(n_grp, 1);

  itemp = 0;
  for (i = 0; i < n_grp; i++)
    sigma2_village[i][0] = sigma2_villageAll[itemp++];
  /* PdoubleMatrix(sigma2_village, n_grp, 1); */

  itemp = 0;
  for (i = 0; i < n_grp; i++)
    sigma2_district[i][0] = sigma2_districtAll[itemp++];

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = 0; j < *n_covW; j++)
      delta_district[i][j] = beta_district[itemp++];
  }
  /* PdoubleMatrix(delta_district, n_grp, *n_covW); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = *n_district; j < n_dimV; j++)
      delta_village[i][j] = beta_village[itemp++];
    for (j = 0; j < *n_district; j++) {
      delta_village[i][j] = norm_rand() * sqrt(sigma2_district[i][0]);
      for (k = 0; k < *n_covW; k++)
	delta_village[i][j] += delta_district[i][k] * W[j][k];
    }
  }
  /* PdoubleMatrix(delta_village, n_grp, n_dimV); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = *n_village; j < n_dim; j++)
      delta[i][j] = beta[itemp++];
    for (j = 0; j < *n_village; j++) {
      delta[i][j] = norm_rand() * sqrt(sigma2_village[i][0]);
      for (k = 0; k < n_dimV; k++)
	delta[i][j] += delta_village[i][k] * V[j][k];
    }
  }
  /* PdoubleMatrix(delta, n_grp, n_dim); */


  /** Prior **/
  double **delta0 = doubleMatrix(n_grp, n_dim);
  double ***A0delta = doubleMatrix3D(n_grp, n_dim, n_dim);
  double **delta0_village = doubleMatrix(n_grp, n_dimV);
  double ***A0delta_village = doubleMatrix3D(n_grp, n_dimV, n_dimV);
  double **delta0_district = doubleMatrix(n_grp, *n_covW);
  double ***A0delta_district = doubleMatrix3D(n_grp, *n_covW, *n_covW);

  itemp = 0; 
  for (i = 0; i < n_grp; i++) { 
    for (j = *n_village; j < n_dim; j++)
      delta0[i][j] = beta0[itemp++];
    for (j = 0; j < *n_village; j++) {
      delta0[i][j] = 0;
      for (k = 0; k < n_dimV; k++)
	delta0[i][j] += delta_village[i][k] * V[j][k];
    }  
  }
  /* PdoubleMatrix(delta0, n_grp, n_dim); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) { 
    for (j = *n_district; j < n_dimV; j++)
      delta0_village[i][j] = beta0_village[itemp++];
    for (j = 0; j < *n_district; j++) {
      delta0_village[i][j] = 0;
      for (k = 0; k < *n_covW; k++)
	delta0_village[i][j] += delta_district[i][k] * W[j][k];
    }  
  }
  /* PdoubleMatrix(delta0, n_grp, n_dim); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {     
    for (j = 0; j < *n_covW; j++)
      delta0_district[i][j] = beta0_district[itemp++];
  }
  /* PdoubleMatrix(delta0_district, n_grp, *n_covW); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = *n_village; j < n_dim; j++) {
      for (i = *n_village; i < n_dim; i++)
	A0delta[k][i][j] = A0betaAll[itemp++];
      for (i = 0; i < *n_village; i++)
	A0delta[k][i][j] = 0;
    }
    for (j = 0; j < *n_village; j++) {
      for (i = 0; i < n_dim; i++) {
	if (i == j) 
	  A0delta[k][i][j] = sigma2_village[k][0];
	else 
	  A0delta[k][i][j] = 0;
      } 
    }
  }
  /* PdoubleMatrix3D(A0delta, n_grp, n_dim, n_dim); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = *n_district; j < n_dimV; j++) {
      for (i = *n_district; i < n_dimV; i++)
	A0delta_village[k][i][j] = A0beta_villageAll[itemp++];
      for (i = 0; i < *n_district; i++)
	A0delta_village[k][i][j] = 0;
    }
    for (j = 0; j < *n_district; j++) {
      for (i = 0; i < n_dimV; i++) {
	if (i == j) 
	  A0delta_village[k][i][j] = sigma2_district[k][0];
	else 
	  A0delta_village[k][i][j] = 0;
      } 
    }
  }
  /* PdoubleMatrix3D(A0delta, n_grp, n_dim, n_dim); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = 0; j < *n_covW; j++) {
      for (i = 0; i < *n_covW; i++)
	A0delta_district[k][i][j] = A0beta_districtAll[itemp++];
    }
  }
  /* PdoubleMatrix3D(A0delta_district, n_grp, *n_covW, *n_covW); */

  /** Proposal precisoin **/
  double ***deltaPro = doubleMatrix3D(n_grp, n_dim, n_dim);

  itemp = 0;
  for (k = 0; k < n_grp; k++)
    for (j = 0; j < *n_village; j++) {
      for (i = 0; i < *n_village; i++)
	if (i == j) 
	  deltaPro[k][i][j] = alphaPro[itemp++];
	else 
	  deltaPro[k][i][j] = 0;
      for (i = *n_village; i < n_dim; i++)
	deltaPro[k][i][j] = 0;
    }

  itemp = 0;
  for (k = 0; k < n_grp; k++)
    for (j = *n_village; j < n_dim; j++) {
      for (i = *n_village; i < n_dim; i++) 
	if (i == j) 
	  deltaPro[k][i][j] = betaPro[itemp++];
	else
	  deltaPro[k][i][j] = 0;
      for (i = 0; i < *n_village; i++) 
	deltaPro[k][i][j] = 0;
    }
  /* PdoubleMatrix3D(deltaPro, n_grp, n_dim, n_dim); */

  /** Counter to calculate acceptance ratio for delta */
  int **deltaCounter = intMatrix(n_grp, 1);  

  itemp = 0; 
  for (i = 0; i < n_grp; i++) 
    deltaCounter[i][0] = 0; 

  /** MCMC **/
  if (*verbose) 
    Rprintf("\n *** Starting posterior sampling... ***\n");
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *n_village; j++)
      vitemp[j] = 0;
    for (j = 0; j < *tmax; j++)
      treatSum[j] = 0;
    /* Sample Zstar */
    for (i = 0; i < *n_samp; i++) {
      if (treat[i] > 0) { 
	if (Y[i] == (*J+1)) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0; Xdelta0[i] = 0;
	  for (j = 0; j < n_dim; j++) { 
	    Xdelta[i] += X[i][j]*delta[treat[i]][j];
	    Xdelta0[i] += X[i][j]*delta[0][j];
	  }
	  if ((ceiling[treat[i]-1] == 1) && (Y[i] == *J)) { /* ceiling effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J-1, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp2 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J-1;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = *J;
	    }
	  } else if ((floor[treat[i]-1] == 1) && (Y[i] == 1)) { /* floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(1, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = 0;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 1;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 0;
	    }
	  } else { /* no ceiling and floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = Y[i] - 1;
	    } else { 
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = Y[i];
	    }
	  }
	}
	treatSum[treat[i]-1]++;
      }
      vitemp[village[i]]++;
    }

    /* Sample delta for sensitive items */
    for (k = 1; k < n_grp; k++) {
      BinomLogit(Zstar[k-1], Xtemp[k-1], delta[k], treatSum[k-1], 1, n_dim, delta0[k], A0delta[k], deltaPro[k], 1, 
		 deltaCounter[k]);
    } 

    /* Sample delta for control items */
    BinomLogit(Y0, X, delta[0], *n_samp, *J, n_dim, delta0[0], A0delta[0], deltaPro[0], 1, deltaCounter[0]);
    
    /* Update the village level model */
    for (k = 0; k < n_grp; k++) {
      bNormalReg(delta[k], V, delta_village[k], sigma2_village[k], *n_village, n_dimV, 1, 1, delta0_village[k], A0delta_village[k], 
		 1, s0_village[k], nu0_village[k], 0, 0);
    }

    /* Update the district level model */
    for (k = 0; k < n_grp; k++) {
      bNormalReg(delta_village[k], W, delta_district[k], sigma2_district[k], *n_district, *n_covW, 1, 1, delta0_district[k], 
		 A0delta_district[k], 1, s0_district[k], nu0_district[k], 0, 0);
    }

    /* Update the prior information */
    for (i = 0; i < n_grp; i++) { 
      for (j = 0; j < *n_village; j++) {
	delta0[i][j] = 0;
	for (k = 0; k < n_dimV; k++)
	  delta0[i][j] += delta_village[i][k] * V[j][k];
      } 
      for (j = 0; j < *n_district; j++) {
	delta0_village[i][j] = 0;
	for (k = 0; k < *n_covW; k++)
	  delta0_village[i][j] += delta_district[i][k] * W[j][k];
      } 
    }

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	/* individual level fixed effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = *n_village; j < n_dim; j++) {
	    allresults[itemp++] = delta[i][j];
	  }      
	}
	/* individual level random effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_village; j++) {
	    allresults[itemp++] = delta[i][j];
	  }      
	}
	/* village level fixed effects coefficients */
	for (i = 0; i < n_grp; i++) {
	  for (j = *n_district; j < n_dimV; j++) {
	    allresults[itemp++] = delta_village[i][j];
	  }      
	}
	/* village level random effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_district; j++) {
	    allresults[itemp++] = delta_village[i][j];
	  }      
	}
	/* district level fixed effects coefficients */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_covW; j++) {
	    allresults[itemp++] = delta_district[i][j];
	  }      
	}
	/* village level sigmas */
	for (i = 0; i < n_grp; i++) {
	  allresults[itemp++] = sigma2_village[i][0];
	}
	/* district level sigmas */
	for (i = 0; i < n_grp; i++) {
	  allresults[itemp++] = sigma2_district[i][0];
	}
	/* acceptance ratios */
	for (i = 0; i < n_grp; i++) 
	  allresults[itemp++] = ((double) deltaCounter[i][0] / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("\n%3d percent done.\n    Metropolis acceptance ratios\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("      Control item model:");
	Rprintf(" %3g", fprec((double) deltaCounter[0][0] / (double) (main_loop + 1), 3));
	Rprintf("\n");
	for (i = 1; i < n_grp; i++) {
	  Rprintf("      Sensitive item model %3d:", i);
	  Rprintf(" %3g", fprec((double) deltaCounter[i][0] / (double) (main_loop + 1), 3));
	  Rprintf("\n");
	}
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(vitemp);
  free(treatSum);
  FreeintMatrix(Zstar, *tmax);
  free(Y0);
  FreeMatrix(X, *n_samp);
  Free3DMatrix(Xtemp, *tmax, *n_samp);
  free(Xdelta);
  free(Xdelta0);
  FreeintMatrix(gtmp, *tmax);
  FreeMatrix(V, *n_village);
  free(vitemp1);
  FreeMatrix(W, *n_district);
  FreeMatrix(delta, n_grp);
  FreeMatrix(delta_village, n_grp);
  FreeMatrix(delta_district, n_grp);
  FreeMatrix(sigma2_village, n_grp);
  FreeMatrix(sigma2_district, n_grp);
  FreeMatrix(delta0, n_grp);
  FreeMatrix(delta0_village, n_grp);
  FreeMatrix(delta0_district, n_grp);
  Free3DMatrix(A0delta, n_grp, n_dim);
  Free3DMatrix(A0delta_village, n_grp, n_dimV);
  Free3DMatrix(A0delta_district, n_grp, *n_covW);
  Free3DMatrix(deltaPro, n_grp, n_dim);
  FreeintMatrix(deltaCounter, n_grp); 

}



/** 
  Item Count Technique Binomial Mixed Effects Regression for the Multiple Sensitive Item Design
  
  This is the 4-level multilevel model 

**/

void ictregBinomMulti4Level(int *Y,             /* outcome vector */
			    int *J,             /* # of control items */
			    int *n_samp,        /* sample size */
			    int *n_draws,       /* # of MCMC draws */
			    int *treat,         /* treatment indicator vector: 0, 1, or 2 */
			    int *tmax,          /* number of sensitive items */
			    double *Xall,       /* individual-level covariates in a vector form */
			    double *Vall,       /* village-level covariates in a vector form */
			    double *Wall,       /* district-level covariates in a vector form */
			    double *Zall,       /* province-level covariates in a vector form */
			    double *beta,       /* individual-level coefs for control and sensitive items */
			    double *beta_village,  /* village-level coefs for control and sensitive items */
			    double *beta_district, /* district-level coefs for control and sensitive items */
			    double *beta_province, /* province-level coefs for control and sensitive items */
			    int *n_cov,         /* # of individual-level covariates */
			    int *n_covV,        /* # of village-level covariates */
			    int *n_covW,        /* # of district-level covariates */
			    int *n_covZ,        /* # of province-level covariates */
			    double *beta0,         /* prior mean for beta */
			    double *beta0_village, /* prior mean for beta_village */
			    double *beta0_district, /* prior mean for beta_district */
			    double *beta0_province, /* prior mean for beta_province */
			    double *A0betaAll,  /* prior precision for beta */
			    double *A0beta_villageAll,  /* prior precision for beta_village */
			    double *A0beta_districtAll,  /* prior precision for beta_district */
			    double *A0beta_provinceAll,  /* prior precision for beta_province */
			    double *betaPro,    /* proposal precision for beta */
			    int *ceiling,       /* ceiling effects */
			    int *floor,         /* floor effects */
			    int *village,       /* village indicator, 0, 1, ..., G-1 */
			    int *n_village,     /* number of villages, G */
			    int *district,      /* district indicator, 0, 1, ..., D-1 */
			    int *n_district,    /* number of districts, D */
			    int *province,      /* province indicator, 0, 1, ..., P-1 */
			    int *n_province,    /* number of districts, P */
			    double *alphaPro,   /* proposal precision for beta */
			    double *sigma2_villageAll,  /* variances for village level random effects  */
			    double *sigma2_districtAll, /* variances for district level random effects  */
			    double *sigma2_provinceAll, /* variances for province level random effects  */
			    int *nu0_village,           /* prior df for sigma2_village */
			    double *s0_village,         /* prior scale for sigma2_vilalge */
			    int *nu0_district,          /* prior df for sigma2_district */
			    double *s0_district,        /* prior scale for sigma2_district */
			    int *nu0_province,          /* prior df for sigma2_province */
			    double *s0_province,        /* prior scale for sigma2_province */
			    int *burnin,        /* number of burnins */
			    int *keep,          /* keep every *th draw */
			    int *verbose,       /* want to print progress? */
			    double *allresults  /* storage for all results */
			    ) {
  
  int i, j, k, main_loop, itemp, itempP = ftrunc((double) *n_draws/10);
  int progress = 1;
  int itempK = 1;
  double dtemp, dtemp1, dtemp2, dtemp3;

  /** get random seed **/
  GetRNGstate();

  /** Individual level Data **/
  int n_grp = *tmax + 1;
  int n_dim = *n_cov + *n_village;
  int **Zstar = intMatrix(*tmax, *n_samp);
  int *Y0 = intArray(*n_samp);
  double **X = doubleMatrix(*n_samp, n_dim); /* village dummies plus individual covariates */
  double ***Xtemp = doubleMatrix3D(*tmax, *n_samp, n_dim);
  double *Xdelta = doubleArray(*n_samp);
  double *Xdelta0 = doubleArray(*n_samp);
  int **gtmp = intMatrix(*tmax, *n_samp);
  int *vitemp = intArray(*n_village);
  int *treatSum = intArray(*tmax);
  
  itemp = 0;
  for (j = *n_village; j < n_dim; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = Xall[itemp++];
  
  for (j = 0; j < *n_village; j++) 
    vitemp[j] = 0; 
  for (i = 0; i < *tmax; i++)
    treatSum[i] = 0;

  itemp = 0;
  for (i = 0; i < *n_samp; i++) {
    Y0[i] = Y[i];
    /* village dummies */
    for (j = 0; j < *n_village; j++)
      X[i][j] = 0;
    X[i][village[i]] = 1;
    if (treat[i] > 0) {
      /* village indicator for each treatment group */
      gtmp[treat[i]-1][treatSum[treat[i]-1]] = village[i];
      /* X matrix for each treatment group */
      for (j = 0; j < n_dim; j++) 
	Xtemp[treat[i]-1][treatSum[treat[i]-1]][j] = X[i][j];
      treatSum[treat[i]-1]++;
      if ((ceiling[treat[i]-1] == 1) && (Y[i] == (*J+1))) {
	error("ceiling effects are allowed in Bayesian models only when no treated observation takes Y = J+1\n");	
      }
      if ((floor[treat[i]-1] == 1) && (Y[i] == 0)) {
	error("floor effects are allowed in Bayesian models only when no treated observation takes Y = 0\n");
      }
    }
    vitemp[village[i]]++;
  }
  /* PdoubleMatrix(X, *n_samp, n_dim); */

  /** Village level data **/
  int n_dimV = *n_covV + *n_district;
  int *vitemp1 = intArray(*n_district);
  double **V = doubleMatrix(*n_village, n_dimV); /* district dummies plus village level covariates */

  for (i = 0; i < *n_district; i++) 
    vitemp1[i] = 0;

  itemp = 0;
  for (j = *n_district; j < n_dimV; j++)
    for (i = 0; i < *n_village; i++) 
      V[i][j] = Vall[itemp++];

  for (i = 0; i < *n_village; i++) {
    for (j = 0; j < *n_district; j++)
      V[i][j] = 0;
    V[i][district[i]] = 1;
    vitemp1[district[i]]++;
  }
  /* PdoubleMatrix(V, *n_village, n_dimV); */
  
  /** District level data **/
  int n_dimW = *n_covW + *n_province;
  int *vitemp2 = intArray(*n_province);
  double **W = doubleMatrix(*n_district, n_dimW); /* province dummies plus district level covariates */

  for (i = 0; i < *n_province; i++) 
    vitemp2[i] = 0;

  itemp = 0;
  for (j = *n_province; j < n_dimW; j++)
    for (i = 0; i < *n_district; i++) 
      W[i][j] = Wall[itemp++];

  for (i = 0; i < *n_district; i++) {
    for (j = 0; j < *n_province; j++)
      W[i][j] = 0;
    W[i][province[i]] = 1;
    vitemp2[province[i]]++;
  }
  /* PdoubleMatrix(W, *n_district, n_dimW); */

  /** Province level data **/
  double **Z = doubleMatrix(*n_province, *n_covZ);

  itemp = 0;
  for (j = 0; j < *n_covZ; j++)
    for (i = 0; i < *n_province; i++) 
      Z[i][j] = Zall[itemp++];
  /* PdoubleMatrix(Z, *n_province, *n_covZ); */

  /** parameters for fixed and random effects with starting values **/
  double **delta = doubleMatrix(n_grp, n_dim);
  double **delta_village = doubleMatrix(n_grp, n_dimV);
  double **delta_district = doubleMatrix(n_grp, n_dimW);
  double **delta_province = doubleMatrix(n_grp, *n_covZ);
  double **sigma2_village = doubleMatrix(n_grp, 1);
  double **sigma2_district = doubleMatrix(n_grp, 1);
  double **sigma2_province = doubleMatrix(n_grp, 1);

  itemp = 0;
  for (i = 0; i < n_grp; i++)
    sigma2_village[i][0] = sigma2_villageAll[itemp++];
  /* PdoubleMatrix(sigma2_village, n_grp, 1); */

  itemp = 0;
  for (i = 0; i < n_grp; i++)
    sigma2_district[i][0] = sigma2_districtAll[itemp++];

  itemp = 0;
  for (i = 0; i < n_grp; i++)
    sigma2_province[i][0] = sigma2_provinceAll[itemp++];

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = 0; j < *n_covZ; j++)
      delta_province[i][j] = beta_province[itemp++];
  }
  /* PdoubleMatrix(delta_district, n_grp, *n_covW); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = *n_province; j < n_dimW; j++)
      delta_district[i][j] = beta_district[itemp++];
    for (j = 0; j < *n_province; j++) {
      delta_district[i][j] = norm_rand() * sqrt(sigma2_province[i][0]);
      for (k = 0; k < *n_covZ; k++)
	delta_district[i][j] += delta_province[i][k] * Z[j][k];
    }
  }
  /* PdoubleMatrix(delta_village, n_grp, n_dimV); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = *n_district; j < n_dimV; j++)
      delta_village[i][j] = beta_village[itemp++];
    for (j = 0; j < *n_district; j++) {
      delta_village[i][j] = norm_rand() * sqrt(sigma2_district[i][0]);
      for (k = 0; k < n_dimW; k++)
	delta_village[i][j] += delta_district[i][k] * W[j][k];
    }
  }
  /* PdoubleMatrix(delta_village, n_grp, n_dimV); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {
    for (j = *n_village; j < n_dim; j++)
      delta[i][j] = beta[itemp++];
    for (j = 0; j < *n_village; j++) {
      delta[i][j] = norm_rand() * sqrt(sigma2_village[i][0]);
      for (k = 0; k < n_dimV; k++)
	delta[i][j] += delta_village[i][k] * V[j][k];
    }
  }
  /* PdoubleMatrix(delta, n_grp, n_dim); */


  /** Prior **/
  double **delta0 = doubleMatrix(n_grp, n_dim);
  double ***A0delta = doubleMatrix3D(n_grp, n_dim, n_dim);
  double **delta0_village = doubleMatrix(n_grp, n_dimV);
  double ***A0delta_village = doubleMatrix3D(n_grp, n_dimV, n_dimV);
  double **delta0_district = doubleMatrix(n_grp, n_dimW);
  double ***A0delta_district = doubleMatrix3D(n_grp, n_dimW, n_dimW);
  double **delta0_province = doubleMatrix(n_grp, *n_covZ);
  double ***A0delta_province = doubleMatrix3D(n_grp, *n_covZ, *n_covZ);

  itemp = 0; 
  for (i = 0; i < n_grp; i++) { 
    for (j = *n_village; j < n_dim; j++)
      delta0[i][j] = beta0[itemp++];
    for (j = 0; j < *n_village; j++) {
      delta0[i][j] = 0;
      for (k = 0; k < n_dimV; k++)
	delta0[i][j] += delta_village[i][k] * V[j][k];
    }  
  }
  /* PdoubleMatrix(delta0, n_grp, n_dim); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) { 
    for (j = *n_district; j < n_dimV; j++)
      delta0_village[i][j] = beta0_village[itemp++];
    for (j = 0; j < *n_district; j++) {
      delta0_village[i][j] = 0;
      for (k = 0; k < n_dimW; k++)
	delta0_village[i][j] += delta_district[i][k] * W[j][k];
    }  
  }
  /* PdoubleMatrix(delta0, n_grp, n_dim); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) { 
    for (j = *n_province; j < n_dimW; j++)
      delta0_district[i][j] = beta0_district[itemp++];
    for (j = 0; j < *n_province; j++) {
      delta0_district[i][j] = 0;
      for (k = 0; k < *n_covZ; k++)
	delta0_district[i][j] += delta_province[i][k] * Z[j][k];
    }  
  }
  /* PdoubleMatrix(delta0, n_grp, n_dim); */

  itemp = 0; 
  for (i = 0; i < n_grp; i++) {     
    for (j = 0; j < *n_covZ; j++)
      delta0_province[i][j] = beta0_province[itemp++];
  }
  /* PdoubleMatrix(delta0_district, n_grp, *n_covW); */


  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = *n_village; j < n_dim; j++) {
      for (i = *n_village; i < n_dim; i++)
	A0delta[k][i][j] = A0betaAll[itemp++];
      for (i = 0; i < *n_village; i++)
	A0delta[k][i][j] = 0;
    }
    for (j = 0; j < *n_village; j++) {
      for (i = 0; i < n_dim; i++) {
	if (i == j) 
	  A0delta[k][i][j] = sigma2_village[k][0];
	else 
	  A0delta[k][i][j] = 0;
      } 
    }
  }
  /* PdoubleMatrix3D(A0delta, n_grp, n_dim, n_dim); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = *n_district; j < n_dimV; j++) {
      for (i = *n_district; i < n_dimV; i++)
	A0delta_village[k][i][j] = A0beta_villageAll[itemp++];
      for (i = 0; i < *n_district; i++)
	A0delta_village[k][i][j] = 0;
    }
    for (j = 0; j < *n_district; j++) {
      for (i = 0; i < n_dimV; i++) {
	if (i == j) 
	  A0delta_village[k][i][j] = sigma2_district[k][0];
	else 
	  A0delta_village[k][i][j] = 0;
      } 
    }
  }
  /* PdoubleMatrix3D(A0delta, n_grp, n_dim, n_dim); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = *n_province; j < n_dimW; j++) {
      for (i = *n_province; i < n_dimW; i++)
	A0delta_district[k][i][j] = A0beta_districtAll[itemp++];
      for (i = 0; i < *n_province; i++)
	A0delta_district[k][i][j] = 0;
    }
    for (j = 0; j < *n_province; j++) {
      for (i = 0; i < n_dimW; i++) {
	if (i == j) 
	  A0delta_district[k][i][j] = sigma2_province[k][0];
	else 
	  A0delta_district[k][i][j] = 0;
      } 
    }
  }
  /* PdoubleMatrix3D(A0delta, n_grp, n_dim, n_dim); */

  itemp = 0;
  for (k = 0; k < n_grp; k++) {
    for (j = 0; j < *n_covZ; j++) {
      for (i = 0; i < *n_covZ; i++)
	A0delta_province[k][i][j] = A0beta_provinceAll[itemp++];
    }
  }
  /* PdoubleMatrix3D(A0delta_district, n_grp, *n_covW, *n_covW); */

  /** Proposal precisoin **/
  double ***deltaPro = doubleMatrix3D(n_grp, n_dim, n_dim);

  itemp = 0;
  for (k = 0; k < n_grp; k++)
    for (j = 0; j < *n_village; j++) {
      for (i = 0; i < *n_village; i++)
	if (i == j) 
	  deltaPro[k][i][j] = alphaPro[itemp++];
	else 
	  deltaPro[k][i][j] = 0;
      for (i = *n_village; i < n_dim; i++)
	deltaPro[k][i][j] = 0;
    }

  itemp = 0;
  for (k = 0; k < n_grp; k++)
    for (j = *n_village; j < n_dim; j++) {
      for (i = *n_village; i < n_dim; i++) 
	if (i == j) 
	  deltaPro[k][i][j] = betaPro[itemp++];
	else
	  deltaPro[k][i][j] = 0;
      for (i = 0; i < *n_village; i++) 
	deltaPro[k][i][j] = 0;
    }
  /* PdoubleMatrix3D(deltaPro, n_grp, n_dim, n_dim); */

  /** Counter to calculate acceptance ratio for delta */
  int **deltaCounter = intMatrix(n_grp, 1);  

  itemp = 0; 
  for (i = 0; i < n_grp; i++) 
    deltaCounter[i][0] = 0; 

  /** MCMC **/
  if (*verbose) 
    Rprintf("\n *** Starting posterior sampling... ***\n");
  for (main_loop = 0; main_loop < *n_draws; main_loop++) {
    for (j = 0; j < *n_village; j++)
      vitemp[j] = 0;
    for (j = 0; j < *tmax; j++)
      treatSum[j] = 0;
    /* Sample Zstar */
    for (i = 0; i < *n_samp; i++) {
      if (treat[i] > 0) { 
	if (Y[i] == (*J+1)) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	  Y0[i] = *J;
	} else if (Y[i] == 0) {
	  Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	  Y0[i] = 0;
	} else {
	  Xdelta[i] = 0; Xdelta0[i] = 0;
	  for (j = 0; j < n_dim; j++) { 
	    Xdelta[i] += X[i][j]*delta[treat[i]][j];
	    Xdelta0[i] += X[i][j]*delta[0][j];
	  }
	  if ((ceiling[treat[i]-1] == 1) && (Y[i] == *J)) { /* ceiling effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J-1, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp2 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(*J, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J-1;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = *J;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = *J;
	    }
	  } else if ((floor[treat[i]-1] == 1) && (Y[i] == 1)) { /* floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(1, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp3 = exp(- log1p(exp(Xdelta[i])) + dbinom(0, *J, 1 / (1 + exp(-Xdelta0[i])), 1)); 
	    dtemp = unif_rand();
	    if (dtemp < (dtemp1 / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = 0;
	    } else if (dtemp < ((dtemp1 + dtemp2) / (dtemp1 + dtemp2 + dtemp3))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 1;
	    } else {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = 0;
	    }
	  } else { /* no ceiling and floor effects */
	    dtemp1 = exp(Xdelta[i] - log1p(exp(Xdelta[i])) + 
			 dbinom(Y[i] - 1, *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    dtemp2 = exp(- log1p(exp(Xdelta[i])) + dbinom(Y[i], *J, 1 / (1 + exp(-Xdelta0[i])), 1));
	    if (unif_rand() < (dtemp1 / (dtemp1 + dtemp2))) {
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 1;
	      Y0[i] = Y[i] - 1;
	    } else { 
	      Zstar[treat[i]-1][treatSum[treat[i]-1]] = 0;
	      Y0[i] = Y[i];
	    }
	  }
	}
	treatSum[treat[i]-1]++;
      }
      vitemp[village[i]]++;
    }

    /* Sample delta for sensitive items */
    for (k = 1; k < n_grp; k++) {
      BinomLogit(Zstar[k-1], Xtemp[k-1], delta[k], treatSum[k-1], 1, n_dim, delta0[k], A0delta[k], deltaPro[k], 1, 
		 deltaCounter[k]);
    } 

    /* Sample delta for control items */
    BinomLogit(Y0, X, delta[0], *n_samp, *J, n_dim, delta0[0], A0delta[0], deltaPro[0], 1, deltaCounter[0]);
    
    /* Update the village level model */
    for (k = 0; k < n_grp; k++) {
      bNormalReg(delta[k], V, delta_village[k], sigma2_village[k], *n_village, n_dimV, 1, 1, delta0_village[k], A0delta_village[k], 
		 1, s0_village[k], nu0_village[k], 0, 0);
    }

    /* Update the district level model */
    for (k = 0; k < n_grp; k++) {
      bNormalReg(delta_village[k], W, delta_district[k], sigma2_district[k], *n_district, n_dimW, 1, 1, delta0_district[k], 
		 A0delta_district[k], 1, s0_district[k], nu0_district[k], 0, 0);
    }

    /* Update the province level model */
    for (k = 0; k < n_grp; k++) {
      bNormalReg(delta_district[k], Z, delta_province[k], sigma2_province[k], *n_province, *n_covZ, 1, 1, delta0_province[k], 
		 A0delta_province[k], 1, s0_province[k], nu0_province[k], 0, 0);
    }

    /* Update the prior information */
    for (i = 0; i < n_grp; i++) { 
      for (j = 0; j < *n_village; j++) {
	delta0[i][j] = 0;
	for (k = 0; k < n_dimV; k++)
	  delta0[i][j] += delta_village[i][k] * V[j][k];
      } 
      for (j = 0; j < *n_district; j++) {
	delta0_village[i][j] = 0;
	for (k = 0; k < n_dimW; k++)
	  delta0_village[i][j] += delta_district[i][k] * W[j][k];
      } 
      for (j = 0; j < *n_province; j++) {
	delta0_district[i][j] = 0;
	for (k = 0; k < *n_covZ; k++)
	  delta0_district[i][j] += delta_province[i][k] * Z[j][k];
      } 
    }

    /* Store the results */
    if (main_loop >= *burnin) {
      if (itempK == *keep) {
	/* individual level fixed effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = *n_village; j < n_dim; j++) {
	    allresults[itemp++] = delta[i][j];
	  }      
	}
	/* individual level random effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_village; j++) {
	    allresults[itemp++] = delta[i][j];
	  }      
	}
	/* village level fixed effects coefficients */
	for (i = 0; i < n_grp; i++) {
	  for (j = *n_district; j < n_dimV; j++) {
	    allresults[itemp++] = delta_village[i][j];
	  }      
	}
	/* village level random effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_district; j++) {
	    allresults[itemp++] = delta_village[i][j];
	  }      
	}
	/* district level fixed effects coefficients */
	for (i = 0; i < n_grp; i++) {
	  for (j = *n_province; j < n_dimW; j++) {
	    allresults[itemp++] = delta_district[i][j];
	  }      
	}
	/* district level random effects */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_province; j++) {
	    allresults[itemp++] = delta_district[i][j];
	  }      
	}
	/* province level fixed effects coefficients */
	for (i = 0; i < n_grp; i++) {
	  for (j = 0; j < *n_covZ; j++) {
	    allresults[itemp++] = delta_province[i][j];
	  }      
	}
	/* village level sigmas */
	for (i = 0; i < n_grp; i++) {
	  allresults[itemp++] = sigma2_village[i][0];
	}
	/* district level sigmas */
	for (i = 0; i < n_grp; i++) {
	  allresults[itemp++] = sigma2_district[i][0];
	}
	/* province level sigmas */
	for (i = 0; i < n_grp; i++) {
	  allresults[itemp++] = sigma2_province[i][0];
	}
	/* acceptance ratios */
	for (i = 0; i < n_grp; i++) 
	  allresults[itemp++] = ((double) deltaCounter[i][0] / (double) (main_loop + 1));
	itempK = 1;
      } else {
	itempK++;
      }
    }

    /* printing */
    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("\n%3d percent done.\n    Metropolis acceptance ratios\n", progress*10);
	itempP += ftrunc((double) *n_draws/10); 
	progress++;
	Rprintf("      Control item model:");
	Rprintf(" %3g", fprec((double) deltaCounter[0][0] / (double) (main_loop + 1), 3));
	Rprintf("\n");
	for (i = 1; i < n_grp; i++) {
	  Rprintf("      Sensitive item model %3d:", i);
	  Rprintf(" %3g", fprec((double) deltaCounter[i][0] / (double) (main_loop + 1), 3));
	  Rprintf("\n");
	}
	R_FlushConsole(); 
      }
    }

    /* allow for user interrupt */
    R_CheckUserInterrupt();
  }

  /** write out the random seed **/
  PutRNGstate();
 
  /** freeing memory **/
  free(vitemp);
  free(treatSum);
  FreeintMatrix(Zstar, *tmax);
  free(Y0);
  FreeMatrix(X, *n_samp);
  Free3DMatrix(Xtemp, *tmax, *n_samp);
  free(Xdelta);
  free(Xdelta0);
  FreeintMatrix(gtmp, *tmax);
  FreeMatrix(V, *n_village);
  free(vitemp1);
  FreeMatrix(W, *n_district);
  free(vitemp2);
  FreeMatrix(Z, *n_province);
  FreeMatrix(delta, n_grp);
  FreeMatrix(delta_village, n_grp);
  FreeMatrix(delta_district, n_grp);
  FreeMatrix(delta_province, n_grp);
  FreeMatrix(sigma2_village, n_grp);
  FreeMatrix(sigma2_district, n_grp);
  FreeMatrix(sigma2_province, n_grp);
  FreeMatrix(delta0, n_grp);
  FreeMatrix(delta0_village, n_grp);
  FreeMatrix(delta0_district, n_grp);
  FreeMatrix(delta0_province, n_grp);
  Free3DMatrix(A0delta, n_grp, n_dim);
  Free3DMatrix(A0delta_village, n_grp, n_dimV);
  Free3DMatrix(A0delta_district, n_grp, n_dimW);
  Free3DMatrix(A0delta_province, n_grp, *n_covZ);
  Free3DMatrix(deltaPro, n_grp, n_dim);
  FreeintMatrix(deltaCounter, n_grp); 

}


/**
   A Random Walk Metroplis Sampler for Binomial Logistic Regression 
   with Independent Normal Prior
   
   proposal distribution is the univariate normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.
**/

void BinomLogit(int *Y,        /* outcome variable: 0, 1, ..., J */
		double **X,    /* (N x K) covariate matrix */
		double *beta,  /* K coefficients */
		int n_samp,    /* # of obs */
		int n_size,    /* # of size, J */
		int n_cov,     /* # of covariates, K */
		double *beta0, /* K prior mean vector */
		double **A0,   /* (K x K) prior precision */
		double **Var, /* K proposal precision */
		int n_gen,     /* # of MCMC draws */
		int *counter   /* # of acceptance */
		) {
  
  int i, j, main_loop;
  double numer, denom, Xbeta, Xprop;
  double *prop = doubleArray(n_cov);

  for (main_loop = 0; main_loop < n_gen; main_loop++) {
    /** Sample from the proposal distribution **/
    rMVN(prop, beta, Var, n_cov);
    
    /** Calculating the ratio (log scale) **/
    /* prior */
    numer = dMVN(prop, beta0, A0, n_cov, 1);
    denom = dMVN(beta, beta0, A0, n_cov, 1);   
    
    /* likelihood */
    for (i = 0; i < n_samp; i++) {
      Xbeta = 0;
      Xprop = 0;
      for (j = 0; j < n_cov; j++) {
	Xbeta += X[i][j]*beta[j];
	Xprop += X[i][j]*prop[j];
      }
      denom += dbinom(Y[i], n_size, 1 / (1 + exp(-Xbeta)), 1); 
      numer += dbinom(Y[i], n_size, 1 / (1 + exp(-Xprop)), 1); 
    }
      
    /** Rejection **/
    if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
      counter[0]++;
      for (j = 0; j < n_cov; j++) {
	beta[j] = prop[j];
      }
    }
  }
  
  free(prop);
} /* end of BinomLogit */


/**
   A Random Walk Metroplis Sampler for Binomial Logistic Mixed Effects 
   Regression with Independent Normal Prior and Normal random effects.
   
   proposal distribution for fixed effects is the normal whose mean is
   the current value and variance is given by the input. 

   proposal distribution for random effects is the multivariate normal
   whose mean is the current value and variance is given by the
   current value of covariance matrix multiplied by the input tuning
   parameter. 

**/

void BinomLogitMixed(int *Y,          /* outcome variable: 0, 1, ..., J */
		     double **X,      /* (N x K) covariate matrix for
					 fixed effects */
		     double ***Z,     /* covariates for random effects 
				         organized by groups */
		     int *grp,        /* group indicator, 0, 1, ..., G-1 */
		     double *beta,    /* K coefficients for fixed effects */
		     double **gamma,  /* (G x L) matrix of random effects */
		     double **Psi,    /* LxL precision matrix for random effecs */
		     int n_samp,      /* # of obs */
		     int J,           /* size of binomial, J */
		     int n_fixed,     /* # of fixed effects, K */
		     int n_random,    /* # of random effects, L */
		     int n_grp,       /* # of groups, G */
		     double *beta0,   /* K dimensional prior mean vector */
		     double **A0,     /* (K x K) prior precision */
		     int tau0,        /* prior df for Psi */
		     double **T0,     /* prior scale for Psi */
		     double **tune_fixed,  /* K proposal variance-covariance matrix */
		     double *tune_random, /* tuning constant for random effects of each group */
		     int n_gen,        /* # of MCMC draws */
		     int *acc_fixed,   /* # of acceptance for fixed effects */
		     int *acc_random   /* # of acceptance for random effects */
		     ) {
  
  int i, j, k, main_loop;
  int *vitemp = intArray(n_grp);
  double numer, denom;
  /* proposal values */
  double *beta1 = doubleArray(n_fixed);
  double *gamma1 = doubleArray(n_random);
  /* prior for gamma = 0 */
  double *gamma0 = doubleArray(n_random);
  /* data holders */
  double *Xbeta = doubleArray(n_samp);
  double *Xbeta1 = doubleArray(n_samp);
  double *Zgamma = doubleArray(n_samp);
  double *Zgamma1 = doubleArray(n_samp);
  /* matrix holders */
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  for (j = 0; j < n_fixed; j++)
    beta1[j] = beta[j];

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /** initializing Xbeta and Zgamma **/
  for (j = 0 ; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    Xbeta[i] = 0; Zgamma[i] = 0;
    for (j = 0; j < n_fixed; j++) { 
      Xbeta[i] += X[i][j] * beta[j];
    }
    Xbeta1[i] = Xbeta[i];
    for (j = 0; j < n_random; j++) {
      Zgamma[i] += Z[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
    }
    Zgamma1[i] = Zgamma[i];
    vitemp[grp[i]]++;
  }

  /** MCMC Sampler starts here **/
  for (main_loop = 0; main_loop < n_gen; main_loop++) {

    /** STEP 1: Update Fixed Effects **/
    rMVN(beta1, beta, tune_fixed, n_fixed);
    /* Calculating the ratio (log scale) */
    /* prior */
    numer = dMVN(beta1, beta0, A0, n_fixed, 1);
    denom = dMVN(beta, beta0, A0, n_fixed, 1);   
    /* likelihood */
    for (i = 0; i < n_samp; i++) {
      Xbeta1[i] = Xbeta[i];
      for (j = 0; j < n_fixed; j++) {
	Xbeta1[i] -= X[i][j] * (beta[j] - beta1[j]);
      }
      denom += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma[i])), 1);
      numer += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta1[i]-Zgamma[i])), 1);
    }
    if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
      acc_fixed[0]++;
      for (j = 0; j < n_fixed; j++)
	beta[j] = beta1[j];
      for (i = 0; i < n_samp; i++) 
	Xbeta[i] = Xbeta1[i];
    }
    /** STEP 2: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_grp; j++) {
      for (i = 0; i < n_random; i++)
	for (k = 0; k < n_random; k++) {
	  if (i == k) {
	    mtemp[i][i] = tune_random[j];
	  } else {
	    mtemp[i][k] = 0;
	  }
	}
      rMVN(gamma1, gamma[j], mtemp, n_random);
      /* Calculating the ratio (log scale) */
      /* prior */
      numer = dMVN(gamma1, gamma0, Psi, n_random, 1);
      denom = dMVN(gamma[j], gamma0, Psi, n_random, 1); 
      /* likelihood for group j */
      for (k = 0 ; k < n_grp; k++)
	vitemp[k] = 0;
      for (i = 0; i < n_samp; i++) {
	if (grp[i] == j) {
	  Zgamma1[i] = 0;
	  for (k = 0; k < n_random; k++)
	    Zgamma1[i] += Z[j][vitemp[j]][k]*gamma1[k];
	  denom += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma[i])), 1);
	  numer += dbinom(Y[i], J, 1 / (1 + exp(-Xbeta[i]-Zgamma1[i])), 1);
	}
	vitemp[grp[i]]++;
      }
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	acc_random[j]++;
	for (k = 0; k < n_random; k++)
	  gamma[j][k] = gamma1[k];
	for (i = 0; i < n_samp; i++) {
	  if (grp[i] == j) {
	    Zgamma[i] = Zgamma1[i];
	  }      
	}
      }
    }
    
    /** STEP 3: Update Psi **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (i = 0; i < n_grp; i++)
      for (j = 0; j < n_random; j++)
	for (k = 0; k < n_random; k++)
	  mtemp[j][k] += gamma[i][j] * gamma[i][k];
    dinv(mtemp, n_random, mtemp1);
    rWish(Psi, mtemp1, (tau0+n_grp), n_random);
  }

  /* freeing memory */
  free(beta1);
  free(gamma1);
  free(gamma0);
  free(vitemp);
  free(Xbeta);
  free(Xbeta1);
  free(Zgamma);
  free(Zgamma1);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
} /* end of mixed effects logit */



/*** 
   A Gibbs Sampler for Binary Student-t Regression by   
   Chuanhai Liu
***/ 

void RobitGibbs(int *Y,        /* binary outcome variable */
		double **X,    /* covariate matrix */
		double *beta,  /* coefficients */
		int n_samp,    /* # of obs */ 
		int n_cov,     /* # of covariates */
		int prior,     /* Should prior be included in X? */
		double *beta0, /* prior mean */
		double **A0,   /* prior precision */
		int df,        /* degrees of freedom */
		int n_gen      /* # of gibbs draws */
		) {
  
  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double *tau = doubleArray(n_samp);
  double **mtemp = doubleMatrix(n_cov, n_cov);

  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    for (i = 0; i < n_samp; i++){
      dtemp = 0;
      for (j = 0; j < n_cov; j++) 
	dtemp += X[i][j]*beta[j]; 
      if(Y[i] == 0)
	W[i] = TruncT(dtemp-1000,0,dtemp,df,1,1);
      else
	W[i] = TruncT(0,dtemp+1000,dtemp,df,1,1);
      X[i][n_cov] = W[i];
      tau[i] = rgamma(((double)df + 1)/2, 
		      ((double)df + (W[i] - dtemp)*(W[i] - dtemp))/2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0;i < n_samp; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k]*tau[i];
    dtemp = rgamma(((double) df)/2, ((double) df)/2);
    for(i = n_samp;i < n_samp+n_cov; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k]*dtemp;
   
    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);

    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k];
    rMVN(beta, mean, V, n_cov);
 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(W);
  free(tau);
  free(mean);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}


/*** 
     Bayesian Normal Regression: see Chap.14 of Gelman et al. (2004) 
       both proper and improper priors (and their combinations)
       allowed for beta and sig2. 
       Model: Y_i ~ N(X_i^T beta, sig2)
       Prior: 
         if conjugate = 1,
           Prior for beta: p(beta | sig2) = N(beta | beta0, sig2 * A0)
	   Prior for sig2: p(sig2) = inv-Chi2(sig2 | nu0, s0)
	 if conjugate = 0 (semi-conjugate prior),
	   Prior for beta: p(beta) = N(beta | beta0, A0^{-1})
	   Prior for sig2: p(sig2) = inv-Chi2(sig2 | nu0, s0)
       If conjugate = 1, sig2 is sampled from its marginal
      and beta is sampled from its conditional given sig2.
       If conjugate = 0, beta is updated using its conditional given
      sig2 and sig2 is updated using its conditional given beta.
      In this case, starting values for beta and sig2 must be provided.
***/
void bNormalReg(double *Y,   
		double **X,
		double *beta,  /* coefficients */
		double *sig2,  /* variance */
		int n_samp,    /* sample size */
		int n_cov,     /* # of covariates */
		int addprior,  /* Should prior on beta be incorporated
				  into D? */
		int pbeta,     /* Is prior proper for beta? */
		double *beta0, /* prior mean for normal */
		double **A0,   /* prior precision for normal; can be
				  set to zero to induce improper prior
				  for beta alone
			       */
		int psig2,     /* 0: improper prior for sig2
				  p(sig2|X) \propto 1/sig2
				  1: proper prior for sig2
				  p(sigma2|X) = InvChi2(nu0, s0)
			       */
		double s0,     /* prior scale for InvChi2 */
		int nu0,       /* prior d.f. for InvChi2 */
		int sig2fixed,  /* 1: sig2 fixed, 0: sig2 sampled */
		   int conjugate
		) {
  /* model parameters */

  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double **mtemp = doubleMatrix(n_cov, n_cov);
  double temp;

  /* storage parameters and loop counters */
  int i, j, k;  

  double **D = doubleMatrix(n_samp + n_cov, n_cov+1);
  for (i = 0; i < n_samp; i++) {
    D[i][n_cov] = Y[i];
    for (j = 0; j < n_cov; j++)
      D[i][j] = X[i][j];
  }

  /* slice sampler */
  double f, y;
  double w;
  w = .1;
  int m;
  m = 100;
  
  double L, R;
  int J, K;

  double x;

  /* if semi-conjugate model, divide y and X by sig2 */
  if (pbeta) {
    if (!conjugate) {
      for (i = 0; i < n_samp; i++){
	for (j = 0; j <= n_cov; j++) {
	  D[i][j] /= sqrt(sig2[0]);
	}
      }
    }
  }

/*   PdoubleMatrix(A0, n_cov, n_cov); */

  /* read the proper prior for beta as additional data points */
  if (addprior) {
    if (pbeta) {
      dcholdc(A0, n_cov, mtemp);
    } else {
      for (i = 0; i < n_cov; i++)
	for (j = 0; j < n_cov; j++)
	  mtemp[i][j] = 0;
    }
    for(i = 0; i < n_cov; i++) {
      D[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	D[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	D[n_samp+i][j] = mtemp[i][j];
      }
    }
  } 

/*   PdoubleMatrix(D, n_samp + 1, n_cov + 1); */
  
  /* SS matrix */
  for(j = 0; j <= n_cov; j++) {
    for(k = 0; k <= n_cov; k++) {
      SS[j][k]=0;
    }
  }
  for(i = 0;i < n_samp + n_cov; i++) {
    for(j = 0;j <= n_cov; j++) {
      for(k = 0; k <= n_cov; k++) {
	SS[j][k] += D[i][j]*D[i][k];
      }
    }
  }
  
  /* if semi-conjugate model, y and X are scaled back*/
  if (pbeta) {
    if (!conjugate) {
      for (i = 0; i < n_samp; i++){
	for (j = 0; j <= n_cov; j++) {
	  D[i][j] *= sqrt(sig2[0]);
	}
      }
    }
  } 
  
  /* SWEEP SS matrix */
  for(j = 0; j < n_cov; j++)
    SWP(SS, j, n_cov+1);

  if (pbeta) {
    if (!conjugate) {/* semi-conjugate prior case.
			prior for beta is proper */

      for(j = 0; j < n_cov; j++)
	mean[j] = SS[j][n_cov];

      /* draw beta from its conditional given sig2 */
      for(j = 0; j < n_cov; j++)
	for(k = 0; k < n_cov; k++)
	  V[j][k] = -SS[j][k];

      rMVN(beta, mean, V, n_cov);

      /* draw sig2 from its conditional given beta */
      /** sum of squared residuals  **/
      SS[n_cov][n_cov] = 0;
      for (i = 0; i < n_samp; i++) {
	temp = 0;
	for (j = 0; j < n_cov; j++) {
	  temp += D[i][j] * beta[j];
	}
	SS[n_cov][n_cov] += (D[i][n_cov] - temp) * (D[i][n_cov] - temp);
      }
      /* draw sig2 */
      if (!sig2fixed) {
	if (psig2) {
	  sig2[0] = ( SS[n_cov][n_cov] + nu0 * s0) / rchisq((double)n_samp+nu0);
	} else {
	  if (n_samp > 4) { /* if the inverse Chi-squared distribution is proper,
			       sample sig2 from the truncated inverse Chi-squared */
	    sig2[0] = TruncInvChisq(n_samp, SS[n_cov][n_cov] / n_samp, 10, 0);
	  } else { /* if the inverse Chi-squared distribution is improper */
	    /* slice sampler code with prior for sig2 being sig2^{-1} with support [0, 10]  */
	    f = pow(sig2[0], - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * sig2[0]));
	    y = runif(0, f);
	    /** "stepping out" procedure **/
	    L = sig2[0] - w * runif(0, 1);
	    if (L < 0) L = 0;
	    R = L + w;
	    if (R > 10) R = 10;
	    J = floor(m * runif(0, 1));
	    K = (m - 1) - J;
	    f = pow(L, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * L));
	    while (J > 0 && y < f && L > 0) {
	      L = L - w;
	      if (L < 0) {
		L = 0;
		break;
	      }
	      J = J - 1;
	      f = pow(L, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * L));
	    }
	    f = pow(R, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * R));
	    while (K > 0 && y < f && R < 10) {
	      R = R + w;
	      if (R > 10) {
		R = 10;
		break;
	      }
	      K = K - 1;
	      f = pow(R, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * R));
	    }
	    /** "shrinkage" procedure **/
	    do {
	      x = runif(L, R);
	      f = pow(x, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * x));
	      if (x < sig2[0]) {
		L = x;
	      } else {
		R = x;
	      }
	    } while (y > f);
	    sig2[0] = x;
	  }

/* 	  sig2[0] = SS[n_cov][n_cov] / rchisq((double)n_samp); */
/* 	  /\*  rejection sampling  *\/ */
/* 	  while (sig2[0] > 10) { */
/* 	    sig2[0] = SS[n_cov][n_cov] / rchisq((double)n_samp); */
/* 	  } */
	  
	}
      }
    
    }
  }

  if (!pbeta || conjugate) { /*  conjugate case or improper prior for beta */

    /* draw sig2 from its marginal dist */
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (!sig2fixed) {
      if (psig2) {  /* proper prior for sig2 */
	if (pbeta)   /* proper prior for beta */
	  sig2[0]=(SS[n_cov][n_cov]+nu0*s0)/rchisq((double)n_samp+nu0);
	else        /* improper prior for beta */
	  sig2[0]=(n_samp*SS[n_cov][n_cov]/(n_samp-n_cov)+nu0*s0)/rchisq((double)n_samp+nu0);
      } else         /* improper prior for sig2 */
	sig2[0]=SS[n_cov][n_cov]/rchisq((double)n_samp-n_cov);
    }

    /* draw beta from its conditional given sig2 */
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2[0];
    rMVN(beta, mean, V, n_cov);

  }

  /* freeing memory */
  free(mean);
  FreeMatrix(D, n_samp + n_cov);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}


void bNormalReg1(double **D,
		double *beta,  /* coefficients */
		double *sig2,  /* variance */
		int n_samp,    /* sample size */
		int n_cov,     /* # of covariates */
		int addprior,  /* Should prior on beta be incorporated
				  into D? */
		int pbeta,     /* Is prior proper for beta? */
		double *beta0, /* prior mean for normal */
		double **A0,   /* prior precision for normal; can be
				  set to zero to induce improper prior
				  for beta alone
			       */
		int psig2,     /* 0: improper prior for sig2
				  p(sig2|X) \propto 1/sig2
				  1: proper prior for sig2
				  p(sigma2|X) = InvChi2(nu0, s0)
			       */
		double s0,     /* prior scale for InvChi2 */
		int nu0,       /* prior d.f. for InvChi2 */
		int sig2fixed,  /* 1: sig2 fixed, 0: sig2 sampled */
		   int conjugate
		) {
  /* model parameters */

  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double **mtemp = doubleMatrix(n_cov, n_cov);
  double temp;

  /* storage parameters and loop counters */
  int i, j, k;  

  /* slice sampler */
  double f, y;
  double w;
  w = .1;
  int m;
  m = 100;
  
  double L, R;
  int J, K;

  double x;

  /* if semi-conjugate model, divide y and X by sig2 */
  if (pbeta) {
    if (!conjugate) {
      for (i = 0; i < n_samp; i++){
	for (j = 0; j <= n_cov; j++) {
	  D[i][j] /= sqrt(sig2[0]);
	}
      }
    }
  }

/*   PdoubleMatrix(A0, n_cov, n_cov); */

  /* read the proper prior for beta as additional data points */
  if (addprior) {
    if (pbeta) {
      dcholdc(A0, n_cov, mtemp);
    } else {
      for (i = 0; i < n_cov; i++)
	for (j = 0; j < n_cov; j++)
	  mtemp[i][j] = 0;
    }
    for(i = 0; i < n_cov; i++) {
      D[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	D[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	D[n_samp+i][j] = mtemp[i][j];
      }
    }
  } 

/*   PdoubleMatrix(D, n_samp + 1, n_cov + 1); */
  
  /* SS matrix */
  for(j = 0; j <= n_cov; j++) {
    for(k = 0; k <= n_cov; k++) {
      SS[j][k]=0;
    }
  }
  for(i = 0;i < n_samp + n_cov; i++) {
    for(j = 0;j <= n_cov; j++) {
      for(k = 0; k <= n_cov; k++) {
	SS[j][k] += D[i][j]*D[i][k];
      }
    }
  }
  
 /* if semi-conjugate model, y and X are scaled back*/
  if (pbeta) {
    if (!conjugate) {
      for (i = 0; i < n_samp; i++){
	for (j = 0; j <= n_cov; j++) {
	  D[i][j] *= sqrt(sig2[0]);
	}
      }
    }
  }
  
  
  /* SWEEP SS matrix */
  for(j = 0; j < n_cov; j++)
    SWP(SS, j, n_cov+1);

  if (pbeta) {
    if (!conjugate) {/* semi-conjugate prior case.
			prior for beta is proper */

      for(j = 0; j < n_cov; j++)
	mean[j] = SS[j][n_cov];

      /* draw beta from its conditional given sig2 */
      for(j = 0; j < n_cov; j++)
	for(k = 0; k < n_cov; k++)
	  V[j][k] = -SS[j][k];

      rMVN(beta, mean, V, n_cov);

      /* draw sig2 from its conditional given beta */
      /** sum of squared residuals  **/
      SS[n_cov][n_cov] = 0;
      for (i = 0; i < n_samp; i++) {
	temp = 0;
	for (j = 0; j < n_cov; j++) {
	  temp += D[i][j] * beta[j];
	}
	SS[n_cov][n_cov] += (D[i][n_cov] - temp) * (D[i][n_cov] - temp);
      }
      /* draw sig2 */
      if (!sig2fixed) {
	if (psig2) {
	  sig2[0] = ( SS[n_cov][n_cov] + nu0 * s0) / rchisq((double)n_samp+nu0);
	} else {
	  if (n_samp > 4) { /* if the inverse Chi-squared distribution is proper,
			       sample sig2 from the truncated inverse Chi-squared */
	    sig2[0] = TruncInvChisq(n_samp, SS[n_cov][n_cov] / n_samp, 10, 0);
	  } else { /* if the inverse Chi-squared distribution is improper */
	    /* slice sampler code with prior for sig2 being sig2^{-1} with support [0, 10]  */
	    f = pow(sig2[0], - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * sig2[0]));
	    y = runif(0, f);
	    /** "stepping out" procedure **/
	    L = sig2[0] - w * runif(0, 1);
	    if (L < 0) L = 0;
	    R = L + w;
	    if (R > 10) R = 10;
	    J = floor(m * runif(0, 1));
	    K = (m - 1) - J;
	    f = pow(L, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * L));
	    while (J > 0 && y < f && L > 0) {
	      L = L - w;
	      if (L < 0) {
		L = 0;
		break;
	      }
	      J = J - 1;
	      f = pow(L, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * L));
	    }
	    f = pow(R, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * R));
	    while (K > 0 && y < f && R < 10) {
	      R = R + w;
	      if (R > 10) {
		R = 10;
		break;
	      }
	      K = K - 1;
	      f = pow(R, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * R));
	    }
	    /** "shrinkage" procedure **/
	    do {
	      x = runif(L, R);
	      f = pow(x, - (double)(n_samp + 2) / 2) * exp(- SS[n_cov][n_cov] / (2 * x));
	      if (x < sig2[0]) {
		L = x;
	      } else {
		R = x;
	      }
	    } while (y > f);
	    sig2[0] = x;
	  }

/* 	  sig2[0] = SS[n_cov][n_cov] / rchisq((double)n_samp); */
/* 	  /\*  rejection sampling  *\/ */
/* 	  while (sig2[0] > 10) { */
/* 	    sig2[0] = SS[n_cov][n_cov] / rchisq((double)n_samp); */
/* 	  } */
	  
	}
      }
    
    }
  }

  if (!pbeta || conjugate) { /*  conjugate case or improper prior for beta */

    /* draw sig2 from its marginal dist */
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (!sig2fixed) {
      if (psig2) {  /* proper prior for sig2 */
	if (pbeta)   /* proper prior for beta */
	  sig2[0]=(SS[n_cov][n_cov]+nu0*s0)/rchisq((double)n_samp+nu0);
	else        /* improper prior for beta */
	  sig2[0]=(n_samp*SS[n_cov][n_cov]/(n_samp-n_cov)+nu0*s0)/rchisq((double)n_samp+nu0);
      } else         /* improper prior for sig2 */
	sig2[0]=SS[n_cov][n_cov]/rchisq((double)n_samp-n_cov);
    }

    /* draw beta from its conditional given sig2 */
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2[0];
    rMVN(beta, mean, V, n_cov);

  }

  /* freeing memory */
  free(mean);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}


/***
     A Gibbs sampler for ordered probit regression for ideal point estimation.
     works only if n_cov = 2, alpha is sampled from the marginal, and beta
     is sampled from the conditional truncated below at 0.
 ***/


void endorseoprobitMCMC(int *Y,        /* ordinal outcome variable: 0, 1,
					   dots, J-1 */
			 double **X,    /* covariate matrix */
			 double *beta,  /* coefficients */
			 double *tau,   /* J cut points: the first
					   cutpoint is set to 0 and the last
					   cutpoint is set to tau_{J-1}+1000 */
			 int n_samp,    /* # of obs */ 
			 int n_cov,     /* # of covariates */
			 int n_cat,     /* # of categories: J */
			 int prior,     /* Should prior be included in X? */
			 double *beta0, /* prior mean */
			 double **A0,   /* prior precision */
			 int mda,       /* use marginal data augmentation? */
			 int mh,        /* use metropolis-hasting step? */
			 double *prop,  /* J-2 proposal variances for MH step */
			 int *accept,   /* counter for acceptance */
			 int n_gen      /* # of gibbs draws */
			 ) {

  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_samp);           /* means for each obs */
  double *mbeta = doubleArray(n_cov);           /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double *Wmax = doubleArray(n_cat);  /* max of W in each categry: 0, 1,
					 ..., J-1 */
  double *Wmin = doubleArray(n_cat);  /* min of W in each category: 0, 1, 
					 ..., J-1 */
  
  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp;
  double *dvtemp = doubleArray(n_cat); dvtemp[0] = tau[0];
  double **mtemp = doubleMatrix(n_cov, n_cov);
  
  /* marginal data augmentation */
  double sig2; sig2 = 1;
  int nu0; nu0 = 1;
  double s0; s0 = 1;

/*   PdoubleMatrix(A0, n_cov, n_cov); */

  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    for (i = 0; i < n_samp; i++){
      mean[i] = 0;
      for (j = 0; j < n_cov; j++) {
	mean[i] += X[i][j]*beta[j]; 
      }
    }

    /* Sampling tau with MH step */
    if (mh) {
      for (j = 1; j < (n_cat-1); j++) {
	dvtemp[j] = TruncNorm(dvtemp[j-1], tau[j+1], tau[j], prop[j-1], 1);
      }
      dtemp = 0; dvtemp[n_cat-1] = dvtemp[n_cat-2] + 1000;
      for (j = 1; j < (n_cat-1); j++) 
	dtemp = dtemp + log(pnorm(tau[j+1]-tau[j], 0, sqrt(prop[j-1]), 1, 0) -
			    pnorm(dvtemp[j-1]-tau[j], 0, sqrt(prop[j-1]), 1, 0)) -
	  log(pnorm(dvtemp[j+1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0) -
	      pnorm(tau[j-1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0));
      for (i = 0; i < n_samp; i++) {
	if (Y[i] == (n_cat-1))  
	  dtemp = dtemp + pnorm(dvtemp[n_cat-2]-mean[i], 0, 1, 0, 1) -
	    pnorm(tau[n_cat-2]-mean[i], 0, 1, 0, 1);
	else if (Y[i] > 0) 
	  dtemp = dtemp + log(pnorm(dvtemp[Y[i]]-mean[i], 0, 1, 1, 0) -
			      pnorm(dvtemp[Y[i]-1]-mean[i], 0, 1, 1, 0)) -
	    log(pnorm(tau[Y[i]]-mean[i], 0, 1, 1, 0) -
		pnorm(tau[Y[i]-1]-mean[i], 0, 1, 1, 0));
      }
      if (unif_rand() < exp(dtemp)) {
	accept[0]++;
	for (j = 1; j < n_cat; j++) {
	  tau[j] = dvtemp[j];
	}
      }
    } 

    /* Sampling the Latent Variable */
    if (!mh) {
      Wmin[0] = tau[0]; Wmax[0] = tau[0]-10;
      for (j = 1; j < n_cat; j++) {
	Wmin[j] = tau[j];
	Wmax[j] = tau[j-1];
      }
    }

    if (mda) /* marginal data augmentation */ 
      sig2 = s0/rchisq((double)nu0);

    for (i = 0; i < n_samp; i++){
      if (Y[i] == 0) {
	W[i] = TruncNorm(mean[i]-1000,0,mean[i],1,0);
      } else if (Y[i] < 0) {
	W[i] = mean[i] + norm_rand();
      } else {
	W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],mean[i],1,0);
      }
      if (!mh) {
	Wmax[Y[i]] = fmax2(Wmax[Y[i]], W[i]);
	Wmin[Y[i]] = fmin2(Wmin[Y[i]], W[i]);
      }
      X[i][n_cov] = W[i]*sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;

    for(i = 0; i < n_samp+n_cov; i++)
      for(j = 0; j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];


/*     PdoubleMatrix(X, n_samp + n_cov, n_cov); */
/*     Rprintf("Oprobit Sweep matrix %3g \n", SS[n_cov][n_cov]); */

    
    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);
    
/*     /\* check scale parameters *\/ */
/*     Rprintf("Prior df %3g, Prior scale %3g, Sweep %3g \n", */
/* 	    nu0, s0, SS[n_cov][n_cov]); */


    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mbeta[j] = SS[j][n_cov];
    if (mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++)
	V[j][k]=-SS[j][k]*sig2;


/*     Rprintf("   Variance-covariance matrix\n"); */
/*     for (j = 0; j < n_cov; j++) { */
/*       for (k = 0; k < n_cov; k++) */
/* 	Rprintf("   %5g", V[j][k]); */
/*       Rprintf("\n"); */
/*     } */


    beta[0] = mbeta[0] + norm_rand()*sqrt(V[0][0]);
/*     Rprintf("   Sampled alpha: %5g\n", beta[0]); */
    
    beta[1] = TruncNorm(0, 10000,
			mbeta[1] + (V[0][1] / V[0][0])*(beta[0] - mbeta[0]),
			V[1][1] - (V[0][1]*V[0][1] / V[0][0]), 0);
/*     Rprintf("   Sampled beta: %5g\n", beta[1]); */

/*     rMVN(beta, mbeta, V, n_cov); */

    /* rescaling the parameters */
    if (mda) {
      for (j = 0; j < n_cov; j++) 
	beta[j] /= sqrt(sig2);
      for (i = 0; i < n_samp; i++)
	X[i][n_cov] /= sqrt(sig2);
    }

    /* sampling taus without MH-step */
    if (!mh) { 
      for (j = 1; j < n_cat-1; j++) {
	tau[j] = runif(fmax2(tau[j-1], Wmax[j]), 
		       fmin2(tau[j+1], Wmin[j+1]));
      }
      tau[n_cat-1] = tau[n_cat-2] + 1000;
    }
    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  
  /* freeing memory */
  FreeMatrix(SS, n_cov+1);
  free(mean);
  free(mbeta);
  FreeMatrix(V, n_cov);
  free(W);
  free(Wmax);
  free(Wmin);
  free(dvtemp);
  FreeMatrix(mtemp, n_cov);
}



/*** 
   A Gibbs Sampler for Binary Probit Regression With and Without
   Marginal Data Augmentation
   
   Marginal Data Augmentation: see p.318 of Imai and van Dyk (2005)
   Journal of Econometrics.
      Prior mean for beta will be set to zero. 
      Improper prior allowed (set A0 to be a matrix of zeros).
***/ 

void bprobitGibbs(int *Y,        /* binary outcome variable */
		  double **X,    /* covariate matrix */
		  double *beta,  /* coefficients */
		  int n_samp,    /* # of obs */ 
		  int n_cov,     /* # of covariates */
		  int prior,     /* Should prior be included in X? */
		  double *beta0, /* prior mean */
		  double **A0,   /* prior precision */
		  int mda,       /* Want to use marginal data augmentation? */ 
		  int n_gen      /* # of gibbs draws */
		  ) {
  
  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double **mtemp = doubleMatrix(n_cov, n_cov);

  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp;
  
  /* marginal data augmentation */
  double sig2 = 1;
  int nu0 = 1;
  double s0 = 1;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	if (!mda)
	  X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /* marginal data augmentation */
    if (mda) sig2 = s0/rchisq((double)nu0);
    
    for (i = 0; i < n_samp; i++){
      dtemp = 0;
      for (j = 0; j < n_cov; j++) 
	dtemp += X[i][j]*beta[j]; 
      if(Y[i] == 0) 
	W[i] = TruncNorm(dtemp-1000,0,dtemp,1,0);
      else 
	W[i] = TruncNorm(0,dtemp+1000,dtemp,1,0);
      X[i][n_cov] = W[i]*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0;i < n_samp; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];
    for(i = n_samp;i < n_samp+n_cov; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];

    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);

    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2;
    rMVN(beta, mean, V, n_cov);
 
    /* rescaling the parameters */
    if(mda) 
      for (j = 0; j < n_cov; j++) beta[j] /= sqrt(sig2);
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(W);
  free(mean);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}

