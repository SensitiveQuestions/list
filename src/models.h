void BinomLogit(int *Y, double **X, double *beta, int n_samp, int n_size, int n_cov, 
		double *beta0, double **A0, double **InvVar, int n_gen, int *counter);

void BinomLogitMixed(int *Y, double **X, double ***Z, int *grp, double *beta, 
		     double **gamma, double **Psi, int n_samp, int J, int n_fixed,
		     int n_random, int n_grp, double *beta0, double **A0, int tau0,
		     double **T0, double **tune_fixed, double *tune_random, int n_gen,
		     int *acc_fixed, int *acc_random);

void RobitGibbs(int *Y, double **X, double *beta, int n_samp, int n_cov, int prior, double *beta0, double **A0,	int df,	int n_gen);

void endorseoprobitMCMC(int *Y, double **X, double *beta, 
			 double *tau, int n_samp, int n_cov, int n_cat, 
			 int prior, double *beta0, double **A0, int mda, 
			 int mh, double *prop, int *accept, int n_gen);
 
void bNormalReg(double *Y, double **X, double *beta, double *sig2, 
		int n_samp, int n_cov, int addprior, int pbeta, 
		double *beta0, double **A0, int psig2, double s0, 
		int nu0, int sig2fixed, int conjugate);

void bNormalReg1(double **D, double *beta, double *sig2, 
		int n_samp, int n_cov, int addprior, int pbeta, 
		double *beta0, double **A0, int psig2, double s0, 
		int nu0, int sig2fixed, int conjugate);

void bprobitGibbs(int *Y, double **X, double *beta,
		  int n_samp, int n_cov, int prior, double *beta0,
		  double **A0, int mda,	int n_gen);
