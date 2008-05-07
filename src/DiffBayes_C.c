#include <R.h>
#include <Rmath.h>

void choleski(double a[], double out[]) ; 
void inverse(double a[], double out[]) ;
void transpose(double a[], double out[]) ;
void product(double a[], double b[], double out[]) ;
void rinvwish(int n, double D[], double out[]) ; 

/******************************************/
/* main function, computing the stats ... */
/******************************************/

void DiffBayes_C (int *n, int *nsimul, int *standard, double patient[], double mean_c[], double A[], double result[])
{
	int i ; 
	double zx = {0}, zy = {0}, rxy = {0}, zstar = {0}; 
	double C[ 4] = {0}, D[ 4] = {0}, z[ 2] = {0}, T[ 4] = {0} ;
	double sigma[ 4] = {0}, mu[ 2] = {0}; 
	
	//constant part
	inverse(A, C) ; 
	choleski(C, D) ;
	
	for(i = 0; i < *nsimul; i++)
	{
		//random generation of sigma: the sum-of-squares and cross-product matrix for this iteration
		rinvwish(*n, D, sigma) ; 
		
		//estimation of mu: the means for this iteration
		choleski(sigma, T) ; 
		
		GetRNGstate();
		z[0] = rnorm(0, 1) ;
		z[1] = rnorm(0, 1) ;
		PutRNGstate();
		
		mu[0] = mean_c[0] + ((T[0]*z[0]) + (T[2]*z[1])) / sqrt(*n) ; 
		mu[1] = mean_c[1] + ((T[1]*z[0]) + (T[3]*z[1])) / sqrt(*n) ;
		
		//computation of the test (standardised vs. not)
		if(*standard)
		{
			zx = (patient[0] - mu[0]) / sqrt(sigma[0]) ;
			zy = (patient[1] - mu[1]) / sqrt(sigma[3]) ;
			rxy = sigma[1] / sqrt(sigma[0] * sigma[3]) ; 
			zstar = (zx - zy) / sqrt(2 - (2*rxy)) ;	
		}
		else 
		{	
			zstar = ((patient[0] - mu[0]) - (patient[1] - mu[1])) / sqrt(sigma[0] + sigma[3] - (2*sigma[1])) ; 
		}
		
		result[i] = pnorm(zstar, 0, 1, 1, 0) ; 
	}
}


/*******************************************/
/* function for random numbers generation  */
/* from an inverse Wishart distribution    */
/*******************************************/


	//random values from an inverse wishart distribution
	void rinvwish(int n, double D[], double out[])
	{
		double x = 0, y = 0, z = 0;
		double G[ 4]={0}, B[ 4]={0}, W[ 4]={0}; 
		double Bprime[ 4]={0}, Dprime[ 4]={0} ; 
		
		//generating random values
		GetRNGstate();
		x = rchisq(n) ;
		y = rnorm(0, 1); 
		z  = rchisq(n - 1) ;
		PutRNGstate(); 
		
		G[0] = sqrt(x) ; 
		G[1] = 0 ;
		G[2] = y ;
		G[3] = sqrt(z) ;
		
		transpose(D, Dprime) ; 
		product(G, Dprime, B) ; 
		transpose(B, Bprime) ; 
		product(Bprime, B, W) ;
			
		inverse(W, out) ;
	}
	
	
/******************************************************/
/* basic functions to make 2*2 matrices computations  */
/******************************************************/

	//Choleski's factorisation 
	void choleski(double a[], double out[])
	{
		out[0] = sqrt(a[0]); 
		out[1] = a[1] / out[0] ; 
		out[2] = 0; 
		out[3] = sqrt(a[3] - (out[1]*out[1]));
	}

	//matrix inversion
	void inverse(double a[], double out[])
	{
		double coef; 
		coef = 1 / ((a[0]*a[3]) - (a[2]*a[1])) ; 
	
		out[0] = coef * a[3] ; 
		out[1] = coef * (-a[1]) ; 
		out[2] = coef * (-a[2]) ; 
		out[3] = coef * a[0] ;
	}
	
	//matrix transposition
	void transpose(double a[], double out[])
	{
		out[0] = a[0] ; 
		out[1] = a[2] ; 
		out[2] = a[1] ; 
		out[3] = a[3] ;
	}
	
	//matrix multiplication
	void product(double a[], double b[], double out[])
	{
		out[0] = (a[0]*b[0]) + (a[2]*b[1]) ; 
		out[1] = (a[1]*b[0]) + (a[3]*b[1]) ; 
		out[2] = (a[0]*b[2]) + (a[2]*b[3]) ; 
		out[3] = (a[1]*b[2]) + (a[3]*b[3]) ; 
	}
	
