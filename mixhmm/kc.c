#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "include/nrutil.h"


/*****************************************************************************

$Revision: 1.1 $
$LastChangedDate: 2007-08-20 10:48:02 -0400 (Mon, 20 Aug 2007) $

This file contains C subroutines that collectively form the kc.pm module, which 
provides fast computation for commonly used mathematical and statistical 
functions in the form of a Perl module. All calculations are done in double 
precision. Please report bugs to kai@mail.med.upenn.edu.

*****************************************************************************/



/******************************************************************************
the following section contains definitions for macros
******************************************************************************/

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define PI 3.141592653579893


/*******************************************************************************
the following section contains helper functions to handle Perl-C interfaces
*******************************************************************************/

FILE *fh_stdout ()
{
	return stdout;
}

FILE *fh_stderr ()
{
	return stderr;
}


/***********************************************************************************
the following section contains statistical distributions, including PDF and CDF
***********************************************************************************/


double cdf_binomial (int n, int k, double p)
/* returns the probability of k-1 or less successes in n trials when the
        probability of a success on a single trial is p.  Here n > k > 0 must be non-negative
        integers, and p > 0.
*/
{
	double betai(double a, double b, double x);
	return 1-betai (k, n-k+1, p);
}

double cdf_chi2 (double df, double x)
/*cdf_chi2(n,x) returns the cumulative chi-squared distribution with n degrees of freedom for n > 0; 
chi2(df,x) = gammap(df/2,x/2).  Returns 0 unless x>0.
*/
{
	double gammp(double a, double x);
	double output;
	output = gammp(df/2.0, x/2.0);
	return output;
}

double cdf_poisson (int k, double x)
/* cumulative distribution function for the Poisson distribution
the value is the probability that the number of Poisson random events occurring will be between 0 and k -1 inclusive, when expected mean number is x
*/
{
	double gammp(double a, double x);
	int factorial (int i);
	double sum = 0;
	int i;
	if (k <= 12) {
		for (i=0; i<k; i++) {
			sum += exp(-x)*pow(x, i)/factorial(i);
		}
		return sum;
	} else {
		return 1-gammp (k, x);
	}
}

double cdf_normal (double x, double mu, double sigma)
/* Returns the cumulative probability density function for a normal distribution with mean as mu and standard deviation as sigma
cumulative normal distribution
                  x    2
         1      /   -t  / 2
     ---------- |  e         dt
     sqrt(2 pi) /
                 -inf
*/
{
	return (1+erf((x-mu)/(sigma*sqrt(2))))/2;
}

double cdf_stdnormal (double x)
//Returns the cumulative density function for a standard normal distribution (see http://en.wikipedia.org/wiki/Normal_distribution)
{
	return (1+erf(x/sqrt(2)))/2;
}


double cdf_f (double df1, double df2, double x)
/*F(n1,n2,x) returns the cumulative F distribution with n1 numerator and n2 denominator degrees
        of freedom, for n1, n2 > 0; returns 0 unless x > 0.
*/
{
	double betai(double a, double b, double x);
	double output;
	if (x <= 0) return 0;
	output = 1 - betai (df2/2, df1/2, df2/(df2+df1*x));
	return output;
}

double cdf_t (double df, double x)
/*returns the reverse cumulative (upper-tail) Student's t distribution with n > 0
        degrees of freedom.  This function returns probability T>t.
*/
{
	double betai(double a, double b, double x);
	double output;
	output = 1 - betai (0.5*df, 0.5, df/(df+x*x)) / 2;		//it is important to use 0.5 instead of 1/2, otherwise it will be treated as ZERO!
	return output;
}



double pdf_stdnormal (double x)
//Returns the probability density function for a standard normal distribution (see http://en.wikipedia.org/wiki/Normal_distribution)
{
	return exp(-x*x/2)/sqrt(2*PI);
}

double pdf_normal (double x, double mu, double sigma)
//Returns the probability density function for a normal distribution with mean as mu and standard deviation as sigma
{
	return exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) / (sigma * sqrt(2*PI));
}

double pdf_binomial (int n, int k, double p)
/* calculate the PDF for binomial distribution
*/
{
	double bico (int n, int k);
	if (k > n) return 0;
	if (p == 0) {
		return k==0 ? 1 : 0;
	} else if (p == 1) {
		return k==n ? 1: 0;
	} else {
		return exp(log (bico(n, k)) + k * log(p) + (n-k) * log(1-p));
	}
}

double pdf_beta (double a, double b, double x)
/* calculate the PDF of beta distribution (http://en.wikipedia.org/wiki/Beta_distribution)
Note the difference between "beta distribution" and "beta function"! the beta function is used in the calculation below.
a, b > 0; returns 0 unless 0 <= x <= 1.

*/
{
	double beta(double z, double w);
	if (x<0 | x>1) return 0;			//the beta distribution has support interval of [0,1]
	return pow(x, a-1) * pow(1-x, b-1) / beta (a, b);
}

double bico (int n, int k)
/*calculate the binomial coefficient (n, k)
*/
{
	double gammln(double x);
        return floor(0.5+exp(gammln(n+1)-gammln(k+1)-gammln(n-k+1)));
}


double invnormal (double x)
//inverse standard normal cumulative distribution function, or quantile function, can be expressed in terms of the inverse error function
{
	double inverff (double x);
	if (x>1 || x<0) nrerror ("Error in subroutine argument for invnormal: x (probability value) should be between 0 and 1");
	return sqrt(2) * inverff (2*x-1);
}


void reg_linear (double x[], double y[], int ndata, double *a, double *b, double *F, double *P)
/* perform linear regression of two variables with F and P value reported
*/
{
	int i;
	double t,sxoss,syoss,sx=0.0,sy=0.0, st2;
	double mss=0.0,rss=0.0,tss=0.0;			//model sum of squares, residual sum of squares, total sum of squares
	*b=0.0;
	
	for (i=0;i<ndata;i++) {
		sx += x[i];
		sy += y[i];
	}
	sxoss=sx/ndata;
	syoss=sy/ndata;
	for (i=0;i<ndata;i++) {
		t=x[i]-sxoss;
		st2 += t*t;
		*b += t*y[i];
	}
	*b /= st2;					//Solve for a, b
	*a=(sy-sx*(*b))/ndata;

	for (i=0;i<ndata;i++) {
		tss += (y[i]-syoss)*(y[i]-syoss);
		rss += (y[i]-*a-*b*x[i])*(y[i]-*a-*b*x[i]);
	}
	mss = tss-rss;
	*F = mss / (rss/(ndata-2));
	*P = 1-cdf_f(1,ndata-2,*F);
}



/********************************************************************************
the following section contains summary statistics
*********************************************************************************/
double mean(double *data, int n)
/*Given array data[0..n-1], returns its mean as ave and its variance as var.
see mean2 for running statistics (continuously updating statistics) calculation for mean
*/
{
	int j;
	double ave = 0.0;
	for (j=0;j<n;j++) ave += data[j];
	ave /= n;
	return ave;
}

double mean2(double *data, int n)
/*Given array data[0..n-1], returns its mean as ave and its variance as var.
unlike the mean subroutine, this subroutine continuously update the mean by calculating running statistics
*/
{
	double ave;
	int i;
	ave = data[0];
	for (i=1; i<n; i++) {
		ave += (data[i] - ave) / (i + 1);
	}
	return ave;
}

void avevar(double *data, long n, double *ave, double *var)
/* Given array data[0..n-1], returns its mean as ave and its variance as var.
this subroutine uses idea from Numerical Recipe.
it uses a two-pass formula to correct roundoff errors for variance calculation
see avevar2 subroutine for running statistics (continuously updating statistics) for mean and avariance calculation
*/
{
	long j;
	double s,ep;
	for (*ave=0.0,j=0;j<n;j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=0;j<n;j++) {
		s=data[j]-(*ave);
		ep += s;
		*var += s*s;
	}
	*var=(*var-ep*ep/n)/(n-1); 			//Corrected two-pass formula
}

void avevar2(double *data, int n, double *ave, double *var)
/* stably updating mean and variance
D.H.D. West, Updating mean and variance estimates: an improved method, Comm ACM 22:9, 532 (1979)  
it assumes that you get a weight and a data value (W_i and X_i) that you use to update the estimates XBAR and S2:

    SUMW = W_1
    M = X_1
    T = 0
    For i=2,3,...,n 
    {
       Q = X_i - M
       TEMP = SUM + W_i    // typo: He meant SUMW (I think so)
       R = Q*W_i/TEMP
       M = M + R
       T = T + R*SUMW*Q
       SUMW = TEMP
    }
    XBAR = M
    S2 = T*n/((n-1)*SUMW)
*/
{
	double sumw, m, t, q, temp, r;
	int i;
	
	sumw = 1;			//all weights are treated as 1
	m = data[0];
	t = 0;
	for (i=1; i<n; i++) {
		q = data[i] - m;
		temp = sumw + 1;	//temp is the sum of all weights for previous records
		r = q / temp;		//r is the contribution of new entry to the mean
		m += r;			//m is the current mean
		t += r * sumw * q;
		sumw = temp;
	}
	*ave = m;
	*var = t*n/((n-1)*sumw);
}

void moment(double data[], int n, double *ave, double *adev, double *sdev, double *var, double *skew, double *curt)
/*Given an array of data[0..n-1], this routine returns its mean ave, average deviation adev,
standard deviation sdev, variance var, skewness skew, and kurtosis curt.
Note: it seems that the definition of skewness and kurtosis differ between softwares. The values reported by Numerical Recipe differ from those reported by Microsoft Excel and by STATA (which does not agree with Excel either). Excel seems to agree with many online calculators
I decided to take the more commonly referred definition, where
skewness: b1=sum(xi-mean)^3/(N-1)*s^3 (for Normal Distribution, b1=0)
kurtosis: b2=sum(xi-mean)^4/(N-1)*s^4 (for Normal Distribution, b2=3)
excess kurtosis: b3=b2-3
*/
{
        void nrerror(char error_text[]);
        int j;
        double ep=0.0,s,p;

        if (n <= 1) nrerror("n must be at least 2 in moment");
        s=0.0;
        for (j=0;j<n;j++) s += data[j];
        *ave=s/n;
        *adev=(*var)=(*skew)=(*curt)=0.0;
        for (j=0;j<n;j++) {
                *adev += fabs(s=data[j]-(*ave));
                ep += s;
                *var += (p=s*s);
                *skew += (p *= s);
                *curt += (p *= s);
        }
        *adev /= n;
        *var=(*var-ep*ep/n)/(n-1);			//corrected formula. normally ep=0
        *sdev=sqrt(*var);
        if (*var) {
                *skew /= ((n-1)*(*var)*(*sdev));
                *curt=(*curt)/((n-1)*(*var)*(*var));
        } else
        	nrerror("No skew/kurtosis when variance = 0 (in moment)");
}

/****************************************************************************
the following section contains subroutines for performing statistcal tests
****************************************************************************/

void ttest_onesample (double data[], long n, double expected, double *t, double *p)
/* this subroutine perform one sample t-test to test whether mean is different from expected value
*/
{
	void avevar(double data[], long n, double *ave, double *var);
	double cdf_t (double df, double x);
	double ave, var;
	
	avevar (data, n, &ave, &var);
	*t = (ave - expected) / (sqrt(var) / sqrt(n));
	*p = 2 * (1 - cdf_t (n-1, fabs(*t)));
}
	

void ttest_ev(double data1[], long n1, double data2[], long n2, double *t, double *p)
/* this subroutine calculate 2-sided t-test assuming equal variance of two populations. the pooled variance estimate are used
*/
{
	void avevar(double data[], long n, double *ave, double *var);
	double cdf_t (double df, double x);
	double var1,var2,pvar,df,ave1,ave2;

	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	df=n1+n2-2;
	pvar=(((n1-1)*var1+(n2-1)*var2)/df) * (1.0/n1 + 1.0/n2);	//pooled variance
	*t=(ave1-ave2)/sqrt(pvar);
	*p = 2 * (1 - cdf_t (df, fabs(*t)));
}

void ttest_uev(double data1[], long n1, double data2[], long n2, double *t, double *p)
/* this subroutine calculates the 2-sided t-test assuming unequal variance of two populations
*/
{
	void avevar(double data[], long n, double *ave, double *var);
	double cdf_t (double df, double x);
	double var1,var2,df,ave1,ave2;
	
	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	df=(var1/n1+var2/n2)*(var1/n1+var2/n2)/((var1/n1)*(var1/n1)/(n1-1)+(var2/n2)*(var2/n2)/(n2-1));
	*t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
	*p = 2 * (1 - cdf_t (df, fabs(*t)));
}

void ftest(double data1[], long n1, double data2[], long n2, double *f, double *p)
/* this subroutine test whether two "normally distributed" populations have equal variance
PLEASE NOTICE THAT THE CALCULATED f VALUE IS ALWAYS LARGER THAN 1 ! SO THE *f VALUE MAY NOT ACTUALLY BE THE STANDARD DEVIATION OF FIRST DATA DIVIDED BY SECOND DATA
*/
{
	void avevar(double data[], long n, double *ave, double *var);
	double cdf_f (double df1, double df2, double x);
	double var1,var2,ave1,ave2,df1,df2;
	
	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	if (var1 > var2) {
		*f=var1/var2;
		df1=n1-1;
		df2=n2-1;
	} else {
		*f=var2/var1;
		df1=n2-1;
		df2=n1-1;
	}
	*p = 2.0 * cdf_f (df1, df2, *f);
	if (*p > 1.0) *p=2.0-*p;
}

void chi2test(double bins1[], double bins2[], int nbins, int knstrn, double *df, double *chi2, double *p)
/* Given the arrays bins1[0..nbins-1] and bins2[0..nbins-1], containing two sets of binned
data, perform chi-squared test.
the knstrn is dependent on whether there is a priori constraint on the equal 
number of total elements in two arrays. If yes, then knstrn=1; otherwise knstrn=0.
*/ 
{
	double cdf_chi2 (double df, double x);
	int j;

	*df=nbins-knstrn;
	*chi2=0.0;
	for (j=0;j<nbins;j++)
		if (bins1[j] == 0.0 && bins2[j] == 0.0)
			--(*df);
		else {
			*chi2 += (bins1[j]-bins2[j])*(bins1[j]-bins2[j]) / (bins1[j]+bins2[j]);
		}
	*p = 1 - cdf_chi2 (*df, *chi2);
}

void chi2test_onesample(double bins[], double ebins[], int nbins, int knstrn, double *df, double *chi2, double *p)
/* Given the array bins[1..nbins] containing the observed numbers of events, and an array
ebins[1..nbins] containing the expected numbers of events, perform chi-squared test.
*/
{
	double cdf_chi2 (double df, double x);
	void nrerror(char error_text[]);
	int j;

	*df=nbins-knstrn;
	*chi2=0.0;
	for (j=0;j<nbins;j++) {
		if (ebins[j] <= 0.0) nrerror("ERROR: invalid expected number encountered in chi2test_onesample: cell is less than or equal to zero");
		*chi2 += (bins[j]-ebins[j])*(bins[j]-ebins[j]) / ebins[j];
	}
	*p = 1 - cdf_chi2 (*df, *chi2);
}

void chi2test_trend (double bins[], double *chi2, double *p)
/* calculate Cochran-Armitage trend association test for 2x3 contingency table
*/
{
	double row1, row2, col1, col2, col3, total, df;
	double e0, e1;	/*expected value in the first and second cell in contingency table*/
	double cdf_chi2 (double df, double x);
	
	row1 = bins[0]+bins[1]+bins[2];
	row2 = bins[3]+bins[4]+bins[5];
	col1 = bins[0]+bins[3];
	col2 = bins[1]+bins[4];
	col3 = bins[2]+bins[5];
	total = row1+row2;
	
	if (! (row1*row2*(total * (col2 + 4*col1) - (col2 + 2*col1)*(col2 + 2*col1)))) {
		*chi2 = -1;
		*p = -1;
		return;
	}
	
	e0 = row1*col1/total;
	e1 = row1*col2/total;
	
	*chi2 = (bins[1] - e1) + 2 * (bins[0] - e0);
	*chi2 *= *chi2;
	*chi2 *= total * total * (total-1) / row1 / row2 / (total * (col2 + 4*col1) - (col2 + 2*col1)*(col2 + 2*col1));
	*p = 1 - cdf_chi2 (1, *chi2);
}

void chi2test_2by2table(double bins[], double *chi2, double *p)
/* calculate chi2test for a 2x2 contingency table
*/
{
	double row1, row2, col1, col2, total, df;
	double cdf_chi2 (double df, double x);
	
	row1 = bins[0]+bins[1];
	row2 = bins[2]+bins[3];
	col1 = bins[0]+bins[2];
	col2 = bins[1]+bins[3];
	total = row1+row2;
	
	if (! (row1*row2*col1*col2)) {
		*chi2 = -1;
		*p = -1;
		return;
	}

/*	This method (calculating expected value for each cell before calculating chi2 value) is too complicated so I use a shortcut formula instead
	void chi2test_onesample(double bins[], double ebins[], int nbins, int knstrn, double *df, double *chi2, double *p);
	double *ebins; 	
 	ebins = dvector (0, 3);
	ebins[0] = row1*col1/total;
	ebins[1] = row1*col2/total;
	ebins[2] = row2*col1/total;
	ebins[3] = row2*col2/total;
	
	chi2test_onesample (bins, ebins, 4, 3, &df, chi2, p);
	free_dvector (ebins, 0, 3);
*/

	*chi2 = (bins[0]*bins[3]-bins[1]*bins[2]);
	*chi2 *= *chi2;
	*chi2 *= total / row1 / row2 / col1 / col2;
	*p = 1 - cdf_chi2 (1, *chi2);
}


void chi2test_2by3table(double bins[], double *chi2, double *p)
/* calculate chi2test for a 2x3 contingency table (for example, the table used in genotypeic association test, or commonly referred to as "2df test")
this table has 2 rows and 3 columns, the input in bins[] array are a11, a12, a13, a21, a22, a23, respectively
*/
{
	void chi2test_onesample(double bins[], double ebins[], int nbins, int knstrn, double *df, double *chi2, double *p);
	double *ebins;
	double row1, row2, col1, col2, col3, total, df;
	
	row1 = bins[0]+bins[1]+bins[2];
	row2 = bins[3]+bins[4]+bins[5];
	col1 = bins[0]+bins[3];
	col2 = bins[1]+bins[4];
	col3 = bins[2]+bins[5];
	total = row1+row2;

	if (! (row1*row2*col1*col2*col3)) {
		*chi2 = -1;
		*p = -1;
		return;
	}
	
	ebins = dvector (0, 5);
	ebins[0] = row1*col1/total;
	ebins[1] = row1*col2/total;
	ebins[2] = row1*col3/total;
	ebins[3] = row2*col1/total;
	ebins[4] = row2*col2/total;
	ebins[5] = row2*col3/total;
	//printf ("bins=%f %f %f %f %f %f\n", bins[0], bins[1], bins[2], bins[3], bins[4], bins[5]);
	//printf ("ebins=%f %f %f %f %f %f\n", ebins[0], ebins[1], ebins[2], ebins[3], ebins[4], ebins[5]);
	
	chi2test_onesample (bins, ebins, 6, 4, &df, chi2, p);
	free_dvector (ebins, 0, 5);
}

void kstest(double data1[], long n1, double data2[], long n2, double *d, double *prob)
/* perform the Kolmogorov-Smirnov test to examine whether two data sets are drawn from the same distribution
the code is adapted from numerical recipe, but it does not handle ties well. I may change it in the future
*/
{
	double probks(double alam);
	void quick_sort(int elements, double *arr);
	long j1=0,j2=0;
	double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	quick_sort(n1,data1);
	quick_sort(n2,data2);
	en1=n1;
	en2=n2;
	*d=0.0;
	while (j1 < n1 && j2 < n2) {
		if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;
		if (d2 <= d1) fn2=j2++/en2;
		if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
	}
	en=sqrt(en1*en2/(en1+en2));
	*prob=probks((en+0.12+0.11/en)*(*d));
}

void kstest_onesample(double data[], long n, double (*func)(double), double *d, double *prob)
/* perform the Kolmogorov-Smirnov test to examine whether two data sets are drawn from the same distribution
*/
{
	double probks(double alam);
	void quick_sort(int elements, double *arr);
	long j;
	double dt,en,ff,fn,fo=0.0;

	quick_sort(n,data);
	en=n;
	*d=0.0;
	for (j=0;j<n;j++) {
		fn=j/en;
		ff=(*func)(data[j]);
		dt=FMAX(fabs(fo-ff),fabs(fn-ff));
		if (dt > *d) *d=dt;
		fo=fn;
	}
	en=sqrt(en);
	*prob=probks((en+0.12+0.11/en)*(*d));
}

#define EPS1 0.001
#define EPS2 1.0e-8
double probks(double alam)
{
	int j;
	double a2,fac=2.0,sum=0.0,term,termbf=0.0;

	a2 = -2.0*alam*alam;
	for (j=1;j<=1000;j++) {
		term=fac*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf=fabs(term);
	}
	return 1.0;
}
#undef EPS1
#undef EPS2

double fisher_exact_1sided (int a, int b, int c, int d)
/* calculate 1-sided Fisher's exact test
*/
{
	double lnfactorial (double x);
	int lownum=a;
	int temp, i, newa, newb, newc, newd;
	double prob = 0, firstprob, currentprob;
	
	if (a<b & a<c & a<d) {
		lownum=a;
	} else if (b<c & b<d) {
		temp=b; b=a; a=temp;
		temp=d; d=c; c=temp;
	} else if (c<d) {
		temp=c; c=a; a=temp;
		temp=d; d=b; b=temp;
	} else {
		temp=d; d=a; a=temp;
		temp=c; c=b; b=temp;
	}

	firstprob = lnfactorial(a+b)+lnfactorial(a+c)+lnfactorial(b+d)+lnfactorial(c+d)-lnfactorial(a+b+c+d);

	for (i=0; i<a; i++) {
		newa = i;
		newb = a+b-i;
		newc = a+c-i;
		newd = -a+d+i;
		currentprob = exp(firstprob - lnfactorial(newa) - lnfactorial(newb) - lnfactorial(newc) - lnfactorial(newd));
		prob += currentprob;
	}

	if (a*d < b*c) {			//the a cell is under-represented; calculate the sum of prob where the cell range from 0 to a
		prob += exp(firstprob - lnfactorial(a) - lnfactorial(b) - lnfactorial(c) - lnfactorial(d));
		if (prob > 1) prob = 1; if (prob < 0) prob = 0;
		return prob;
	} else {
		if (prob > 1) prob = 1; if (prob < 0) prob = 0;
		return 1-prob;
	}
}

double fisher_exact_2sided (int a, int b, int c, int d)
/* calculate 2-sided Fisher's exact test. I first re-arrange the contingency table and make sure that the first cell is the smallest cell
next determine whether the first cell is under-represented or over-represented and then decide the direction of iteration
*/
{
	double lnfactorial (double x);
	int lownum=a;
	int temp, i, newa, newb, newc, newd, maxa, mina;
	double prob = 0, firstprob, currentprob, tableprob;
	
	if (a<b & a<c & a<d) {
		lownum=a;
	} else if (b<c & b<d) {
		temp=b; b=a; a=temp;
		temp=d; d=c; c=temp;
	} else if (c<d) {
		temp=c; c=a; a=temp;
		temp=d; d=b; b=temp;
	} else {
		temp=d; d=a; a=temp;
		temp=c; c=b; b=temp;
	}

	firstprob = lnfactorial(a+b)+lnfactorial(a+c)+lnfactorial(b+d)+lnfactorial(c+d)-lnfactorial(a+b+c+d);
	tableprob = exp(firstprob - lnfactorial(a) - lnfactorial(b) - lnfactorial(c) - lnfactorial(d));			//prob of the original table

	maxa = a+b;
	if (a+c<maxa) maxa=a+c;
	if (a*d < b*c) {			//the a cell is under-represented; calculate the sum of prob where the cell range from a to large numbers until prob is less than tableprob
		for (i=a+1; i<=maxa; i++) {
			newa = i;
			newb = a+b-i;
			newc = a+c-i;
			newd = -a+d+i;
			currentprob = exp(firstprob - lnfactorial(newa) - lnfactorial(newb) - lnfactorial(newc) - lnfactorial(newd));
			if (currentprob < tableprob) break;
			prob += currentprob;
		}
	} else {				//the a cell is over-represented; calculate the sum of prob where the cell range from a to smaller number until prob is less than tableprob
		for (i=a-1; i>=0; i--) {
			newa = i;
			newb = a+b-i;
			newc = a+c-i;
			newd = -a+d+i;
			currentprob = exp(firstprob - lnfactorial(newa) - lnfactorial(newb) - lnfactorial(newc) - lnfactorial(newd));
			if (currentprob < tableprob) break;
			prob += currentprob;
		}
	}
	if (prob > 1) prob = 1; if (prob < 0) prob = 0;
	return 1-prob;
}


void quick_sort(int elements, double *arr) 
/* use quick sort algorithm to sort an array of double numbers
code adapted from http://alienryderflex.com/quicksort/
*/
{
#define  MAX_LEVELS  300
	double piv;
	int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	
	beg[0]=0; end[0]=elements;
	while (i>=0) {
		L=beg[i]; R=end[i]-1;
		if (L<R) {
			piv=arr[L];
			while (L<R) {
				while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
				while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L]; 
			}
			arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
			if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
				swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
				swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; 
			}
		} else {
			i--; 
		}
	}
#undef MAX_LEVELS
}


/*******************************************************************************
the following section contains mathematical functions
********************************************************************************/

double gammln(double x)
/*calculate the natural logarithm of the Gamma function (http://en.wikipedia.org/wiki/Gamma_function)

The gamma function is defined by
                +inf
                /   - t  (z - 1)
     gamma(z) = |  e    t       dt
                /
                0

If z is a positive integer, then

    Gamma(z)=(z-1)!,

this subroutine was modified from jdhedden's formula that I found at http://www.nr.com/forum/showthread.php?t=606
*/
{
    double tmp, ser;
    tmp = x + 4.5 - (x - 0.5) * log(x + 4.5);
    ser = 1.000000000190015 + (76.18009172947146 / x) - (86.50532032941677 / (x + 1.0)) + (24.01409824083091 / (x + 2.0)) - (1.231739572450155 / (x + 3.0)) + (0.1208650973866179e-2 / (x + 4.0)) - (0.5395239384953e-5 / (x + 5.0));
    return (log(2.5066282746310005 * ser) - tmp);
}

double lnfactorial (double x)
/*calculate the natural logarithm of the factorial of x, this is equivalent to gammln(x+1)
*/
{
	return gammln (x+1.0);
}

int factorial (int x)
{
	int i;
	int total=1;
	if (x==0) return 1;
	if (x>12) nrerror ("ERROR: the factorial subroutine cannot process numbers larger than 12 in a typical modern 32bit computer");		//when x=13, the results becomes wrong due to overflow
	for (i=1; i<=x; i++) {
		total *= i;
	}
	return total;
}

double gammp(double a, double x)
/* calculate incomplete gamma function (http://en.wikipedia.org/wiki/Incomplete_gamma_function)

                             x
                      1     /   - t  (a - 1)
gammainc(a, x) = ---------- |  e    t       dt
                  gamma(a)  /
                             0

this subroutine is modified from Numerical Recipe
there are two ways to calculate incomplete gamma function: by series representation (gsea) and by continued fraction representation (gcf)
when x<a+1, gsea converges faster; otherwise gcf converges faster
	
*/
{
        void gcf(double *gammcf, double a, double x, double *gln);
        void gser(double *gamser, double a, double x, double *gln);
        void nrerror(char error_text[]);
        double gamser,gammcf,gln;

        if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
        if (x < (a+1.0)) {
                gser(&gamser,a,x,&gln);
                return gamser;
        } else {
                gcf(&gammcf,a,x,&gln);
                return 1.0-gammcf;
        }
}

double gammq(double a, double x)
{
	double gammp(double a, double x);
	return 1.0 - gammp(a, x);
}



void gser(double *gamser, double a, double x, double *gln)
/* the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
*/
{
        double gammln(double x);
        void nrerror(char error_text[]);
        int n;
        double sum,del,ap;

        *gln=gammln(a);
        if (x <= 0.0) {
                if (x < 0.0) nrerror("x less than 0 in routine gser");
                *gamser=0.0;
                return;
        } else {
                ap=a;
                del=sum=1.0/a;
                for (n=1;n<=ITMAX;n++) {
                        ++ap;
                        del *= x/ap;
                        sum += del;
                        if (fabs(del) < fabs(sum)*EPS) {
                                *gamser=sum*exp(-x+a*log(x)-(*gln));
                                return;
                        }
                }
                nrerror("a too large, ITMAX too small in routine gser");
                return;
        }
}

void gcf(double *gammcf, double a, double x, double *gln)
/* calculate the incomplete gamma function Q(a, x) evaluated by its continued fraction representation as gammcf.

*/
{
        double gammln(double xx);
        void nrerror(char error_text[]);
        int i;
        double an,b,c,d,del,h;

        *gln=gammln(a);
        b=x+1.0-a;
        c=1.0/FPMIN;
        d=1.0/b;
        h=d;
        for (i=1;i<=ITMAX;i++) {
                an = -i*(i-a);
                b += 2.0;
                d=an*d+b;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=b+an/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                del=d*c;
                h *= del;
                if (fabs(del-1.0) < EPS) break;
        }
        if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
        *gammcf=exp(-x+a*log(x)-(*gln))*h;
}


double beta(double z, double w)
/* calculate the beta function (http://en.wikipedia.org/wiki/Beta_function)
The beta function is defined by
                   1
                  /   (z - 1)       (w - 1)
     beta(z, w) = |  t       (1 - t)       dt
                  /
                   0

z and w are non-negative shape parameters
*/
{
        double gammln(double x);
        return exp(gammln(z)+gammln(w)-gammln(z+w));
}

double betai(double a, double b, double x)
/* calculate the incomplete beta function 
                             x
                      1     /   (a - 1)       (b - 1)
betai(a, b, x) = ---------- |  t       (1 - t)       dt
                 beta(a, b) /
                             0
*/
{
        double betacf(double a, double b, double x);
        double gammln(double xx);
        void nrerror(char error_text[]);
        double bt;

        if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
        if (x == 0.0 || x == 1.0) bt=0.0;
        else
                bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
        if (x < (a+1.0)/(a+b+2.0))
                return bt*betacf(a,b,x)/a;
        else
                return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double betacf(double a, double b, double x)
/* continued fraction for incomplete beta function by modified Lentz’s method
*/
{
        void nrerror(char error_text[]);
        int m,m2;
        double aa,c,d,del,h,qab,qam,qap;

        qab=a+b;
        qap=a+1.0;
        qam=a-1.0;
        c=1.0;
        d=1.0-qab*x/qap;
        if (fabs(d) < FPMIN) d=FPMIN;
        d=1.0/d;
        h=d;
        for (m=1;m<=ITMAX;m++) {
                m2=2*m;
                aa=m*(b-m)*x/((qam+m2)*(a+m2));
                d=1.0+aa*d;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=1.0+aa/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                h *= d*c;
                aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
                d=1.0+aa*d;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=1.0+aa/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                del=d*c;
                h *= del;
                if (fabs(del-1.0) < EPS) break;
        }
        if (m > ITMAX) nrerror("a or b too big, or ITMAX too small in betacf");
        return h;
}



double inverff (double x)
/*inverse error function (see http://en.wikipedia.org/wiki/Error_function)
the formula is adapted from http://www.theorie.physik.uni-muenchen.de/~serge/erf-approx.pdf
*/
{
	double y;
	double a = 8887.0/63473.0;
	y = -2.0/PI/a - log (1-x*x)/2.0 + sqrt ((2.0/PI/a + log(1-x*x)/2.0)*(2.0/PI/a + log(1-x*x)/2.0) - log(1-x*x)/a);
	y = sqrt (y);
	if (x<0) y = -y;
	return y;
}

double erf(double x)
/* calculate the error function (http://en.wikipedia.org/wiki/Error_function)
The error function is defined by
                        x     2
                  2    /    -t
     erf(x) = -------- |  e    dt
              sqrt(pi) /
                        0
*/
{
        double gammp(double a, double x);
        return (x < 0.0) ? (-gammp(0.5,x*x)) : gammp(0.5,x*x);
}




/*******************************************************************************
random number generators
*******************************************************************************/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
/* random number generator that generate [0,1] uniform distribution
According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted
for the above values.
*/
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC



double random_exp(double lamda, long *idum)
/* random number generator for exponential distribution, with mean of lamda
*/
{
	float ran3(long *idum);
	float dum;

	do
		dum=ran3(idum);
	while (dum == 0.0);
	return -log(dum) / lamda;
}

float random_stdnormal (long *idum)
/* random number generator for standard normal distribution, with mean of zero and variance of 1
*/
{
	float ran3(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran3(idum)-1.0;
			v2=2.0*ran3(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}






/*********************************************************************************
the following section contains third-party subroutines with slight modifications
*********************************************************************************/


/*
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  
//
// Written by Jan Wigginton
*/
/*
// This code has been modified by Kai to conform to ISO C90 standard to prevent warning messages during compilation
*/
double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
	int obs_homc, obs_homr, rare_copies, genotypes, i, mid, curr_hets, curr_homr, curr_homc;
	double *het_probs;
	double sum, p_hwe;
        void nrerror(char error_text[]);

	if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
		fprintf (stderr, "NOTICE: the three genotype counts are %d  %d %d\n", obs_hets, obs_hom1, obs_hom2);
		nrerror ("FATAL ERROR - SNP-HWE: Current genotype configuration includes a negative count");
	}
	
	obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
	
	rare_copies = 2 * obs_homr + obs_hets;
	genotypes   = obs_hets + obs_homc + obs_homr;
	
	het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	if (het_probs == NULL) {
		nrerror ("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities" );
	}
   

	for (i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	/* start at midpoint */
	mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
	
	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1))
		mid++;

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
	                       / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];
	
		/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		curr_homr++;
		curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
		                    /((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];
		
		/* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
		curr_homr--;
		curr_homc--;
	}

	for (i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	p_hwe = 0.0;
	/*  p-value calculation for p_hwe  */
	for (i = 0; i <= rare_copies; i++) {
		if (het_probs[i] > het_probs[obs_hets])
			continue;
		p_hwe += het_probs[i];
	}
	
	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
	
	free(het_probs);
	
	return p_hwe;
}
