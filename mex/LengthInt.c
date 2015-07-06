/* Subroutine to integrate the tangent vector equations of a planar geometric curve
whos bend angle is parameterized by a bspline curve.  
	
	Name: modelInt.c
	Language: C
	Subroutines called: none
	
    b1        = first term of the basis function recursion relation
    b2        = second term of the basis function recursion relation
    c         = order of the B-spline basis function
    d1[]      = array containing the derivative of the basis functions
                d1[1]) contains the derivative of the basis function associated with B1 etc.
    d2[]      = array containing the derivative of the basis functions
                d2[1] contains the derivative of the basis function associated with B1 etc.
    f1        = first term of the first derivative of the basis function recursion relation
    f2        = second term of the first derivative of the basis function recursion relation
    f3        = third term of the first derivative of the basis function recursion relation
    f4        = fourth term of the first derivative of the basis function recursion relation
    npts      = number of defining polygon vertices
    n[]       = array containing the basis functions
                n[1]) contains the basis function associated with B1 etc.
    nplusc    = constant -- npts + c -- maximum knot value
    s1        = first term of the second derivative of the basis function recursion relation
    s2        = second term of the second derivative of the basis function recursion relation
    s3        = third term of the second derivative of the basis function recursion relation
    s4        = fourth term of the second derivative of the basis function recursion relation
    t         = parameter value
    temp[]    = temporary array
    x[]       = knot vector
*/	

#include	<stdio.h>
#include        <stdlib.h>
#include	<math.h>
#include        "mex.h"
#include        "matrix.h"

/* Input Arguments */
#define P_IN    prhs[0]
#define	T_IN	prhs[1]
#define	C_IN	prhs[2]
#define	X_IN	prhs[3]

/* Output Arguments */
#define	X_OUT	plhs[0]



#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/*----------------------------------------------------*/
void dbasis(double N[],
	     double D1[],
	     double *c,
	     double T,
	     int npts,
	     double x[])

{
  int nplusc,numofpts,curridx;
  int i,j,k;
  double t;
  double b1,b2;
  double f1,f2,f3,f4;
  double s1,s2,s3,s4;
  double *temp,*temp1,*temp2;
    
  nplusc = npts + *c;
  numofpts = npts;
  
  /* Don't forget that C index starts from 0! */

  t = T;
  
  /* allows for number of defining polygon vertices = numofpts */
  temp = mxCalloc((nplusc-1),sizeof(double));		
  temp1 = mxCalloc((nplusc-1),sizeof(double));
  
  /* calculate the first order basis functions n[j] */
  for (j = 0; j <= nplusc-1-1; j++){
    if (( t >= x[j]) && (t < x[j+1]))
      temp[j] = 1;
    else
      temp[j] = 0;
  }
  
  if (t == (double)x[nplusc-1]){		/*    pick up last point	*/
    temp[numofpts-1] = 1;
  }
  
  /* calculate the higher order basis functions */
    
  for (k = 2; k <= *c; k++){
    for (j = 0; j <= nplusc-k-1; j++){
      if (temp[j] != 0)    /* if the lower order basis function is zero skip the calculation */
	b1 = ((t-x[j]) * (temp[j])) / (x[j+k-1] - x[j]);
      else
	b1 = 0;
      
      if (temp[j+1] != 0)     /* if the lower order basis function is zero skip the calculation */
	b2 = ((x[j+k]-t) * (temp[j+1])) / (x[j+k] - x[j+1]);
      else
	b2 = 0;
      
      /*       calculate first derivative */
      
      if (temp[j] != 0)       /* if the lower order basis function is zero skip the calculation */
	f1 = temp[j]/(x[j+k-1] - x[j]);
      else
	f1 = 0;
      
      if (temp[j+1] != 0)     /* if the lower order basis function is zero skip the calculation */
	f2 = -temp[j+1]/(x[j+k] - x[j+1]);
      else
	f2 = 0;
      
      if (temp1[j] != 0)      /* if the lower order basis function is zero skip the calculation */
	f3 = ((t - x[j]) * (temp1[j])) / (x[j+k-1] - x[j]);
      else
	f3 = 0;
      
      if (temp1[j+1] != 0)    /* if the lower order basis function is zero skip the calculation */
	f4 = ((x[j+k] - t) * (temp1[j+1])) / (x[j+k] - x[j+1]);
      else
	f4 = 0;
      
      
      temp[j] = b1 + b2;
      temp1[j] = f1 + f2 + f3 + f4;
    }
  }
  
  /* put in i,j array	*/
  for (j = 0; j <= numofpts-1; j++) {
    N[j] = temp[j];
    D1[j] = temp1[j];
  }
  
  mxFree(temp);
  
  return;
} 

/*----------------------------------------------------*/
double tangfunc1(double T,
		 double p[],
		 double *c,
		 int npts,
		 double x[])
{
  void dbasis(double N[],double D1[], double *c, double T, int npts, double x[]);
  double a=0.0, b=0.0, V;
  double *N, *D1;
  int n;
  N = mxCalloc((npts-1),sizeof(double));
  D1 = mxCalloc((npts-1),sizeof(double));
  
  dbasis(N,D1,c,T,npts,x);
  
  /* Perform inner product between x-spline control points and basis function*/
  for (n=0; n <= npts-1; n++){
    a += D1[n] * p[n];
    b += D1[n] * p[n + npts];
  }
  
  mxFree(N);
  mxFree(D1);
  V = sqrt( pow(a,2) + pow(b,2));
  return V;
}


/*----------------------------------------------------*/
static double *ss = NULL;

void cleanup(void) {
  mxFree(ss);
}
/*----------------------------------------------------*/
/* #define FUNC(x,a1,a2,a3,a4) ((*func)(x,a1,a2,a3,a4))*/

void trapzd(double (*func)(double, double *, double *, int, double *), 
	      double a, 
	      double b, 
	      int n,
	      double *ss,
	      double p[],
	      double *c,
	      int npts,
	      double xx[])
	      

/*This routine computes the nth stage of refinement of an extended 
trapezoidal rule. func is input as a pointer to the function to be 
integrated between limits a and b, also input. When called with n=1,
the routine returns the crudest estimate of b a f(x)dx. Subsequent calls 
with n=2,3,... (in that sequential order) will improve the accuracy by 
adding 2n-2 additional interior points. */

{
  double x,tnm,sum,del;
  /*static double ss;*/
  int it,j;
  if (n == 1) {
    (*ss) = 0.5*(b-a)*(((*func)(a,p,c,npts,xx)) + ((*func)(b,p,c,npts,xx)) );
    return;
  } 
  else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm = it;
    del = (b-a)/tnm; /*This is the spacing of the points to be added.*/
    x = a + 0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += ((*func)(x,p,c,npts,xx));
    (*ss) = 0.5*((*ss) + (b-a)*sum/tnm); /*This replaces s by its refined value.*/
    return; /*ss*/
  }
}
/*----------------------------------------------------*/

#define EPS 1.0e-5
#define JMAX 20
double qsimp(double (*func)(double, double *, double *, int, double *), 
	     double a,
	     double b,
	     double *ss,
	     double p[],
	     double *c,
	     int npts,
	     double xx[])
/*Returns the integral of the function func from a to b. The parameters 
EPS can be set to the desired fractional accuracy and JMAX so that 2 to 
the power JMAX-1 is the maximum allowed number of steps. Integration is 
performed by Simpson's rule. */
{
  void trapzd(double (*func)(double, double *, double *, int, double *), 
		double a, 
		double b, 
		int n,
		double *ss,
		double p[],
		double *c,
		int npts,
		double xx[]);
  /*void nrerror(char error_text[]);*/
  int j;
  double s,st,ost=0.0,os=0.0;

  for (j=1;j<=JMAX;j++) {
    trapzd(func,a,b,j,ss,p,c,npts,xx);
    st = *ss;
    s=(4.0*st-ost)/3.0; 
    if (j > 5) /*Avoid spurious early convergence.*/
      if (fabs(s-os) < EPS*fabs(os) ||
	  (s == 0.0 && os == 0.0))
	return s;
    os=s;
    ost=st;
  }
  /*nrerror("Too many steps in routine qsimp");*/
  return s; /*Never get here.*/
  
}
/*----------------------------------------------------*/

void LengthInt(double X[],
	       double p[],
	       double *c,
	       double T[],
	       int npts,
	       double x[],
	       double *ss)
		   
{
  double tangfunc(double T, double p[], double *c, int npts, double x[]);
  double qsimp(double (*func)(double, double *, double *, int, double *), 
	       double a,
	       double b,
	       double *ss,
	       double p[],
	       double *c,
	       int npts,
	       double xx[]);
  double A, B;
  
  A = T[0];
  B = T[1];  
  
  X[0] = qsimp(tangfunc1,A,B,ss,p,c,npts,x);

  return;  
}

/*----------------------------------------------------*/
/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
  double *X; 
  double *c,*T,*x,*p;
  int npts; 
  int mm ; 
  static double *ss;
  
  /* Check for proper number of arguments */
    
  if (nrhs != 4) { 
    mexErrMsgTxt("Five input arguments required."); 
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments."); 
  } 
  
  
  /* Assign pointers to the input parameters */ 
  c = mxGetPr(C_IN); 
  T = mxGetPr(T_IN);
  x = mxGetPr(X_IN);
  p = mxGetPr(P_IN);

  /* p is a N x 2 matrix of control points for [x y]*/
  mm = mxGetM(P_IN);
  npts = mm;

  /* Create a matrix for the return argument */ 
  X_OUT = mxCreateDoubleMatrix(1,1, mxREAL); 
    
   /* Assign pointers to the output parameters */ 
  X  = mxGetPr(X_OUT);
  
  ss = mxGetPr(mxCreateDoubleScalar(0));
  mexMakeMemoryPersistent(ss);

  /* Do the actual computations in a subroutine */
  LengthInt(X,p,c,T,npts,x,ss); 
  mexAtExit(cleanup);
  return;
  
}
