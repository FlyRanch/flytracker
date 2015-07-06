/* Subroutine to generate B-spline basis functions and their derivatives for uniform open knot vectors.
	C code for An Introduction to NURBS
	by David F. Rogers. Copyright (C) 2000 David F. Rogers,
	All rights reserved.
	
	Name: dbasis.c
	Language: C
	Subroutines called: none
	Book reference: Section 3.10, Ex. 3.18, Alg. p. 283


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

#define	C_IN	prhs[0]
#define	T_IN	prhs[1]
#define	Npts_IN	prhs[2]
#define	X_IN	prhs[3]


/* Output Arguments */

#define	N_OUT	plhs[0]
#define	D1_OUT	plhs[1]
#define	D2_OUT	plhs[2]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif


static void dbasisP(double N[],
		    double D1[],
		    double D2[],
		    double *c,
		    double T[],
		    double *npts,
		    double x[],
		    int sizeofT)

{

  int nplusc,numofpts,curridx;
  int i,j,k;
  double t;
  double b1,b2;
  double f1,f2,f3,f4;
  double s1,s2,s3,s4;
  double *temp,*temp1,*temp2;
    
  nplusc = *npts + *c;
  numofpts = *npts;
  
  /* Don't forget that C index starts from 0! */

  for (i = 0; i <= sizeofT - 1 ; i++){
    t = T[i];
  
    /* allows for number of defining polygon vertices = numofpts */
    temp = mxCalloc((nplusc-1),sizeof(double));		
    temp1 = mxCalloc((nplusc-1),sizeof(double));		
    temp2 = mxCalloc((nplusc-1),sizeof(double));		
    
    /*    zero the temporary arrays 
    for (j = 0; j <= nplusc-1; j++){
      temp[j] = 0;
      temp1[j] = 0;
      temp2[j] = 0;
      }*/
    
    /* calculate the first order basis functions n[j] */
    
    for (j = 0; j <= nplusc-1-1; j++){
      if (( t >= x[j]) && (t < x[j+1]))
	temp[j] = 1;
      else
	temp[j] = 0;
    }
    
    if (t == (double)x[numofpts-1]){		/*    pick up last point	*/
      temp[numofpts-1] = 1;
      temp[numofpts] = 0;
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
	
	/*       calculate second derivative */
	
	if (temp1[j] != 0)      /* if the lower order basis function is zero skip the calculation */
	  s1 = (2 * (temp1[j])) / (x[j+k-1] - x[j]);
	else
	  s1 = 0;
	
	if (temp1[j+1] != 0)      /* if the lower order basis function is zero skip the calculation */
	  s2 = (-2 * (temp1[j+1])) / (x[j+k] - x[j+1]);
	else 
	  s2 = 0;
	
	if (temp2[j] != 0)       /* if the lower order basis function is zero skip the calculation */
	  s3 = ((t - x[j]) * (temp2[j])) / (x[j+k-1] - x[j]);
	else
	  s3 = 0;
	
	if (temp2[j+1] != 0)    /* if the lower order basis function is zero skip the calculation */
	  s4 = ((x[j+k] - t) * (temp2[j+1])) / (x[j+k] - x[j+1]);
	else
	  s4 = 0;
	
	temp[j] = b1 + b2;
	temp1[j] = f1 + f2 + f3 + f4;
	temp2[j] = s1 + s2 + s3 + s4;
      }
    }
    
    /* put in i,j array	*/
    for (j = 0; j <= numofpts-1; j++) {
      /*mxSetPr(N[i+(j-1)*numofpts],temp[j]);
	mxSetPr(D1[i+(j-1)*numofpts],temp1[j]);
	mxSetPr(D2[i+(j-1)*numofpts],temp2[j]);*/
      curridx = i + j*(sizeofT);
      N[curridx] = temp[j];
      D1[curridx] = temp1[j];
      D2[curridx] = temp2[j];
    }
    
    mxFree(temp);
    mxFree(temp1);
    mxFree(temp2);

  }/* Loop over T */
  return;
} 

/* c = 4; T = linspace(-.5,.5,20); npts = 6; x = [0 oknot(npts,c,1,0)] */

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
  double *N,*D1,*D2; 
  double *c,*T,*npts,*x;
  int sizeofT; 
  int rowT,colT,nn; 
  
  /* Check for proper number of arguments */
    
  if (nrhs != 4) { 
    mexErrMsgTxt("Two input arguments required."); 
  } else if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments."); 
  } 
  
  
  /* Assign pointers to the input parameters */ 
  c = mxGetPr(C_IN); 
  T = mxGetPr(T_IN);
  npts = mxGetPr(Npts_IN);
  x = mxGetPr(X_IN);

  /* Basis function matrix is DIM = [length(t) npts]. t is a row or column
   vector*/
  rowT = mxGetM(T_IN);
  colT = mxGetN(T_IN);
  sizeofT = MAX(rowT,colT);

  nn = *npts;

  /* Create a matrix for the return argument */ 
  N_OUT = mxCreateDoubleMatrix(sizeofT,nn, mxREAL); 
  D1_OUT = mxCreateDoubleMatrix(sizeofT,nn, mxREAL); 
  D2_OUT = mxCreateDoubleMatrix(sizeofT,nn, mxREAL); 
  
   /* Assign pointers to the output parameters */ 
  N  = mxGetPr(N_OUT);
  D1 = mxGetPr(D1_OUT);
  D2 = mxGetPr(D2_OUT);
  
  /*N = plhs[0];
    D1= plhs[1];
    D2= plhs[2];*/
  
  
  /* Do the actual computations in a subroutine */
  dbasisP(N,D1,D2,c,T,npts,x,sizeofT); 
  return;
  
}
