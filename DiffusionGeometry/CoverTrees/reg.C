#include "mex.h"
#include <cmath>

void f(double* X,int dim,int m) {

  int N=pow(m,dim);
  for(int i=0;i<dim;i++) {
    int p=pow(m,i);
    for(int j=0;j<N;j++) {
      X[i+j*dim]=double((j/p)%m);
      //mexPrintf("X[%d,%d]=%g\n",i,j,X[i+j*dim]);
    }
  }
}

mwSize dims[2];
int ndims=2;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[]) {
  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    mexErrMsgTxt("Two inputs required.");
  } 
 
  /* Get first field of input */
  const mxArray* tmp=prhs[0];
  if(!mxIsClass(tmp,"int32")) {
    mexErrMsgTxt("First field of input must be scalar int32\n");
  }
  int dim=*(int*)mxGetData(tmp);

  /* Get second field of input */
  tmp=prhs[1];
  if(!mxIsClass(tmp,"int32")) {
    mexErrMsgTxt("Second field of input must be scalar int32\n");
  }
  int m=*(int*)mxGetData(tmp);

  int N=pow(m,dim);

  mexPrintf("dim=%d m=%d N=%d\n",dim,m,N);
  
  dims[0]=dim;
  dims[1]=N;
  mxArray* fout=mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
  double* X=(double*)mxGetData(fout);

  f(X,dim,m);

  plhs[0]=fout;

}
