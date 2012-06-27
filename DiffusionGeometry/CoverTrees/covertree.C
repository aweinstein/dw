//mex -O -DMEX covertree.C Cover.C DisjointLists.C DLList.C Vectors.C sort.C

#include "mex.h"
#include "Cover.H"
#include "string.h"
#include "Vectors.H"


mwSize dims[2];
int ndims=2;

mwSize diag_dims[]={1,4};
int diag_ndims=2;

const char* fnames_in[]={
  "vectors",
  "theta",
  "maxdescend"
};

const char* fnames_out_0[]={
  "theta",
  "params", 
  "level_parents",
  "children",
  "radii",
  "level_counters",
  "level_offsets",
  "levels"
};

const char* fnames_out_1[]={
  "dist_ncalls_to_get",
  "dist_ncalls_to_set",
  "child_ncalls_to_get",
  "child_ncalls_to_set"
};
  
/*
params are:   root=cover.getRoot();
              N=points->getNumber();
	      Nc=cover.getCoverNumber();
              number_duplicates=cover.getNumDuplicates()
              minlevel=cover.getMinLevel();
              maxlevel=cover.getMaxLevel();
	      nlevels=cover.getNumLevels();
              maxdescend=cover.getMaxDescend();
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[]) {

  //Check for proper number of inputs
  if(nrhs != 1) {
    mexErrMsgTxt("One input required\n");
  }
  //Check for proper number of outputs

  bool diag_flag=false;
  if (nlhs == 2) {
    diag_flag=true;
  } else {
    if(nrhs !=1 ) {
    mexErrMsgTxt("One or two outputs required.");
    } 
  }
 
  int nfields = mxGetNumberOfFields(prhs[0]);
  if(nfields!=3) {
    mexErrMsgTxt("First input must have three fields.");
  }
 
  /* Get dimensions of first field of input */
  const mxArray* tmp=mxGetField(prhs[0],0,fnames_in[0]);
  mwSize ndims_in=mxGetNumberOfDimensions(tmp);
  const mwSize* dims_in=mxGetDimensions(tmp);
  if(!mxIsClass(tmp,"double")) {
    mexErrMsgTxt("First field of input must be double\n");
  }

  int N=dims_in[ndims_in-1];
  int n=1;
  for(int i=0;i<ndims_in-1;i++) {
    n*=dims_in[i];
  }
  
  mexPrintf("dim n=%d number N=%d\n",n,N);

  double* X=(double*)mxGetData(tmp);

  /* Get second field of input */

  tmp=mxGetField(prhs[0],0,fnames_in[1]);
  if(!mxIsClass(tmp,"double")) {
    mexErrMsgTxt("Second field of input must be double\n");
  }
  double* ptheta=(double*)mxGetData(tmp);  

  /* Get third field of input */

  tmp=mxGetField(prhs[0],0,fnames_in[2]);
  if(!mxIsClass(tmp,"int32")) {
    mexErrMsgTxt("Third field of input must be int32\n");
  }
  int* pmaxdescend=(int*)mxGetData(tmp); 
  int maxdescend=(int)*pmaxdescend;

  /* Create matrix for the return argument. */

  plhs[0] = mxCreateStructMatrix(1, 1, 8, fnames_out_0);
 
  mxArray* fout;
  
  dims[0]=1;
  dims[1]=1;
  fout =mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
  double* p=(double*)mxGetData(fout);
  mxSetField(plhs[0],0,fnames_out_0[0],fout);
  p[0]=ptheta[0];
  
  dims[0]=1;
  dims[1]=8;
  fout =mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
  int* params=(int*)mxGetData(fout);
  mxSetField(plhs[0],0,fnames_out_0[1], fout);

  dims[0]=2;
  dims[1]=N;
  fout =mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
  int* plp=(int*)mxGetData(fout);
  mxSetField(plhs[0],0,fnames_out_0[2], fout);

  dims[0]=4;
  dims[1]=N;
  fout =mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
  int* pchildren=(int*)mxGetData(fout);
  mxSetField(plhs[0],0,fnames_out_0[3], fout);

  int* pdescend_list=(int*)mxMalloc(2*N*sizeof(int));
  int* pdist_flags=(int*)mxMalloc(N*sizeof(int));
  int* pindices_to_dist_flags=(int*)mxMalloc(N*sizeof(int));
  double* pdistances=(double*)mxMalloc(N*sizeof(double));
  int* pcurrent_child_flags=(int*)mxMalloc(N*sizeof(int));
  int* pindices_to_current_child_flags=(int*)mxMalloc(N*sizeof(int));
  int* pcurrent_children=(int*)mxMalloc(N*sizeof(int));

  Vectors vectors(n,N,X);

  Cover cover(*ptheta,
	      maxdescend,
	      &vectors,plp,
	      pchildren,
	      pdescend_list,
	      pdist_flags,
	      pindices_to_dist_flags,
	      pdistances,
	      pcurrent_child_flags,
	      pindices_to_current_child_flags,
	      pcurrent_children);


  params[0]=cover.getRoot();
  params[1]=N;
  params[2]=cover.getCoverNumber();
  params[3]=cover.getNumDuplicates();
  params[4]=cover.getMinLevel();
  params[5]=cover.getMaxLevel();
  int numlevels=cover.getNumLevels();
  params[6]=numlevels;
  params[7]=cover.getMaxDescend();

  dims[0]=1;
  dims[1]=cover.getNumLevels();
  fout=mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
  double* pradii=(double*)mxGetData(fout);
  pradii[0]=cover.getRadius();
  for(int i=1;i<numlevels;i++) {
    pradii[i]=ptheta[0]*pradii[i-1];
  }
  mxSetField(plhs[0],0,fnames_out_0[4], fout);

  //mexPrintf("cover.getDistCounter=%d\n",cover.getDistCtr());

  
  dims[0]=1;
  dims[1]=numlevels;
  fout =mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
  int* plevel_counters=(int*)mxGetData(fout);
  mxSetField(plhs[0],0,fnames_out_0[5], fout);

  dims[0]=1;
  dims[1]=numlevels;
  fout =mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
  int* plevel_offsets=(int*)mxGetData(fout);
  mxSetField(plhs[0],0,fnames_out_0[6], fout);

  dims[0]=1;
  dims[1]=cover.getCoverNumber();
  fout =mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
  int* plevels=(int*)mxGetData(fout);
  mxSetField(plhs[0],0,fnames_out_0[7], fout);

  Levels(&cover,plevel_counters,plevel_offsets,plevels);
  
  if(diag_flag) {
    double* p=0;
    plhs[1]= mxCreateStructMatrix(1, 1, 4, fnames_out_1);
    ndims=2;
    dims[0]=1;
    dims[1]=1;

    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[1],0,fnames_out_1[0],fout);
    p[0]=cover.getDistNCallsToGet();

    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[1],0,fnames_out_1[1],fout);
    p[0]=cover.getDistNCallsToSet();

    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[1],0,fnames_out_1[2],fout);
    p[0]=cover.getChildrenNCallsToGet();

    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[1],0,fnames_out_1[3],fout);
    p[0]=cover.getChildrenNCallsToSet();
  }

  mxFree(pdescend_list);
  mxFree(pdist_flags);
  mxFree(pindices_to_dist_flags);
  mxFree(pdistances);
  mxFree(pcurrent_child_flags);
  mxFree(pindices_to_current_child_flags);
  mxFree(pcurrent_children);  
  
}
