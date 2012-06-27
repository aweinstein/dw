//Compile command: mex -g fwithin.C Cover.C DLList.C DisjointLists.C
//Vectors.C sort.C

#include "mex.h"
#include "Cover.H"
#include "Vectors.H"
#include "Vector.H"

const char* cover_in[]={
  "theta",
  "params", 
  "level_parents",
  "children",
  "radii",
  "level_counters",
  "level_offsets",
  "levels",
};

const char* within_in[]={
  "vectors",
  "within", 
  "level"
};

const char* fnames_out_2[]={
  "dist_ncalls_to_get",
  "dist_ncalls_to_set",
  "child_ncalls_to_get",
  "child_ncalls_to_set"
};
  


mwSize dims[2];
int ndims=2;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[]) {


  /* Check for proper number of arguments. */
  if (nrhs != 3) {
    mexErrMsgTxt("Three inputs required.");
  } 
 
  //the first argument is the array of vectors used to build the tree
  /*
    int nelements=mxGetNumberOfFields(prhs[0]);
  if(nelements!=1) {
    mexErrMsgTxt("Input should have one element.");
  }
  */
  
  //the second argument is the struct returned by covertree
  int nfields = mxGetNumberOfFields(prhs[1]);
  if(nfields!=8) {
    mexErrMsgTxt("Second input should have 8 fields.");
  }

  //the third argument is the struct whose first member is the array of
  //vectors being studied; 
  //whose second member is the distance; 
  //whose third member is the depth
  nfields = mxGetNumberOfFields(prhs[2]);
  if(nfields!=3) {
    mexErrMsgTxt("Third input should have three fields.");
  }
  
  const mxArray* tmp=0;

  //Check for proper number of return arguments; [D] or [D E] or [D E F]
  //Return argument one is a cell array D
  //D{i}=indices of A.vectors within distance of A.vectors(:,i) at right level
  //E{i} is corresponding distances
  //F is diagnostics

  bool dist_flag=false;
  bool diag_flag=false;
  if(nlhs==3) {
    dist_flag=true;
    diag_flag=true;
  } else if (nlhs=3) {
    dist_flag=true;
  } else {
    if(nlhs!=1) {
      mexErrMsgTxt("One, two or three return arguments required\n");
    }
  }

  //Extract appropriate members from first input; 
  //this is what what was passed to covertree

  tmp=prhs[0];
  mwSize ndims_in=mxGetNumberOfDimensions(tmp);
  const mwSize* dims_in=mxGetDimensions(tmp);
  
  int N=dims_in[ndims_in-1];
  int n=1;
  for(int i=0;i<ndims_in-1;i++) {
    n*=dims_in[i];
  }
  
  mexPrintf("n=%d N=%d\n",n,N);


  double* X=(double*)mxGetData(tmp);

  Vectors vectors(n,N,X);

  // Extract appropriate members from second input; 
  //this is what was returned from covertree

  tmp=mxGetField(prhs[1],0,cover_in[0]);
  double* ptheta=(double*)mxGetData(tmp);
  double theta=*ptheta;

  tmp=mxGetField(prhs[1],0,cover_in[1]);
  int* params=(int*)mxGetData(tmp); 

  tmp=mxGetField(prhs[1],0,cover_in[2]);
  int* lp=(int*)mxGetData(tmp); 

  tmp=mxGetField(prhs[1],0,cover_in[3]); 
  int* pchildren=(int*)mxGetData(tmp); 
  DisjointLists children(N,pchildren,false);


  int* pdescend_list=(int*)mxMalloc(2*N*sizeof(int));
  int* pdist_flags=(int*)mxMalloc(N*sizeof(int));
  int* pindices_to_dist_flags=(int*)mxMalloc(N*sizeof(int));
  double* pdistances=(double*)mxMalloc(N*sizeof(double));
  int* pcurrent_child_flags=(int*)mxMalloc(N*sizeof(int));
  int* pindices_to_current_child_flags=(int*)mxMalloc(N*sizeof(int));
  int* pcurrent_children=(int*)mxMalloc(N*sizeof(int));

//Get third input

  //tmp=prhs[2];
  tmp=mxGetField(prhs[2],0,within_in[0]);
  ndims_in=mxGetNumberOfDimensions(tmp);
  dims_in=mxGetDimensions(tmp);
  

  int N2=dims_in[ndims_in-1];
  int n2=1;
  for(int i=0;i<ndims_in-1;i++) {
    n2*=dims_in[i];
  }

  mexPrintf("N2=%d\n",N2);
  
  if(n2!=n) {
    mexPrintf("n2=%d must equal n=%d\n",n2,n);
  }

  double* Y=(double*)mxGetData(tmp);

  tmp=mxGetField(prhs[2],0,within_in[1]);
  double* pdwithin=(double*)mxGetData(tmp); 
  double dwithin=*pdwithin;

  tmp=mxGetField(prhs[2],0,within_in[2]);
  int* pdepth=(int*)mxGetData(tmp); 
  int depth=*pdepth;

  mexPrintf("point=%g dwithin=%g depth=%d\n",*Y,dwithin,depth);
  
  Cover cover(theta,
	      params,
	      &vectors,
	      lp,children,
	      pdescend_list,
	      pdist_flags,pindices_to_dist_flags,
	      pdistances,
	      pcurrent_child_flags,pindices_to_current_child_flags,
	      pcurrent_children);


  ndims=2;
  dims[0]=1;
  dims[1]=N2;

  mxArray* pointer0=mxCreateCellArray(ndims,dims);
  mxArray* pointer1=0;
  Cover::DescendList* pdescendlist=0;
  if(dist_flag) {
   pointer1= mxCreateCellArray(ndims,dims);
   pdescendlist=(Cover::DescendList*)&cover.getDescendList();
  }

  for(int i=0;i<N2;i++) {
    //mexPrintf("loop i=%d\n",i);
    Vector v(n,Y+i*n);
    cover.findWithin(&v,dwithin,depth);
    int count=cover.getDescendList().getCount();
    dims[1]=count;
    mxArray* fout=mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    int* arr=(int*)mxGetData(fout);
    cover.fillArrFromDescendList(arr);      
    /*
    bool test=cover.checkFindWithin(&v,dwithin,depth);    
    if(test) {
      mexPrintf("checkFindWithin passed\n");
    } else {
      mexPrintf("checkFindWithin failed\n");
    }
    */
    mxSetCell(pointer0,i,fout);
    if(dist_flag) {
      fout=mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
      double* dist=(double*)mxGetData(fout);
      for(int i=0;i<count;i++) {
	dist[i]=pdescendlist->getDist(&v,arr[i]);
      }
    }	
    mxSetCell(pointer1,i,fout);
    cover.clearDescendList();
  }

  plhs[0]=pointer0;
  if(dist_flag) {
    plhs[1]=pointer1;
  }
 
  if(diag_flag) {
    double* p=0;
    plhs[2]= mxCreateStructMatrix(1, 1, 4, fnames_out_2);
    ndims=2;
    dims[0]=1;
    dims[1]=1;

    mxArray* fout=0;
    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[2],0,fnames_out_2[0],fout);
    p[0]=cover.getDistNCallsToGet();

    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[2],0,fnames_out_2[1],fout);
    p[0]=cover.getDistNCallsToSet();

    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[2],0,fnames_out_2[2],fout);
    p[0]=cover.getChildrenNCallsToGet();

    fout = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    p=(double*)mxGetData(fout);
    mxSetField(plhs[2],0,fnames_out_2[3],fout);
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

