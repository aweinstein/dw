#include <cmath>
#include "Vectors.H"
#include "Vector.H"

double getDist(const double* x,const double* y,int dim) {
  double distsq=0.0;
  for(int k=0;k<dim;k++) {
    double diff=x[k]-y[k];
    distsq+=diff*diff;
  }
  return sqrt(distsq);
}

double Vectors::getDist(int i, int j) const {
  double* a=X+i*dim;
  double* b=X+j*dim;
  return ::getDist(a,b,dim);
}

double Vectors::getDist(const Point* p, int j) const {
  
  const double* a=(double*)p->getData();
    double* b=X+j*dim;
    return ::getDist(a,b,dim);
} 

void Vectors::printOn(int i,ostream& os) const {
  os << "(";
  for(int j=0;j<dim-1;j++) {
    os << X[i*dim+j] << ",";
  }
  os << X[i*dim+dim-1] << ")";
}

void Vectors::printOn(ostream& os) const {
  int NN=getNumber();
  for(int i=0;i<NN;i++) {
    printOn(i,os);
    os << endl;
  }
}

