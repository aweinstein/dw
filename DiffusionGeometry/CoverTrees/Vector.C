#include <cmath>
#include "Vector.H"

/*double Vector::getDist(const void* y) const {
  double distsq=0.0;
  double* z=(double*)y;
  for(int i=0;i<dim;i++) {
    double diff=x[i]-z[i];
   distsq+=diff*diff;
  }
  return sqrt(distsq);
}  
*/
void Vector::printOn(ostream& os) const {
  os << "(";
  for(int i=0;i<dim-1;i++) {
    os << x[i] << ",";
  }
  os << x[dim-1] << ")";
}

