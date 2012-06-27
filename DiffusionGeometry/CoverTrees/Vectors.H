#ifndef VECTORS_H
#define VECTORS_H

#include "Points.H"
#include "Point.H"
#include "Lookup.C"
#include <iostream>

using namespace std;

class Vectors : public Points {
private:
  int dim;
  double* X;
public:
  Vectors(int idim, int iN, double* iX) : Points(iN), dim(idim), X(iX) {} 
  ~Vectors() {}
  int getDim() const { return dim; }
  void* getData() const { return X; }
  void* getData(int i) const { return X+i*dim; }
  double getDist(int i,int j) const;
  double getDist(const Point*,int) const;
  void printOn(int i,ostream& = cout) const;
  void printOn(ostream& = cout) const;
};

#endif
