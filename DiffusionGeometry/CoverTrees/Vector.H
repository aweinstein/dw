#ifndef VECTOR_H
#define VECTOR_H

#include "Point.H"
#include <iostream>

using namespace std;

class Vector : public Point {
private:
  int dim;
  double* p;
public:
  Vector(int idim, double* pp) : dim(idim), p(pp) {}
  ~Vector() {}
  int getDim() const { return dim; }
  const void* getData() const { return p; }
  //void printOn(ostream& = cout);
};

#endif
