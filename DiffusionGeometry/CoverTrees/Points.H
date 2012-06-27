#ifndef POINTS_H
#define POINTS_H

#include "Point.H"

class Points {
private:
  int N;
public:
  Points(int NN): N(NN) {}
  ~Points() {}
  int getNumber() const { return N; }
  virtual double getDist(int,int) const = 0 ;
  virtual double getDist(const Point*,int) const = 0 ;
  virtual void printOn(int,ostream& = cout) const = 0 ;
  virtual void printOn(ostream& =cout) const = 0;
};

/*
class Points {
private:
  int N;
  int* found;
  int* marks;
  int marked;
  double* dist; 
public:
  Points(int,int*,int*,double*);
  ~Points() {}
  int getNumber() const { return N; }
  virtual double getDist(int i,int j) const = 0;
  virtual double getDist(Point,int j) const = 0;
  void resetForInsertion();
  virtual double getDistForInsertion(int,int);
};
*/

#endif
