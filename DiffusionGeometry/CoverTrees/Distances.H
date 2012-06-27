#ifndef DISTANCES_H
#define DISTANCES_H

class Distances {
private:
  int N1;
  int N2;
public:
  Distances(int iN1,int iN2) : N1(iN1), N2(iN2) {}
  Distances(int iN) : N1(iN), N2(iN) {}
  int getN1() const { return N1; }
  int getN2() const { return N2; }  
  ~Distances() {}
  virtual double getDist(int,int) const = 0;
};



#endif
