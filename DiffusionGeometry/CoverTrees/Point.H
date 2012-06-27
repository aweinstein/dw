#ifndef POINT_H
#define POINT_H 

#include <iostream>
using namespace std;

class Point {
public:
  Point() {}
  ~Point() {}
  virtual const void* getData() const = 0;
};

#endif
