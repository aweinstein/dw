#include "Points.H"

/*
Points::Points(int NN, int* f, int* m, double* d) : N(NN), found(f), marks(m),
  marked(0), dist(d) {}

void Points::resetForInsertion() {
 for(int j=0;j<marked;j++) {
    found[marks[j]]=0;
  }
  marked=0;
}

double Points::getDistForInsertion(int i, int j) {
 if(found[j]) {
    return dist[j];
  }
 //dist_ctr++;
  double d=getDist(i,j);
  found[j]=1;
  dist[j]=d;
  marks[marked++]=j;
  return d;
}  
*/
