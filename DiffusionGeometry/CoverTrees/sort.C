/* C, hand-coded quicksort */
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>

typedef double T;

void sort(int* index, T* data, int N)
{
  int i, j, k;
  T v, t;

  if(N<=1) return;

  // Partition elements
  v = data[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(data[++i] < v && i < N) { }
    while(data[--j] > v) { }
    if(i >= j) break;
    t = data[i]; data[i] = data[j]; data[j] = t;
    k = index[i]; index[i] = index[j]; index[j] = k;
  }
  t = data[i-1]; data[i-1] = data[0]; data[0] = t;
  k = index[i-1]; index[i-1] = index[0]; index[0] = k;
  sort(index,data,i-1);
  sort(index+i,data+i,N-i);
}

