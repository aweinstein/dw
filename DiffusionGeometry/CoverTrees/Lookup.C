#include "Lookup.H"

template<class T> Lookup<T>::Lookup(int N, int* pflags, 
  int* pindex_to_true_flags, T* pdata) :
  number(N), flags(pflags), index_to_true_flags(pindex_to_true_flags),
  number_true_flags(0), data(pdata),
  ncalls_to_get(0), ncalls_to_set(0) {
  for(int i=0;i<number;i++) {
    flags[i]=false;
  }
}

template<class T> void Lookup<T>::setData(int i,T t) {
  ncalls_to_set++;
  flags[i]=true;
  index_to_true_flags[number_true_flags++]=i;
  data[i]=t;
}

template<class T> void Lookup<T>::reset() {
  for(int i=0;i<number_true_flags;i++) {
    flags[index_to_true_flags[i]]=false;
  }
  number_true_flags=0;
}
    
