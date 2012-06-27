#ifndef LOOKUP_H
#define LOOKUP_H

template<class T> class Lookup {
 private:
  int number;
  int* flags;
  int* index_to_true_flags;
  int number_true_flags;
  T* data;
  int ncalls_to_get;
  int ncalls_to_set;
 public:
  Lookup(int N,int*,int* p,T* q);
  ~Lookup() {}
  int getNumber() const { return number; }
  int getFlags(int i) const { return flags[i]; }
  int getNumberTrueFlags() const { return number_true_flags; }
  void setData(int,T);
  const T* getData(int i) { ncalls_to_get++; 
    if(flags[i]) { 
      return data+i; 
    } else {
      return 0;
    }
  }
  const T* readData(int i) { return data+i; }
  void changeData(int i, T t) { data[i]=t; }
  int getNCallsToGet() const { return ncalls_to_get; }
  int getNCallsToSet() const { return ncalls_to_set; }
  void reset();
};

#endif
