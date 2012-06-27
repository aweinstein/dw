#ifndef DLLIST_H
#define DLLIST_H
#include <iostream>

using namespace std;

class DLList { 
private:
  int N;
  int count;
  int* address;
  int first;
  int last;
public:
  DLList(int,int*);
  ~DLList() {}
  int getNumber() const { return N; }
  int getCount() const { return count; }
  int getFirst() const { return first; }
  int getLast() const { return last; }
  int setFirst(int key) { first=key; }
  int setLast(int key) { last=key; }
  int getPrev(int key) const { return address[2*key]; }
  int getNext(int key) const { return address[2*key+1]; }
  void setPrev(int key,int prev) { address[2*key]=prev; }
  void setNext(int key,int next) { address[2*key+1]=next; }
  void insertBefore(int current,int key);
  void insertAfter(int current,int key);
  void remove(int key);
  void append(int key); 
  void prepend(int key);
  void printOn(ostream&) const;
  void clear();
};
  
#endif
