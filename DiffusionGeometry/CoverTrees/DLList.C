#include "DLList.H"
#include <assert.h>

//p had better point to 2*N ints
DLList::DLList(int NN,int* p) : N(NN), count(0), address(p), 
  first(-1), last(-1) { 
  for(int i=0;i<N;i++) {
    address[2*i+1]=address[2*i]=-1;
  }
}

void DLList::append(int key) {
  if(last==-1) { //list is empty
    last=first=key;
    setPrev(key,-1);
    setNext(key,-1);
  } else { //list is not empty
    setNext(last,key);
    setPrev(key,last);
    setNext(key,-1);
    last=key;
  }
  count++;
}

void DLList::prepend(int key) {
  if(first==-1) { // list is empty
    last=first=key;
    setPrev(key,-1);
    setNext(key,-1);
  } else { //list is not empty
    setNext(key,first);
    setPrev(first,key);
    setPrev(key,-1);
    first=key;
  }
  count++;
}

void DLList::insertBefore(int current,int key)  {
  int prev=getPrev(current);
  if(prev!=-1) { //current is not first
    setPrev(key,prev);
    setNext(key,current);
    setPrev(current,key);
    setNext(prev,key);
  } else {
    setPrev(key,-1);
    setNext(key,current);
    setPrev(current,key);
    first=key;
  }
  count++;
} 
void DLList::insertAfter(int current,int key) {
  int next=getNext(current);
  if(next!=-1) { //current is not last
  setPrev(key,current);
  setNext(key,next);
  setNext(current,key);
  setPrev(next,key);
  } else { //current is last
    setPrev(key,current);
    setNext(key,-1);
    setNext(current,key);
    last=key;
  }
  count++;
}
  
void DLList::remove(int key) {
  int prev=getPrev(key);
  int next=getNext(key);
  if(prev==-1) { // key is first
    if(next==-1) { // key is first and last
      setFirst(-1);
      setLast(-1);
    } else { //key is first and not last
      setFirst(next);
      setPrev(next,-1);
    }
  } else { // key is not first
    if(next==-1) { // key not first and is last
      setNext(prev,-1);
      setLast(prev);
    } else { //key is not first and not last
      setPrev(next,prev);
      setNext(prev,next);
    }
  }
  count--;
}

void DLList::clear() {
  int i=first;
  while(i!=-1) {
    int ii=getNext(i);
    remove(i);
    i=ii;
  }
  count=0;
}

void DLList::printOn(ostream& os) const {
  os << "first=" << first << " last=" << last << endl;
  for(int i=first;i!=-1;i=getNext(i)) {
    os << '\t' << i << endl;
  }
}
