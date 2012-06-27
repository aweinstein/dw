#include "DisjointLists.H"
#include <assert.h>
#ifdef MEM_DEBUG
#include "MemoryDebugger.H"
#endif

//p had better point to 4*N ints
DisjointLists::DisjointLists(int NN,int* p,bool val) : N(NN), address(p) { 
  if(val) {
    for(int i=0;i<4*N;i++) {
      address[i]=-1;
    }
  }
}

void DisjointLists::clear() {
  for(int i=0;i<4*N;i++) {
    address[i]=-1;
  }
}

void DisjointLists::append(int parent,int child) {
  int last=getLast(parent);
  if(last==-1) { //no children
    setLast(parent,child);
    setFirst(parent,child);
    setPrev(child,-1);
    setNext(child,-1);
  } else { //list is not empty
    setNext(last,child);
    setPrev(child,last);
    setNext(child,-1);
    setLast(parent,child);
  }
}

void DisjointLists::prepend(int parent,int child) {
  int first=getFirst(parent);
  if(first==-1) { // list is empty
    setFirst(parent,child);
    setLast(parent,child);
    setPrev(child,-1);
    setNext(child,-1);
  } else { //list is not empty
    setNext(child,first);
    setPrev(first,child);
    setPrev(child,-1); 
    setFirst(parent,child);
  }
}

void DisjointLists::insertBefore(int parent,int child,int child_to_insert)  {
  int prev=getPrev(child);
  if(prev!=-1) { //child is not first
    setPrev(child_to_insert,prev);
    setNext(child_to_insert,child);
    setPrev(child,child_to_insert);
    setNext(prev,child_to_insert);
  } else {
    setPrev(child_to_insert,-1);
    setNext(child_to_insert,child);
    setPrev(child,child_to_insert);
    setFirst(parent,child_to_insert);
  }
}
  
void DisjointLists::insertAfter(int parent,int child,int child_to_insert) {
  int next=getNext(child);
  if(next!=-1) { //child is not last
    setPrev(child_to_insert,child);
    setNext(child_to_insert,next);
    setNext(child,child_to_insert);
    setPrev(next,child_to_insert);
  } else { //child is last
    setPrev(child_to_insert,child);
    setNext(child_to_insert,-1);
    setNext(child,child_to_insert);
    setLast(parent,child_to_insert);
  }
}
  
void DisjointLists::remove(int parent,int child) {
  int prev=getPrev(child);
  int next=getNext(child);
  if(prev==-1) { // child is first
    if(next==-1) { // child is first and last
      setFirst(parent,-1);
      setLast(parent,-1);
    } else { //child is first and not last
      setFirst(parent,next);
      setPrev(next,-1);
    }
  } else { // child is not first
    if(next==-1) { // child not first and is last
      setNext(prev,-1);
      setLast(parent,prev);
    } else { //child is not first and not last
      setPrev(next,prev);
      setNext(prev,next);
    }
  }    
}

void DisjointLists::printOn(ostream& os) const {
  for(int i=0;i<N;i++) {
    int first=getFirst(i);
    if(first!=-1)
      os << "parent=" << i << endl;
    for(int j=first;j!=-1;j=getNext(j)) {
      os << "\t child=" << j << endl;
    }
  }
}

bool DisjointLists::isPresent(int parent,int query) const {
  int count=0;
  for(int child=getFirst(parent);child!=-1;child=getNext(child)) {
    count++;
    if(count>=N){ 
      cout << "parent " << parent << " has bad child list" << endl;
      return false;
    }
    if(child==query) {
      return true;
    }
  }
  return false;
}
