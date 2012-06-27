#include "Cover.H"
#include <limits.h>
#ifdef MEM_DEBUG
#include "MemoryDebugger.H"
#endif
#include <cmath>
#include <float.h>
#include <assert.h>

void sort(int*,double*,int);

void Cover::DescendList::reset() {
  clear();
  dist_lookup.reset();
  current_child_lookup.reset();
}

bool Cover::setMinLevelAndBase(double dist) {
  if(radius<=dist) {
    while(radius<=dist) {
      radius/=base;
      minlevel--;
    }
    return true;
  } else {
    return false;
  }
}

Cover::Cover(double ibase,
	     int imaxdescend,
	     const Points* ppoints,
	     int* plp, 
	     int* pchildren, 
	     int* p, 
	     int* p1, 
	     int* p2, 
	     double* p3, 
	     int* r1,
	     int* r2,
	     int* r3) : 
  points(ppoints), num_duplicates(0), 
  lp(plp), 
  children(ppoints->getNumber(),pchildren), 
  descend_list(ibase,
	       this,
	       ppoints,
	       p,
	       p1,
	       p2,
	       p3,
	       r1,
	       r2,
	       r3),
  minlevel(0), maxlevel(0), numlevels(0), root(0), base(ibase), 
  maxdescend(imaxdescend), radius(1.0), cover_number(0) {

  makeNewEntry(0,minlevel,-1); // minlevel=0
  int N=points->getNumber();
  for(int i=1;i<N;i++) {
    insert(i);
  }

  maxlevel=minlevel;
  for(int i=0;i<N;i++) {
    int l=lp[2*i];
    if((l>maxlevel)&&(l!=INT_MAX))  {
      maxlevel=l;
    }
  }
  numlevels=maxlevel-minlevel+1;
}

//This constructor duplicates an existing Cover object in a useful way

/*
  params[0]=root
  params[1]=points->getNumber()
  params[2]=cover_number
  params[3]=num_duplicates
  params[4]=minlevel
  params[5]=maxlevel
  params[6]=numlevels;
  params[7]=maxdescend
*/

Cover::Cover(double ibase,
	     int* params,
             const Points* ppoints,
             int* plp, 
             DisjointLists& rchildren,
             int* paddress,
             int* pdist_flags,
             int* pindices_to_dist_flags,
             double* pdistances,
             int* pcurrent_child_flags,
             int* pindices_to_current_child_flags,
	     int* pcurrent_child) :
  base(ibase),
  points(ppoints),
  root(params[0]), 
  cover_number(params[2]),
  num_duplicates(params[3]), 
  minlevel(params[4]), 
  maxlevel(params[5]), 
  numlevels(params[6]), 
  maxdescend(params[7]),  
  lp(plp), 
  children(rchildren), 
  descend_list(ibase,
	       this,
               ppoints,
               paddress,
               pdist_flags,
	       pindices_to_dist_flags,
               pdistances,
               pcurrent_child_flags,
               pindices_to_current_child_flags,
               pcurrent_child),
  radius(1.0) {
}



Cover::~Cover() {}


void Cover::appendToChildren(int parent,int child) {
  children.append(parent,child);
}

void Cover::prependToChildren(int parent,int child) {
  children.prepend(parent,child);
}

void Cover::insertInChildrenAtRightLevel(int parent,int child,int new_level) {
  int first=getFirstChild(parent);
  int last=getLastChild(parent);
  int j=first;
  if(j==-1) { //no children
    children.prepend(parent,child);
  } else {
    for( ;(j!=-1) && (getLevel(j)<new_level); j=children.getNext(j)) {
    }
    if(j==-1) { //at end
      children.append(parent,child);
    } else {
      children.insertBefore(parent,j,child);
    }
  }
  return;
}

void Cover::makeNewEntry(int i,int l,int p) { //i=index l=level p=parent

  //if duplicate i=index l=INT_MAX p=matched index
  //if too deep  i=index l=INT_MAX p=-1
  int offset=2*i;
  lp[offset++]=l;
  lp[offset++]=p;
  if(l!=INT_MAX) 
    cover_number++;
}
    
void Cover::newRoot(double dist, int newroot) {
  setMinLevelAndBase(dist); 
  makeNewEntry(newroot,getMinLevel(),-1);
  setLevel(root,minlevel+1);
  setParent(root,newroot);
  prependToChildren(newroot,root);
  root=newroot;
}

void Cover::insert(int index) {
  
  descend_list.reset();

  double dist=descend_list.getDist(index,root);
  if(dist==0.0) { //p is root
    makeNewEntry(root,INT_MAX,0);
    num_duplicates++;
    return;
  }

  

  //p is not root
  if(dist>radius/(1.0-base)) { //S_Jbar is empty; make p new root
    newRoot(dist,index);
    return;
  }

 
  //dist<radius/base

  descend_list.initialize(this,index);
  while((descend_list.getFirst()!=-1) 
    && (descend_list.getLevel()-minlevel<maxdescend)) { //will enter
    int rnode=descend_list.descend(index);
    
    if(rnode!=-1) {
      makeNewEntry(rnode,INT_MAX,rnode);
      num_duplicates++;
      return;
    }
  }
  
  if(descend_list.getCount()==0) {
    int parent=descend_list.getParent();
    if(parent!=-1) {
      int new_level=descend_list.getParentLevel()+1;
      makeNewEntry(index,new_level,parent);
      insertInChildrenAtRightLevel(parent,index,new_level);
      setMaxLevel(new_level);
    } else { //parent not found
      newRoot(dist,index);
    }
  } else {
    makeNewEntry(index,INT_MAX,-1);
  }

  return;
}

Cover::DescendList::~DescendList() {
}


 Cover::DescendList::DescendList(double ibase,
				 Cover* pcover,
				 const Points* ppoints,
				 int* paddress,
				 int* pdist_flags,
				 int* pindices_to_dist_flags,
				 double* pdistances,
				 int* pcurrent_child_flags,
				 int* pindices_to_current_child_flags, 
				 int* pcurrent_child) :
				
    DLList(ppoints->getNumber(),paddress), 
    cover(pcover), points(ppoints),
    level(0), base(ibase), radius(1.0), lambda(ibase/(1.0-ibase)),
    parent(0), parent_level(0),
    dist_lookup(ppoints->getNumber(),pdist_flags,pindices_to_dist_flags,pdistances),
    current_child_lookup(ppoints->getNumber(),pcurrent_child_flags,pindices_to_current_child_flags,
			 pcurrent_child)
{ }

void Cover::DescendList::initialize(const Cover* cover,int index) {
  radius=cover->getRadius();
  parent_level=level=cover->getMinLevel();
  parent=-1;
  append(cover->getRoot());
}

int Cover::DescendList::getCurrentChild(Cover* cover,int i) {
  const int* pcurrent_child=current_child_lookup.getData(i);
  if(pcurrent_child) {
    return *pcurrent_child;
  } else {
    int current_child=cover->getFirstChild(i);
    current_child_lookup.setData(i,current_child);
    return current_child;
  }
}

void Cover::DescendList::setCurrentChild(int i,int j) {
  current_child_lookup.changeData(i,j);
  /*assert(current_child_flags[i]);
  current_child[i]=j;
  */
}

double Cover::DescendList::getDist(int i,int k) {
  const double* pdist=dist_lookup.getData(k);
  if(pdist) {
    return *pdist;
  } else {
    double dist=points->getDist(i,k);
    dist_lookup.setData(k,dist);
    return dist;
  }
}

double Cover::DescendList::getDist(const Point* p,int k) {
  const double* pdist=dist_lookup.getData(k);
  if(pdist) {
    return *pdist;
  } else {
    double dist=points->getDist(p,k);
    dist_lookup.setData(k,dist);
    return dist;
  }
}

int Cover::DescendList::descend(int i) {
  //constructs S_{j+1} from S_j where j=level

  double lambdaradius=lambda*radius;
  bool flag=true; //true=not found
  // append to S_{j+1} from P^{-1}[S_j] and update min_dist
  int j=getFirst();
  int k=-1;
  DisjointLists* children=cover->getChildren();
  while(j!=-1) {
    k=getCurrentChild(cover,j);
    int jj=getNext(j);
    int count=0;
    while((k!=-1) && (cover->getLevel(k)==level+1)) {
      //child_ctr++;
      int kk=children->getNext(k);      
      double dist=getDist(i,k);
      if(dist==0.0)
	return k;
      if(dist< lambdaradius) {
	prepend(k);
						     
      }
      k=kk;
    }
    setCurrentChild(j,k);
    double dist=getDist(i,j);
    if(dist==0.0)
      return j;
    if (flag && (dist<radius)){
      parent=j;
      parent_level=level;
      flag=false;
    }
    if(dist>=lambdaradius) { //prune S_j
      remove(j);
      
    }
    j=jj;
  }

  level++;
  radius*=base;

  return -1;

}  

void Cover::fillArrFromDescendList(int* arr) const {
  //arr must point to descend_list.getCount() ints
  int i=0;
  for(int index=descend_list.getFirst();index!=-1;index=descend_list.getNext(index)) {
    arr[i++]=index;
  }
  return;
}

/***************************************************************************/

//findWithin

void Cover::findWithin(const Point* p,double r,int mlevel) {
  //cout << "Cover::findWithin" << endl;
  descend_list.reset();
  descend_list.initialize(this,root);
  int m= mlevel>maxlevel ? maxlevel : mlevel;
  //nearest_ctr=0;
  for(int l=minlevel;l<m;l++) {  //will look for children of level m-1
    descend_list.descendForWithin(p,r);
  }

  /*
    cout << "descend_list.count before pruning=" 
       << descend_list.getCount() << endl;
  */
  for(int i=descend_list.getFirst();i!=-1;i=descend_list.getNext(i)) {
    double dist=descend_list.getDist(p,i);
    if(dist>r) {
      descend_list.remove(i);
     }
  }
  
  /*
    cout << "descend_list.count after pruning=" 
       << descend_list.getCount() << endl;
  */

  /*
    int count=0;
  for(int j=descend_list.getFirst();j!=-1;j=descend_list.getNext(j)) {
    cout << "\tcount=" << count++ << " j=" << j << " flag=" 
	 << descend_list.dist_lookup.getFlags(j) 
	 << " point=" << double(j)/points->getNumber() 
	 << " dist=" << *descend_list.dist_lookup.readData(j) 
	 << endl;
  }
  */
}

void Cover::DescendList::descendForWithin(const Point* p,double r) {

  double lambdaradius=lambda*radius;
  // append to S_{j+1} from P^{-1}[S_j]
  int j=getFirst();
  int k=-1;
  DisjointLists* children=cover->getChildren();
  while(j!=-1) {
    k=getCurrentChild(cover,j);
    int jj=getNext(j);
    int count=0;
    while((k!=-1) && (cover->getLevel(k)==level+1)) {
      //child_ctr++;
      int kk=children->getNext(k);
      double dist=getDist(p,k);
      //cover->nearest_ctr++;
      if(dist< r+lambdaradius) {
	//prepend(k);
	append(k);
      }
      k=kk;
    }
    setCurrentChild(j,k);
    double dist=getDist(p,j);
    if( dist>=r+lambdaradius) {
      remove(j);
    }
    j=jj;
  }

  level++;
  radius*=base;

  return;
}  

bool Cover::checkFindWithin(const Point *p,double r,int maxlevel) const {
  //call after Cover::findWithin
  int count1=0;
  int count2=0;
  for(int i=descend_list.getFirst();i!=-1;i=descend_list.getNext(i)) {
    count1++;
    double dist=points->getDist(p,i);
    int l=getLevel(i);
    if(dist>r) {
      cout << "point i is at dist=" << dist << endl;
      return false;
    }
    if(maxlevel<l) {
      cout << "point i is at level=" << l << endl;
      return false;
    }
  }
  int N=points->getNumber();
  for(int i=0;i<N;i++) {
    double dist=points->getDist(p,i);
    int l=getLevel(i);
    if( (dist<=r) && (l<=maxlevel))
      count2++;
  }
  if(count1!=count2) {
    cout << "count1=" << count1 << " count2=" << count2 << endl;
    return false;
  }
 
  return true;
}


/***************************************************************************/

//findNearest

int Cover::findNearest(const Point* p,int k,int* indices,double* d) {
  //cout << "Cover::findNearestWithin" << endl;
  descend_list.reset();
  descend_list.initialize(this,root);
  for(int l=minlevel;l<maxlevel;l++) {  //will look for children of level m-1
    descend_list.descendForNearest(p,k,indices,d);
  }

  int count=descend_list.getCount();
  sort(indices,d,count);

  k= count<k ? count : k;

  double val=d[k-1];
  
  int K=k;
  for( ;K<count;K++) {
    if(d[K]>val)
      break;
  }

  return K;
}

void Cover::DescendList::descendForNearest(const Point* p,int k,int* indices,
  double* d) {

  // compute U_{j+1}=T_j \cup P^{-1}[T_j]
  int count=0;
  int j=getFirst();
  int l=-1;
  DisjointLists* children=cover->getChildren();
  while(j!=-1) {
    l=getCurrentChild(cover,j);
    int jj=getNext(j);
    int count=0;
    while((l!=-1) && (cover->getLevel(l)==level+1)) {
      int ll=children->getNext(l);
      prepend(l);
      l=ll;
    }
    setCurrentChild(j,l);
    j=jj;
  }

//compute the k-distance

  count=0;
  for(int j=getFirst();j!=-1;j=getNext(j)) {
    indices[count]=j;
    d[count++]=getDist(p,j);
  }
  sort(indices,d,count);

  k= count<k ? count : k;

  double val=d[k-1];

//prune
  double lambdaradius=lambda*radius;
  j=getFirst();
  while(j!=-1) {
    int jj=getNext(j);
    if(getDist(p,j)>val+lambdaradius) {
      remove(j);
    }
    j=jj;
  }


  level++;
  radius*=base;

  return;
}

bool Cover::checkFindNearest(const Point *p, int k,int K,int* indices, 
  double* d,int* Indices, double* D) const {
  
  int Nc=getCoverNumber();
  
  if(K>Nc)
    return false;

  k= k>Nc ? Nc : k;
  if(K<k)
    return false;

  for(int i=k;i<K;i++) {
    if(d[i]!=d[k-1])
      return false;
  }
  
  int N=getNumber();
  int l=getMaxLevel();
  for(int i=0;i<N;i++) {
      Indices[i]=i;
      if(getLevel(i)<=l) {
	D[i]=points->getDist(p,i);
      } else {
	D[i]=DBL_MAX;
      }
  }

  sort(Indices,D,N);

  for(int i=0;i<K;i++) {
    if(d[i]!=D[i]) 
      return false;
  }
  
  return true;
}

/***************************************************************************/
  
void Cover::printOn(int n,ostream& os) {   
  for(int i=0;i<n;i++) {
    os << "X[" << i << "]=";
    points->printOn(i,os);
    os << " level= " << getLevel(i)
       << " parent= " << getParent(i)
       << endl;
  }  
  for(int i=0;i<n;i++) {
    int first=getFirstChild(i);
    os << "children of node " << i << ":" << endl;
    for(int j=first;j!=-1;j=children.getNext(j))
      os << "\t" << j << " level=" << getLevel(j) << endl;
  }
}

#ifdef MEX
bool Cover::checkDistances(double* radii) const {
#else
bool Cover::checkDistances(double* radii,ostream& os) const {
#endif
  radii[0]=radius;
  for(int i=1;i<numlevels;i++) {
    radii[i]=base*radii[i-1];
  }
#ifdef MEX
  mexPrintf("checking distance to parents\n");
#else
  os<< "checking distance to parents" << endl;
#endif
  bool flag=true;
  int N=points->getNumber();
  for(int i=0;i<N;i++) {
    int p=getParent(i);
    if(p!=-1) {
      int l=getLevel(i);
      double val=radii[l-1-minlevel];
      double dist=points->getDist(i,p);
      if(dist>=val) {
	flag=false;
#ifdef MEX
      mexPrintf("level=%d dist=%g\n",l,dist);
      mexPrintf("parent too far i=%d p=%d\n",i,p);
#else 
      cout<< " child " << i << " at level " << getLevel(i) 
	  <<  " is at dist " << dist << " >= " << val
	  << " from parent " << p << " at level " << l << endl;
#endif
      }
    }
  }
#ifdef MEX
  mexPrintf("checking distance between points\n");
#else  
  os<< "checking distance between points" << endl;
#endif

  for(int i=0;i<N;i++) {
    if(!isDuplicate(i)) {
      for(int j=i+1;j<N;j++) {
	if(!isDuplicate(j)) {
	  int li=getLevel(i);
	  int lj=getLevel(j);
	  double dist=points->getDist(i,j);
	  double val = li< lj ? radii[lj-minlevel] : radii[li-minlevel];
	  if(dist<val) {
	    flag=false;
#ifdef MEX 
	    mexPrintf("too close i=%d j=%d dist=%g\n",i,j,dist);
#else
	    os<< i << " at level " << getLevel(i) 
	      << " is at dist " << dist << " from " << j
	      << " at level " << getLevel(j) << " ; should be no less than "
	      << val << endl;
#endif
	  }
	}
      }
    }
  }


  return flag;
}

#ifdef MEX
bool Cover::checkChildren() const {
#else
bool Cover:: checkChildren(ostream& os) const {
#endif
  bool flag=true;
  int N=getNumber();
#ifdef MEX 
  mexPrintf("checking order of children\n");
#else
  os<< "checking order of children" << endl;
#endif

  for(int i=0;i<N;i++) {
    if(!isDuplicate(i)) {
      int first=getFirstChild(i);
      int last=getLastChild(i);
      if(first!=-1) {
	int l=getLevel(first);
	for(int n=children.getNext(first);n!=-1;n=children.getNext(n)) {
	  int k=getLevel(n);
	  if(k<l) {
	    flag=false;
#ifdef MEX 
	    mexPrintf("%d is out of order in children of %d\n",n,i);
#else
	    os<< n << "is out of order in children of " << i << endl;
#endif
	  }
	  if(l<k) 
	    l=k;
	}
      }
    }
  }

  for(int p=0;p<N;p++) {
    if((!isDuplicate(p)) && (getLevel(p)>minlevel)) {
      for(int c=children.getFirst(p);c!=-1;c=children.getNext(c)) {
	if(p!=getParent(c)) {
#ifdef MEX 
	  mexPrintf("problem with parent %d and child %d\n",p,c);
#else
	  os<< "problem with parent " << p << " and child " << c << endl;
#endif
	  flag=false;
	}
      }
    }
  }
  return flag;
}      


Levels::Levels(const Cover* ccover,
  int* plevel_counters, int* plevel_offsets, int* plevels ) 
  : level_counters(plevel_counters) ,level_offsets(plevel_offsets), 
    levels(plevels),
    cover(ccover),
    minlevel(cover->getMinLevel()), maxlevel(cover->getMaxLevel()),
    numlevels(cover->getNumLevels()) {

  int N=cover->getNumber();

  for(int i=0;i<numlevels;i++) {
    level_counters[i]=0;
  }

  for(int i=0;i<N;i++) {
    int l=cover->getLevel(i);
    if(l!=INT_MAX) {
      level_counters[l-minlevel]++;
    }
  }

  level_offsets[0]=0;
  for(int j=1;j<numlevels;j++) {
    level_offsets[j]=level_offsets[j-1]+level_counters[j-1];
  }
  
  for(int i=0;i<N;i++) {
    int l=cover->getLevel(i);
    if(l!=INT_MAX) {
      int j=l-minlevel;
      levels[level_offsets[j]++]=i;
    }
  }
  
 level_offsets[0]=0;
  for(int j=1;j<numlevels;j++) {
    level_offsets[j]=level_offsets[j-1]+level_counters[j-1];
  }
}

bool Levels::checkLevels() const {
  for(int L=0;L<numlevels;L++) {
    int* p=levels+level_offsets[L];
    int K=level_counters[L];
    for(int k=0;k<K;k++) {
      if(cover->getLevel(*(p+k))!=minlevel+L)
	return false;
    }
  }
  return true;
}
 
void Levels::printOn(const Cover* cover,ostream& os) const {

  os << "minlevel=" << minlevel << " maxlevel=" << maxlevel
       << " numlevels=" << numlevels << endl;

  for(int j=0;j<numlevels;j++) {
    os << "level_counters[" << j << "]=" << level_counters[j] << endl;
    for(int k=0;k<level_counters[j];k++) {
      cout << "\tindex=" << levels[level_offsets[j]+k] << endl;
    }
  }
}


void Cover::checkLevels(const Levels* levels,ostream& os) const {
  for(int l=0;l<numlevels;l++) {
    int n=levels->getLevelCounter(l);
    int* p=levels->getLevel(l);
    double val=pow(base,l+minlevel);
    for(int i=0;i<n;i++) {
      int k=p[i];
      for(int j=i+1;j<n;j++) {
	double dist=points->getDist(k,p[j]);
	if(dist<val) {
	  os << "dist(" << i << "," << j << ")=" << dist << endl;
	  os << "\tshould be at least " << val << endl;
	}
      }
    }
  }
}

#include "Lookup.C"
template class Lookup<double>;
