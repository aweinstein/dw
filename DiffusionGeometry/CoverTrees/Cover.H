#ifndef COVER_H
#define COVER_H
#include <iostream>
#include "DLList.H"
#include "DisjointLists.H"
#include "Lookup.H"
#include "Points.H"
#include <float.h>
#include <limits.h>

#ifdef MEX
#include "mex.h"
#endif

using namespace std;

class Levels;

class Cover {
private:
  const Points* points;

  int root;

  int cover_number;
  int num_duplicates;

  int minlevel; //always <=0
  int maxlevel;
  int numlevels; //maxlevel-minlevel+1
  int maxdescend; //the maximum number of levels you will descend
                  
  double base; 
  double radius; //this is base^{-minlevel}; it's dynamic

  int* lp; // level parent
  //if root, level=minlevel (this is dynamic) parent=-1
  //if duplicate, level=INT_MAX parent=index matched
  //if too deep according to maxdescend level=INT_MAX parent=-1
 

  //children
  DisjointLists children;
  
  //int nearest_ctr;

  bool setMinLevelAndBase(double);
  void setMaxLevel(int maxl) { if(maxlevel<maxl) maxlevel=maxl; } 
  
  //auxiliary doubly linked list
public:
  class DescendList : public DLList {
  private:
    Cover* cover;
    const Points* points;
    int level;
    double base;
    double radius; //Not to be confused with Cover::radius
    double lambda;
    int parent;
    int parent_level;

    //in order not to compute distances more than once.
    //keep track of distances for descent
    //keep track of the current child for descent
    //static int child_ctr;

    Lookup<double> dist_lookup;
    Lookup<int> current_child_lookup;
  public:
    double getDist(int,int);
    double getDist(const Point*,int);

    DescendList(double,
		Cover*,
		const Points*,
		int*,
		int*,
		int*,
		double*, 
		int*,
		int*,
		int*);
    ~DescendList();
    int getLevel() const { return level; }
    int getParentLevel() const { return parent_level; }
    double getRadius() const { return radius; }
    int getParent() const { return parent; }
    void setRadius(double r) { radius=r; }
    int getCurrentChild(Cover*,int i);
    void setCurrentChild(int,int);

    int getNumberTrueFlags() const { return dist_lookup.getNumberTrueFlags(); }
    void reset();
    //double getDistForInsertion(const Points*,int,int);
    int descend(int);
    void descendForWithin(const Point*,double);
    void descendForNearest(const Point*,int,int*,double*);
    
    void initialize(const Cover*,int);

    int getDistNCallsToGet() const { 
      return dist_lookup.getNCallsToGet(); }
    int getDistNCallsToSet() const { 
      return dist_lookup.getNCallsToSet(); }
    int getChildrenNCallsToGet() const { 
      return current_child_lookup.getNCallsToGet(); }
    int getChildrenNCallsToSet() const { 
      return current_child_lookup.getNCallsToSet(); }

  };

private:
  DescendList descend_list;

public:
  Cover(double base, 
	int maxdescend,
	const Points* ppoints,  //N=ppoints->getNumber()
	int* plp, // 2*N; level parents
	int* pchildren, // 4*N  children DisjointLists
	int* p, // 2*N space for DLList in descend_list
	int* p1, //N indices_to_dist_flags
	int* p2, //N dist_flags
	double* p3, //N distances
	int* r1, //N current_child_flags
	int* r2, //N indices_to_current_child_flags
	int* r3); //N current_children

  Cover(double ibase,
	int* params,
	const Points* pPoints,
	int* plp, 
	DisjointLists&  rchildren, 
	int* paddress,
	int* pdist_flags,
	int* pindices_to_dist_flags,
	double* pdistances,
	int* pcurrent_child_flags,
	int* pindices_to_current_child_flags,
	int* pcurrent_child);


  ~Cover();
  int getNumber() const { return points->getNumber(); }
  int getCoverNumber() const { return cover_number; }
  int getNumDuplicates() const { return num_duplicates; }
  bool isDuplicate(int i) const { return getLevel(i)==INT_MAX; }
  const Points* getPoints() const { return points; }
  void setLevel(int i,int l) { lp[2*i]=l; }
  void setParent(int i,int p) { lp[2*i+1]=p; }
  void makeNewEntry(int index,int level,int parent);
  double getRadius() const { return radius; }
  int getLevel(int i) const { return lp[2*i]; }
  int getParent(int i) const { return lp[2*i+1]; }
  
  DisjointLists* getChildren() { return &children; }
 
  int getFirstChild(int i) const { return children.getFirst(i); }
  int getLastChild(int i) const { return children.getLast(i); }
  int getMinLevel() const { return minlevel; }
  int getMaxLevel() const { return maxlevel; }
  int getMaxDescend() const { return maxdescend; }
  int getNumLevels() const { return numlevels; }
  int getRoot() const { return root; }
  double getBase() const { return base; }
  void newRoot(double,int);
  void insert(int);
  void appendToChildren(int parent,int child);
  void prependToChildren(int parent,int child);
  void append(int index,int level, int parent,int index_into_children);
  void insertInChildrenAtRightLevel(int,int,int);
  const Cover::DescendList& getDescendList() const { return descend_list; }
  void clearDescendList() { descend_list.clear(); }
  void fillArrFromDescendList(int*) const;
  void findWithin(const Point*,double r,int l);
  bool checkFindWithin(const Point*,double,int) const;
  int findNearest(const Point*,int,int*,double*);
  bool checkFindNearest(const Point *p, int k,int K, int* indices, double* d,
    int* Indices, double* D) const;
#ifdef MEX
  bool checkDistances(double* radii) const;
#else
  bool checkDistances(double* radii,ostream& os = cout) const;
#endif
#ifdef MEX
  bool checkChildren() const;
#else
  bool checkChildren(ostream& os = cout) const;
#endif
  void printOn(int n,ostream& os);
  void checkLevels(const Levels*,ostream&) const;
  int getDistNCallsToGet() const { 
    return descend_list.getDistNCallsToGet(); }
  int getDistNCallsToSet() const { 
    return descend_list.getDistNCallsToSet(); }
  int getChildrenNCallsToGet() const { 
    return descend_list.getChildrenNCallsToGet(); }
  int getChildrenNCallsToSet() const { 
    return descend_list.getChildrenNCallsToSet(); }
};

class Levels {
private:
  const Cover* cover;
  int N;
  int minlevel, maxlevel, numlevels;
  int* level_counters;
  int* level_offsets;
  int* levels;
public:
  Levels(const Cover*,int *,int*,int*);
  int* getLevelOffsets() const { return level_offsets; }
  int getLevelCounter(int l) const { return level_counters[l]; }
  int* getLevel(int l) const { return levels+level_offsets[l]; }
  int* getLevels() const { return levels; }
  bool checkLevels() const;
  void printOn(const Cover*,ostream& = cout) const;
};

#endif
