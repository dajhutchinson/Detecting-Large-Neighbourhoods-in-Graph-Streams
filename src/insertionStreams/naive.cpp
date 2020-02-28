/*
 *  Anaive implementation of neighbourhood detection algorithm for insertion-only Streams
 *  See algorithm 2 in paper
 */

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <set>

using namespace std;

int BYTES=0; // space used atm

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/
using time_point=chrono::high_resolution_clock::time_point;
using vertex=string;

struct edge { // undirected edge
  vertex fst;
  vertex snd;
};

/*------------*
 * SIGNATURES *
 *------------*/

void update(vertex v1, vertex v2, map<vertex,int>& degrees, map<vertex,set<vertex> >& neighbourhoods);
int naive(string edge_file_path, int c, int d, vertex& root, set<vertex>& neighbourhood);
void parse_edge(string str, edge& e);

/*------*
 * BODY *
 *------*/

int main() {
  int d=586; int c=20; // d=max_degree,c=accuracy
  string edge_file_path="../../data/facebook.edges";

  vertex root;
  set<vertex> neighbourhood;

  time_point before=chrono::high_resolution_clock::now();
  int edge_count=naive(edge_file_path,c,d,root,neighbourhood);
  time_point after=chrono::high_resolution_clock::now();

  auto duration=chrono::duration_cast<chrono::microseconds>(after-before).count();

  if(edge_count!=-1) {
    cout<<"root="<<root<<endl<<"# edges checked="<<edge_count<<endl<<"Neighbourhood ("<<neighbourhood.size()<<")"<<endl; // output result
    cout<<"time="<<duration<<"mus"<<endl<<"space used="<<BYTES/1024<<"mb"<<endl<<"{";
    for (set<vertex>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) cout<<*it<<",";
    cout<<"}"<<endl;
  }


  return -1;
}

int naive(string edge_file_path, int c, int d, vertex& root, set<vertex>& neighbourhood) {
  map<vertex,int> degrees;
  map<vertex,set<vertex> > neighbourhoods;
  BYTES+=sizeof(map<vertex,int>)+sizeof(map<vertex,set<vertex> >);

  ifstream stream(edge_file_path);
  string line; edge e; int edge_count=0;
  BYTES+=sizeof(string)+sizeof(edge)+sizeof(int);
  while (getline(stream,line)) {
    edge_count+=1;
    parse_edge(line,e);

    // update degree count & neighbourhoods
    update(e.fst,e.snd,degrees,neighbourhoods);
    update(e.snd,e.fst,degrees,neighbourhoods);

    // check if updated neighbourhoods meet targets
    if (neighbourhoods[e.fst].size()>=d/c) { // if meet targets
      root=e.fst;neighbourhood=neighbourhoods[e.fst];
      return edge_count; // end program
    }

    if (neighbourhoods[e.snd].size()>=d/c) { // if meet targets
      root=e.snd;neighbourhood=neighbourhoods[e.snd];
      return edge_count; // end program
    }

  }

  return -1; //fail
}

// update degree count & neighbourhoods
void update(vertex v1, vertex v2, map<vertex,int>& degrees, map<vertex,set<vertex> >& neighbourhoods) {
  if (degrees.count(v1)) {
    degrees[v1]+=1;
    neighbourhoods[v1].insert(v2);
    BYTES+=sizeof(vertex);
  } else {
    degrees[v1]=1;
    BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*);
    set<vertex> new_set;
    new_set.insert(v2);
    neighbourhoods[v1]=new_set;
    BYTES+=(2*sizeof(vertex))+sizeof(set<vertex>)+sizeof(void*);
  }
}

/*---------*
 * UTILITY *
 *---------*/

// parse ege from stream
void parse_edge(string str, edge& e) {
  string fst="",snd="";
  bool after=false;

  for (char& c:str) {
    if (c==' ') { // seperator
      after=true;
    } else if (after) { // second id
      snd+=c;
    } else { // first id
      fst+=c;
    }
  }

  // Set edge end-points
  e.fst=fst;e.snd=snd;
}
