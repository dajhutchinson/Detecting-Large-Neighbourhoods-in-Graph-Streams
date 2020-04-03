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

void execute_test(int c_min, int c_max, int c_step, int d, string in_file, string out_file);
int naive(string edge_file_path, int c, int d, vertex& root, set<vertex>& neighbourhood);
void update(vertex v1, vertex v2, map<vertex,int>& degrees, map<vertex,set<vertex> >& neighbourhoods);
void parse_edge(string str, edge& e);

/*------*
 * BODY *
 *------*/

int main() {
  int d=5948; int c=20; // d=max_degree,c=accuracy
  string edge_file_path="../../data/gplus.edges";
  string output_file_path="../../results/gplus_naive_results.csv";

  execute_test(2,100,1,d,edge_file_path,output_file_path);

  return -1;
}

// deterministic so dont need repetitions
void execute_test(int c_min, int c_max, int c_step, int d, string in_file, string out_file) {
  ofstream outfile(out_file);
  outfile<<"name,"<<in_file<<endl<<"d,"<<d<<endl<<endl<<endl;
  outfile<<"c,time (microseconds),space (bytes),mean egdes checked"<<endl;

  vertex root; set<vertex> neighbourhood;
  for (int c=c_min;c<=c_max;c+=c_step) {
    cout<<c<<endl;
    BYTES=0;
    vertex* p=&root; p=nullptr;
    neighbourhood.clear();

    ifstream stream(in_file);

    time_point before=chrono::high_resolution_clock::now();
    int edge_count=naive(in_file,c,d,root,neighbourhood);
    time_point after=chrono::high_resolution_clock::now();

    stream.close();
    auto duration=chrono::duration_cast<chrono::microseconds>(after-before).count();
    if (edge_count!=-1) {
      outfile<<c<<","<<duration<<","<<BYTES<<","<<edge_count<<endl;
    } else {
      cout<<c<<" FAIL";
      return;
    }
  }
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
