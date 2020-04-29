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
  //int d=586, n=747, reps=100; int c=3;
  //display_results(c,d,n,"../../data/facebook.edges");

  string out_file_path="results_naive.csv";
  int n,d,reps; string edge_file_path;

  // c=runs, d/c=d2, n=# vertices, NOTE - set d=max degree, n=number of vertices
  //n=12417; d=5948; reps=10; edge_file_path="../../data/gplus.edges"; // NOTE - # edges=1,179,613
  //n=102100; d=104947; reps=10; edge_file_path="../../data/gplus_large.edges"; // NOTE - # edges=30,238,035
  //n=52; d=35; reps=10; edge_file_path="../../data/facebook_small.edges"; // NOTE - # edges=292
  //n=747; d=586; reps=10; edge_file_path="../../data/facebook.edges"; // NOTE - # edges=60,050
  //n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000star.edges"; // NOTE - # edges=999
  //n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000complete.edges"; // NOTE - # edges=499,500
  //execute_test(3,20,1,d,edge_file_path,out_file_path);

  n=52; d=35; reps=10; edge_file_path="../../data/facebook_small.edges"; // NOTE - # edges=292
  out_file_path="0_results_naive_facebook_small.csv";
  execute_test(3,20,1,d,edge_file_path,out_file_path);

  n=747; d=586; reps=10; edge_file_path="../../data/facebook.edges"; // NOTE - # edges=60,050
  out_file_path="0_results_naive_facebook.csv";
  execute_test(3,20,1,d,edge_file_path,out_file_path);

  n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000star.edges"; // NOTE - # edges=999
  out_file_path="0_results_naive_1000star.csv";
  execute_test(3,20,1,d,edge_file_path,out_file_path);

  n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000complete.edges"; // NOTE - # edges=499,500
  out_file_path="0_results_naive_1000complete.csv";
  execute_test(3,20,1,d,edge_file_path,out_file_path);

  n=12417; d=5948; reps=10; edge_file_path="../../data/gplus.edges"; // NOTE - # edges=1,179,613
  out_file_path="0_results_naive_gplus.csv";
  execute_test(3,20,1,d,edge_file_path,out_file_path);

  return -1;
}

// deterministic so dont need repetitions
void execute_test(int c_min, int c_max, int c_step, int d, string in_file, string out_file) {
  ofstream outfile(out_file);
  outfile<<"name,"<<in_file<<endl<<"d,"<<d<<endl<<endl<<endl;
  outfile<<"c,time (microseconds),space (bytes),mean egdes checked"<<endl;

  vertex root; set<vertex> neighbourhood;
  for (int c=c_min;c<=c_max;c+=c_step) {
    cout<<c<<"/"<<c_max<<" "<<in_file<<" naive"<<endl;
    BYTES=0;
    vertex* p=&root; p=nullptr;
    neighbourhood.clear();

    ifstream stream(in_file);

    time_point before=chrono::high_resolution_clock::now();
    int edge_count=naive(in_file,c,d,root,neighbourhood);
    time_point after=chrono::high_resolution_clock::now();
    cout<<"\r                                                \r";

    stream.close();
    auto duration=chrono::duration_cast<chrono::microseconds>(after-before).count();
    cout<<duration/1000000<<"s"<<endl<<endl;
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
    if (edge_count%10000==0) cout<<"\r"<<edge_count;
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
