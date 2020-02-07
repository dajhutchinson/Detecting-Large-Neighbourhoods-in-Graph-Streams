#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

int BYTES=0; // space used atm

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/
using time_point=chrono::high_resolution_clock::time_point;
using vertex=string;

struct edge {
  vertex fst;
  vertex snd;
  bool insertion;
};

/*------------*
 * SIGNATURES *
 *------------*/

void parse_edge(string str, edge& e);
void naive_neighbourhood_detection(int d, ifstream& stream, vector<vertex>& neighbourhood, vertex& root);

 /*------*
  * BODY *
  *------*/

int main() {
  int d=297;
  ifstream stream("data/gplus.edges");
  vector<vertex> neighbourhood; vertex root; // variables for returned values

  time_point before=chrono::high_resolution_clock::now(); // time before execution
  naive_neighbourhood_detection(d,stream,neighbourhood,root);
  time_point after=chrono::high_resolution_clock::now(); // time before execution

  // Print out returned neighbourhood, if one exists
  vertex* p=&root;
  if (p==nullptr) cout<<"NO SUCCESSES"<<endl;
  else {
    cout<<"Neighbourhood for <"<<root<<">"<<endl;
    cout<<"<";
    for (vector<vertex>::iterator i=neighbourhood.begin(); i!=neighbourhood.end(); i++) cout<<*i<<",";
    cout<<">"<<endl;
  }
  stream.close();

  auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count();
  cout<<"Execution Time - "<<duration<<" microseconds ("<<duration/1000<<" milliseconds)"<<endl;
  cout<<"Max space - "<<BYTES<<" bytes ("<<BYTES/1024<<" kb)"<<endl;
  return 0;
}

// naive solution for insertion-only stream
void naive_neighbourhood_detection(int d, ifstream& stream, vector<vertex>& neighbourhood, vertex& root) {
  string line; edge e; map<vertex,vector<vertex> > neighbourhoods;
  BYTES+=sizeof(string)+sizeof(edge)+sizeof(map<vertex,vector<int> >);

  while (getline(stream,line)) {
    parse_edge(line,e);

    // increase fst's neighbourhood
    if (neighbourhoods.count(e.fst)) { // already in map
      neighbourhoods[e.fst].push_back(e.snd);
      BYTES+=sizeof(vertex);
    } else { // create neighbourhood
      vector<vertex> v={e.snd};
      neighbourhoods[e.fst]=v;
      BYTES+=sizeof(vector<vertex>)+sizeof(vertex);
    }

    // increase snd's neighbourhood
    if (neighbourhoods.count(e.snd)) {
      neighbourhoods[e.snd].push_back(e.fst);
      BYTES+=sizeof(vertex);
    } else {
      vector<vertex> v={e.fst};
      neighbourhoods[e.snd]=v;
      BYTES+=sizeof(vector<vertex>)+sizeof(vertex);
    }
  }

  for (map<vertex,vector<vertex> >::iterator i=neighbourhoods.begin(); i!=neighbourhoods.end(); i++) {
    if (i->second.size()>=d) { // success
      root=i->first;
      vector<vertex> v=i->second;
      v.erase(v.begin()+d,v.end());
      neighbourhood=v;
      return;
    }
  }

  // fail
  neighbourhood.clear();
  vertex* p=&root;
  p=nullptr;
  return;

}


// parse ege from stream
void parse_edge(string str, edge& e) {
  int spaces=count(str.begin(),str.end(),' ');

  if (spaces==2) { // insertion deletion edge
    if (str[0]=='D') e.insertion=false;
    else e.insertion=true;
    str=str.substr(2,str.size());
  } else if (spaces==1) { // insertion edge
    e.insertion=true;
  }

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

  // Update edge values
  e.fst=fst;
  e.snd=snd;
}
