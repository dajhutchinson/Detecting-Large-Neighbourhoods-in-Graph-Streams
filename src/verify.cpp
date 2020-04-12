/*
 *  Verifies that a node has a given neighbourhood
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

using namespace std;

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

using vertex=int;

struct edge { // undirected edge
  vertex fst;
  vertex snd;
  int value;
};

/*------------*
 * SIGNATURES *
 *------------*/

// verify
bool verify_insertion_only(vertex target, set<vertex> neighbourhood, string edge_file_path);
bool verify_insertion_deletion(vertex target, set<vertex> neighbourhood, string edge_file_path);

// utility
void parse_edge(string str, edge& e);
vertex identify_endpoint(edge e,vertex target);

/*------*
 * BODY *
 *------*/

int main() {
  vertex vertices1[]={5,4,3,2}; // 4 is not
  set<vertex> neighbourhood1(vertices1,vertices1+3);
  cout<<boolalpha<<verify_insertion_only(1,neighbourhood1,"../data/artifical/1000vertices.edges")<<endl;

  vertex vertices2[]={2,6,8}; // 2 is not
  set<vertex> neighbourhood2(vertices2,vertices2+3);
  cout<<boolalpha<<verify_insertion_deletion(1,neighbourhood2,"../data/artifical/test_deletion.edges")<<endl;
  return 0;
}

/*--------*
 * VERIFY *
 *--------*/

bool verify_insertion_only(vertex target, set<vertex> neighbourhood, string edge_file_path) {
  string line; edge e; vertex v;
  ifstream edge_stream(edge_file_path);

  int edge_counter=0;
  while (getline(edge_stream,line)) {
    edge_counter+=1;
    if (edge_counter%10000==0) cout<<"\r"<<edge_counter;

    parse_edge(line,e);
    vertex v=identify_endpoint(e,target);

    // if edge is connected to target
    if (v!=-1) neighbourhood.erase(v); // no longer need to find
    if (neighbourhood.size()==0) return true; // all neighbours found
  }
  cout<<endl;
  return false;
}

bool verify_insertion_deletion(vertex target, set<vertex> neighbourhood, string edge_file_path) {
  map<vertex,int> counters;
  for (set<vertex>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) counters[*it]=0;

  string line; edge e; vertex v;
  ifstream edge_stream(edge_file_path);

  int edge_counter=0;
  while (getline(edge_stream,line)) { // NOTE has to run through whole stream
    edge_counter+=1;
    if (edge_counter%10000==0) cout<<"\r"<<edge_counter;

    parse_edge(line,e);
    vertex v=identify_endpoint(e,target);

    if (v!=-1) { // if edge is connected to target
      if (neighbourhood.find(v)!=neighbourhood.end()) { // endpoint is in desired neighbourhood
        counters[v]+=e.value;
    }}
  }
  cout<<endl;

  for (set<vertex>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) {
    if (counters[*it]!=1) return false; // if vertex not in final graph neighbourhood
  }

  return true;
}

/*---------*
 * UTILITY *
 *---------*/

// parse ege from stream
void parse_edge(string str, edge& e) {
  int spaces=count(str.begin(),str.end(),' ');

  if (spaces==2) { // insertion deletion edge
    if (str[0]=='D') e.value=-1;
    else e.value=1;
    str=str.substr(2,str.size());
  } else if (spaces==1) { // insertion edge
    e.value=1;
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
  try {
    e.fst=stoi(fst);
    e.snd=stoi(snd);
  } catch (exception ex) {
    e.fst=-1;
    e.snd=-1;
    e.value=0;
  }
}

// returns endpoint which is not target (or -1 if target not on edge)
vertex identify_endpoint(edge e,vertex target) {
  if (e.fst==target) return e.snd;
  else if (e.snd==target) return e.fst;
  else return -1;
}
