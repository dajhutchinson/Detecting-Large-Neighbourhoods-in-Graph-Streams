/*
 *  Implementation of Perfect 1-sparse recovery
 *  pg 8 - A unifying Framework for l0-Sampling algorithm
 *
 *  TODO - implement hash function
 *  Hash function from
 *    https://stackoverflow.com/questions/16284317/obtaining-a-k-wise-independent-hash-function
 *    https://github.com/Cyan4973/xxHash
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>

#include "xxhash.h"

using namespace std;

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/
using vertex = int; // typemap vertex
struct edge { // undirected edge
    vertex fst;
    vertex snd;
    int value;
 };

 /*------------*
  * SIGNATURES *
  *------------*/

void parse_edge(string str, edge& e);

/*------*
 * BODY *
 *------*/

int main() {
  ifstream edge_stream("../../../data/facebook_deletion.edges");

  vertex target=2290;

  int sum=0, weighted_sum=0, squared_sum=0, j=1, v=0; // TODO some of these should be bigger than int
  string line; edge e;
  while (getline(edge_stream,line)) {
    parse_edge(line,e);
    // TODO - run for multiple values of j

    if (e.fst==target) { // identify other node on edge (if edge contains target)
      v=e.snd;
    } else if (e.snd==target) {
      v=e.fst;
    } else {
      v=-1;
    }

    if (v!=-1) { // increment sums
      // if (h(v)<=(n^3)/(2^j))
      sum+=e.value;
      weighted_sum+=e.value*v;
      squared_sum+=e.value*pow(v,2);
    }

  }

  // for each j
  // check 1 sparse
  if (squared_sum==sum*weighted_sum) { // is 1-sparse
    cout<<weighted_sum/sum<<endl;
    return weighted_sum/sum;
  }
  // not 1-sparse
  cout<<-1<<endl;
  return -1;
}

/*-----------*
 * UTILITIES *
 *-----------*/

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
  e.fst=stoi(fst);
  e.snd=stoi(snd);
}
