/*
 *  A bastardisation of 1-sparse recovery st it makes an array 1-sparse (although sometimes failing)
 *  Achieved by applying a hash function to the vertex value & only incrementing the counters if the
 *    hash value falls in a given bound. This bound is tightened until array is 1-sparse (or zero array).
 *
 *  Uses O(log_2 n) space (n=num vertices)
 */

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <vector>

//#include "xxhash.hpp"

using namespace std;

int P=1073741789; // >2^30

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/
using vertex = int; // typemap vertex

struct edge { // undirected edge
    vertex fst;
    vertex snd;
    int value;
 };

struct hash_params { // parameters for hash function
  unsigned long a;
  unsigned long b;
  unsigned long m;
};

/*------------*
 * SIGNATURES *
 *------------*/

int one_sparse_recovery(string edge_file, vertex target, int num_vertices);
bool verify_1_sparse(int phi,int iota,int tau);
void update_counters(int index,int delta,int j,int* phi_s,int* iota_s,int* tau_s);

// Hashing
hash_params generate_hash(int m);
int hash_function(int key, hash_params ps);

// Utility
void parse_edge(string str, edge& e);
int identify_endpoint(edge e,vertex target);

/*------*
 * BODY *
 *------*/

int main() {
  //ifstream edge_stream("../../../data/facebook_deletion.edges");

  string file_path="../../../data/facebook_deletion.edges";
  vertex target=2290;
  int num_vertices=747;

  cout<<one_sparse_recovery(file_path,target,num_vertices)<<endl;
}

// recovers a neighbour of target (-1 if failure)
int one_sparse_recovery(string edge_file, vertex target, int num_vertices) {

  ifstream edge_stream(edge_file);
  int j_lim=log2(pow(num_vertices,3)); // num of js run
  hash_params ps=generate_hash(pow(num_vertices,3)); // choose hash function

  // arrays for checks
  int* phi_s =new int[j_lim+1]; // sum of weights (sum ai)
  int* iota_s=new int[j_lim+1]; // weighted sum of weights (sum ai*i)
  int* tau_s =new int[j_lim+1]; // squared weighted sum of weights (sum ai*(i**2))
  int* hash_lims=new int[j_lim+1]; // hash limit for each value of j

  int num_vertices_3=pow(num_vertices,3);
  for (int i=1; i<j_lim+1; i++) { // initalise to zero arrays & calculate hash limits
    phi_s[i]=0; iota_s[i]=0; tau_s[i]=0;
    hash_lims[i]=num_vertices_3/pow(2,i);
  }

  string line; edge e; int v;

  while (getline(edge_stream,line)) {
    parse_edge(line,e);
    v=identify_endpoint(e,target);

    if (v!=-1) { // target is on edge
      int h=hash_function(v,ps);
      for (int j=1; j<j_lim+1; j++) { // this is the bastardisation
        if (h<=hash_lims[j]) update_counters(v,e.value,j,phi_s,iota_s,tau_s);
      }
    }

  }

  // check 1 sparsity
  for (int j=1; j<j_lim+1; j++) {
    if (verify_1_sparse(phi_s[j],iota_s[j],tau_s[j])) return iota_s[j]/phi_s[j];
  }

  return -1; // no 1-sparse
}

// returns endpoint which is not target (or -1 if target not on edge)
int identify_endpoint(edge e,vertex target) {
  if (e.fst==target) return e.snd;
  else if (e.snd==target) return e.fst;
  else return -1;
}

// update counters with new edge
void update_counters(int index,int delta,int j,int* phi_s,int* iota_s,int* tau_s) {
  phi_s[j] +=delta;
  iota_s[j]+=delta*index;
  tau_s[j] +=delta*pow(index,2);
}

// verify if array is 1_sparse
bool verify_1_sparse(int phi,int iota,int tau) {
  if (pow(iota,2)==phi*tau & phi!=0) return true;
  return false;
}

/*---------*
 * HASHING *
 *---------*/

// generate parameters to use in hash function
hash_params generate_hash(int m) {
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
  uniform_int_distribution<unsigned long> distribution(0,P-1);
  hash_params ps={
    distribution(generator),
    distribution(generator),
    (unsigned long)m
  };
  return ps;
}

// hash a key
int hash_function(int key, hash_params ps) {
  return ((ps.a*key+ps.b)%P)%ps.m;
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
