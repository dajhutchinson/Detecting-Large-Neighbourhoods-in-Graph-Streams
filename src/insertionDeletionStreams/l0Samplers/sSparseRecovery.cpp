/*
 *  Implementation of exact s-sparse recovery algorithm defined in 2.3.2 of
 *    A unifying framework for l0-sampling algorithms.
 *
 *  This recovers at most s of the non-zero elements of a vector (Very rarely it will return more)
 *
 *  This implementation is for recovering part of the neighbourhood to the to a specified vertex
 *    when edges are defined in a stream.
 */

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

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

// s-sparse
set<vertex> s_sparse_recovery(string edge_file, vertex target, int num_vertices, int s, double delta);

// 1-sparse
bool verify_1_sparse(int phi,int iota,int tau);
void update_1_sparse_counters(int index,int delta,int row,int col,long** phi_s,long** iota_s,long** tau_s);

// Hashing
hash_params generate_hash(int m);
int hash_function(int key, hash_params ps);

// Utility
void parse_edge(string str, edge& e);
int identify_endpoint(edge e,vertex target);
long** initalise_zero_2d_array(int num_cols, int num_rows);

/*------*
 * BODY *
 *------*/

int main() {
  int s=100; // sparsity to recover
  double delta=.1; // acceptable failure

  // details of graph to perform on
  string file_path="../../../data/facebook_deletion.edges";
  vertex target=2290;
  int num_vertices=747;

  set<vertex> recovered_neighbourhood=s_sparse_recovery(file_path, target, num_vertices, s, delta);
  cout<<recovered_neighbourhood.size()<<endl;
}

/*-------------------*
 * s SPARSE RECOVERY *
 *-------------------*/

set<vertex> s_sparse_recovery(string edge_file, vertex target, int num_vertices, int s, double delta) {
  int num_cols=2*s;
  int num_rows=log(s/delta);

  // arrays for 1-sparse recovery
  long** phi_s =initalise_zero_2d_array(num_cols,num_rows); // sum of weights (sum ai)
  long** iota_s=initalise_zero_2d_array(num_cols,num_rows); // weighted sum of weights (sum ai*i)
  long** tau_s =initalise_zero_2d_array(num_cols,num_rows); // squared weighted sum of weights (sum ai*(i**2))

  // choose hash function for each row
  hash_params* ps_s=new hash_params[num_rows];
  for (int i=0; i<num_rows; i++) ps_s[i]=generate_hash(num_cols);

  string line; edge e; int v;
  ifstream edge_stream(edge_file);

  while (getline(edge_stream,line)) {
    parse_edge(line,e);
    v=identify_endpoint(e,target); // determine if edge is connected to target, if so what is the other endpoint

    if (v!=-1) { // edge is incident to target
      for (int r=0; r<num_rows; r++) { // decide which sampler in each row to update
        int c=hash_function(v,ps_s[r]); // col to update
        update_1_sparse_counters(v,e.value,r,c,phi_s,iota_s,tau_s);
      }
    }
  }

  // check values of 1_sparse
  set<vertex> sampled_neighbourhood;
  for (int i=0; i<num_cols; i++) {
    for (int j=0; j<num_rows; j++) {
      if (verify_1_sparse(phi_s[i][j],iota_s[i][j],tau_s[i][j])) {
        sampled_neighbourhood.insert(iota_s[i][j]/phi_s[i][j]);
  }}}

  return sampled_neighbourhood;
}

/*-------------------*
 * 1 SPARSE RECOVERY *
 *-------------------*/

// update counters with new edge
void update_1_sparse_counters(int index,int delta,int row,int col,long** phi_s,long** iota_s,long** tau_s) {
  phi_s[col][row] +=delta;
  iota_s[col][row]+=delta*index;
  tau_s[col][row] +=delta*pow(index,2);
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
int identify_endpoint(edge e,vertex target) {
  if (e.fst==target) return e.snd;
  else if (e.snd==target) return e.fst;
  else return -1;
}

// initalise 2d array with zero in every index
long** initalise_zero_2d_array(int num_cols, int num_rows) {
  long** arr = (long**) malloc(num_cols * sizeof(long*)); // allocate cols

  for (int i=0; i<num_cols; i++) {
    arr[i]=(long*) malloc(num_rows * sizeof(long)); // allocate rows
    for (int j=0; j<num_rows; j++) arr[i][j]=0; // set values to 0
  }

  return arr;
}
