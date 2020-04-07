/*
 *  TEST DESCRIPTION
 *    In this implementation the hash function assigns values uniformly at random to a unique value.
 *    This is inefficient as it has to store where every value goes but it does mean that it is k-wise independent.
 *    This is to test whether using a pairwise independent hash function was the problem.
 *    If this comes out as uniform, then the pairwise independent hash function was not sufficient.
 *    NOTE this only concerns the hash function used to choose whether to pass update to an s-sparse recovery
 *        not the one used during s-sparse recovery
 *
 *  An implementation of an L0-Sampler for a specified vertex in a graph stream.
 *    This samples uniformly from the neighbourhood of the vertex.
 *    Edges are provided as a stream of inputs
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

// L0 Sampling
vertex l0sampling(string file_path, vertex target, int num_vertices, double delta, double gamma);

// s-sparse
hash_params* choose_hash_functions(int num_cols, int num_rows);
void update_s_sparse(vertex endpoint, int edge_value, int num_rows, hash_params* ps_s, long** phi_s, long** iota_s, long** tau_s);
set<vertex> recover_neighbourhood(int num_cols, int num_rows, long** phi_s, long** iota_s, long** tau_s);
vertex recover_vertex(set<vertex> neighbourhood, int sparsity, map<int,int> hash_map);

// 1-sparse
bool verify_1_sparse(int phi,int iota,int tau);
void update_1_sparse_counters(int index,int delta,int row,int col,long** phi_s,long** iota_s,long** tau_s);

// Hashing
hash_params generate_hash(int m);
int hash_function(int key, hash_params ps);
// n = num keys, m=possible hash values, hash=map from key to hash value
void generate_random_hash(int n, int m, map<int,int>& hash);

// Utility
void parse_edge(string str, edge& e);
int identify_endpoint(edge e,vertex target);
long*** initalise_zero_3d_array(int depth, int num_cols, int num_rows);
hash_params** initalise_2d_hash_params_array(int num_cols, int num_rows);
void free_3d_long_array(long*** arr, int depth, int num_cols, int num_rows);
void free_2d_hash_params_array(hash_params** arr, int num_cols, int num_rows);
void write_to_file(string outfile_path, map<vertex,int>& edge_count);

/*------*
 * BODY *
 *------*/

int main() {
  // details of graph to perform on
  string file_path="../../../data/artifical/test_deletion_double.edges";
  vertex target=779;
  int num_vertices=1000;

  // constraints on s-sparse recovery
  double delta=.01; // acceptable failure for L0
  double gamma=.01; // acceptable failure for s-sparse recovery

  // check if edges sampled ~ uniformly
  map<vertex,int> edge_count;

  int succ_count=0;
  int num_runs=10000;
  for (int run=0; run<num_runs; run++) {

    vertex sampled_vertex=l0sampling(file_path,target,num_vertices,delta,gamma);

    // output & update results
    if (sampled_vertex==-1) cout<<run<<",FAIL"<<endl; // recover failed
    else { // recovery succeed
      cout<<run<<",SUCCESS("<<sampled_vertex<<")"<<endl;
      succ_count+=1;
      if (edge_count.count(sampled_vertex)) edge_count[sampled_vertex]+=1; // update number of sample occurences
      else edge_count[sampled_vertex]=1;
    }

  }
  cout<<endl<<succ_count<<"/"<<num_runs<<endl;

  write_to_file("uniform_sample_test.csv", edge_count);
}

/*-------------*
 * L0 Sampling *
 *-------------*/

vertex l0sampling(string file_path, vertex target, int num_vertices, double delta, double gamma) {
  // calculate parameters
  int s=1/delta; // sparsity to recover at
  int j=log2(num_vertices); // number of s-sparse recoveries to run
  int num_cols=2*s;
  int num_rows=log(s/gamma);

  // generate hash map for each vertex (k-independent)
  map<int,int> unique_hash_map;
  generate_random_hash(num_vertices,pow(num_vertices,3),unique_hash_map);

  // arrays for 1-sparse recovery
  long*** phi_s =initalise_zero_3d_array(j,num_cols,num_rows); // sum of weights (sum ai)
  long*** iota_s=initalise_zero_3d_array(j,num_cols,num_rows); // weighted sum of weights (sum ai*i)
  long*** tau_s =initalise_zero_3d_array(j,num_cols,num_rows); // squared weighted sum of weights (sum ai*(i**2))

  // choose hash function for each row
  hash_params** ps_s=initalise_2d_hash_params_array(j,num_rows); // each col is for one s-sparse recovery
  for (int i=0; i<j; i++) ps_s[i]=choose_hash_functions(num_cols,num_rows);

  // hash limits for each value of j
  int* hash_lims=new int[j];
  int n3=pow(num_vertices,3);
  for (int i=1; i<=j; i++) hash_lims[i]=n3/pow(2,i);

  // Process stream
  string line; edge e; int v;
  int sparsity_estimate=0; // r in survery paper algorithm 2
  ifstream edge_stream(file_path);

  // run through stream
  while (getline(edge_stream,line)) {
    parse_edge(line,e);
    v=identify_endpoint(e,target); // determine if edge is connected to target, if so what is the other endpoint

    if (v!=-1) { // edge is connected to target vertex
      sparsity_estimate+=e.value; // increment/decrement depending upon insertion or deletion edge
      for (int i=0; i<j; i++) { // for each s-sparse recovery
        int h=unique_hash_map[v];
        if (h<=hash_lims[i]) update_s_sparse(v,e.value,num_rows,ps_s[i],phi_s[i],iota_s[i],tau_s[i]);
    }}

  }

  // Gather sample from j_sample^th s-sparse recovery
  int j_sample=log2(sparsity_estimate)-1; // -1 since 0 indexed
  set<vertex> sampled_neighbourhood=recover_neighbourhood(num_cols,num_rows,phi_s[j_sample],iota_s[j_sample],tau_s[j_sample]);
  vertex sampled_vertex=recover_vertex(sampled_neighbourhood,s,unique_hash_map);

  // free space
  free_3d_long_array(phi_s,j,num_cols,num_rows);
  free_3d_long_array(iota_s,j,num_cols,num_rows);
  free_3d_long_array(tau_s,j,num_cols,num_rows);
  free_2d_hash_params_array(ps_s,j,num_rows);

  return sampled_vertex;
}

/*-------------------*
 * s-SPARSE RECOVERY *
 *-------------------*/

// choose hash function for each row of s-sparse recovery
hash_params* choose_hash_functions(int num_cols, int num_rows) {
  hash_params* ps_s=new hash_params[num_rows];
  for (int i=0; i<num_rows; i++) ps_s[i]=generate_hash(num_cols);
  return ps_s;
}

// update 1-sparse counters of the s-sparse recovery
void update_s_sparse(vertex endpoint, int edge_value, int num_rows, hash_params* ps_s, long** phi_s, long** iota_s, long** tau_s) {
  for (int r=0; r<num_rows; r++) { // decide which sampler in each row to update
    int c=hash_function(endpoint,ps_s[r]); // col to update
    update_1_sparse_counters(endpoint,edge_value,r,c,phi_s,iota_s,tau_s);
  }
}

// recover neighbourhood from s-sparse recovery counters
set<vertex> recover_neighbourhood(int num_cols, int num_rows, long** phi_s, long** iota_s, long** tau_s) {
  set<vertex> neighbourhood;
  for (int c=0; c<num_cols; c++) {
    for (int r=0; r<num_rows; r++) {
      if (verify_1_sparse(phi_s[c][r],iota_s[c][r],tau_s[c][r])) {
        neighbourhood.insert(iota_s[c][r]/phi_s[c][r]);
  }}}
  return neighbourhood;
}

// recover vertex from recovered neighbourhood
vertex recover_vertex(set<vertex> neighbourhood, int sparsity, map<int,int> hash_map) {
  if (neighbourhood.size()>sparsity || neighbourhood.size()==0) return -1; // s-sparse recovery failed
  else { // return vertex in neighbourhood with min hash value
    int min_hash=INT_MAX, min_val=-1;
    for (set<vertex>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) {
      int h_i=hash_map[*it];
      if (h_i<min_hash) { // lowest yet
        min_hash=h_i;
        min_val=*it;
      }
    }
    return min_val;
  }
}

/*-------------------*
 * 1-SPARSE RECOVERY *
 *-------------------*/

// update counters with new edge
void update_1_sparse_counters(int index,int delta,int row,int col,long** phi_s,long** iota_s,long** tau_s) {
  try {
    phi_s[col][row] +=delta;
    iota_s[col][row]+=delta*index;
    tau_s[col][row] +=delta*pow(index,2);
  } catch (exception exp) {
    cout<<"***********************************VERIFY 1 SPARSE COUNTERS FAILURE";
  }
}

// verify if array is 1_sparse
bool verify_1_sparse(int phi,int iota,int tau) {
  try {
    if (pow(iota,2)==phi*tau && phi!=0) return true;
    return false;
  } catch (exception exp) {
    cout<<"***********************************VERIFY 1 SPARSE FAILURE";
    return false;
  }
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

// generate random hash with unique values for all keys
void generate_random_hash(int n, int m, map<int,int>& hash) {
  vector<int> used; // record hash values which have been used
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
  uniform_int_distribution<unsigned long> distribution(0,m);
  for (int i=1; i<=n; i++) {
    // find unique hash_value
    int hash_value=distribution(generator);
    while (find(used.begin(), used.end(), hash_value)!=used.end()) hash_value=distribution(generator);

    hash[i]=hash_value;         // store in hash map
    used.push_back(hash_value); // record hash_value as used
  }
}

/*-----------*
 * UTILITIES *
 *-----------*/

// parse edge from string "v1 v2" or "(I/D) v1 v2"
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

 // allocate space of 3d array of long intergers, all values set of 0
long*** initalise_zero_3d_array(int depth, int num_cols, int num_rows) {
  long*** arr = (long***) malloc(depth * sizeof(long**)); // allocate depth
  for (int i=0; i<depth; i++) {
    arr[i] = (long**) malloc(num_cols * sizeof(long*)); // allocate cols

    for (int j=0; j<num_cols; j++) {
      arr[i][j]=(long*) malloc(num_rows * sizeof(long)); // allocate rows
      for (int k=0; k<num_rows; k++) arr[i][j][k]=0; // set values to 0
    }

  }

  return arr;
}

// allocate space of 2d array of hash parameters
hash_params** initalise_2d_hash_params_array(int num_cols, int num_rows) {
  hash_params** arr = (hash_params**) malloc(num_cols * sizeof(hash_params*)); // allocate cols
  // allocate rows
  for (int i=0; i<num_cols; i++) arr[i]=(hash_params*) malloc(num_rows * sizeof(hash_params));
  return arr;
}

// free space of 3d array of long integers
void free_3d_long_array(long*** arr, int depth, int num_cols, int num_rows) {
  for (int d=0; d<depth; d++) {
    for (int c=0; c<num_cols; c++) {
      long* ptr2=arr[d][c];
      free(ptr2);
    }
    long** ptr=arr[d];
    free(ptr);
  }
}

// free space of 2d array of hash parameters
void free_2d_hash_params_array(hash_params** arr, int num_cols, int num_rows) {
  for (int c=0; c<num_cols; c++) {
    hash_params* ptr=arr[c];
    free(ptr);
  }
}

// write results to a file
void write_to_file(string outfile_path, map<vertex,int>& edge_count) {
  ofstream outfile(outfile_path);
  outfile<<"vertex,count"<<endl;
  for(map<vertex,int>::iterator it=edge_count.begin(); it!=edge_count.end(); it++) { // write num sampled values to file
    outfile<<it->first<<","<<it->second<<endl;
  }
  outfile.close();
}
