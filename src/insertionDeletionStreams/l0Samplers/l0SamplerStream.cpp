/*
 *  DESCRIPTION
 *    An implementation of the L0-Sampler from l0Sampler.cpp on a stream
 *      st multiple samplers run on a single pass of the stream.
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

// L0 Stream
void initialise_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long****& phi_s, long****& iota_s, long****& tau_s);
void free_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long**** phi_s, long**** iota_s, long**** tau_s);

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
void generate_random_hash(int n, int m, map<vertex,int>& hash);

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

  int num_samplers=100000;

  // prepare samplers
  // same s,j,num_cols,num_rows
  int s=1/delta; // sparsity to recover at
  int j=log2(num_vertices); // number of s-sparse recoveries to run
  int num_cols=2*s;
  int num_rows=log(s/gamma);

  cout<<"# Samplers:"<<num_samplers<<endl<<"# s-sparse:"<<j<<endl<<"# cols:"<<num_cols<<endl<<"# rows:"<<num_rows<<endl;

  // allocate space for counters
  long**** phi_s; long**** iota_s; long**** tau_s;
  cout<<"ALLOCATING COUNTERS"<<endl;
  initialise_l0_sampler_counters(num_samplers,j,num_cols,num_rows,phi_s,iota_s,tau_s);
  cout<<"DONE"<<endl;

  // generate unique_hash_map for each sampler
  cout<<"GENERATING UNIQUE HASH MAPS"<<endl;
  map<vertex,int>* unique_hash_maps=new map<vertex,int>[num_samplers];
  for (int i=0; i<num_samplers; i++) {
    if (i%100==0) cout<<i<<",";
    map<vertex,int> hash_map;
    generate_random_hash(num_vertices,pow(num_vertices,3),hash_map);
    unique_hash_maps[i]=hash_map;
  }
  cout<<"DONE"<<endl;

  // generate hashs for s-sparse recovery
  hash_params*** ps_s=(hash_params***) malloc(num_samplers*sizeof(hash_params**));
  for (int i=0; i<num_samplers; i++) {
    ps_s[i]=initalise_2d_hash_params_array(j,num_rows); // each col is for one s-sparse recovery
    for (int k=0; k<j; k++) ps_s[i][k]=choose_hash_functions(num_cols,num_rows);
  }

  // hash limits for each value of j
  int* hash_lims=new int[j];
  int n3=pow(num_vertices,3);
  for (int i=1; i<=j; i++) hash_lims[i]=n3/pow(2,i);

  string line; edge e; int v;
  int sparsity_estimate=0; // r in survery paper algorithm 2
  ifstream edge_stream(file_path);

  cout<<"STREAM STARTING"<<endl;
  int edge_counter=0;
  while (getline(edge_stream,line)) {
    parse_edge(line,e);
    edge_counter+=1;
    if (edge_counter%10000==0) cout<<edge_counter<<",";
    v=identify_endpoint(e,target); // determine if edge is connected to target, if so what is the other endpoint

    if (v!=-1) { // edge is connected to target vertex
      sparsity_estimate+=e.value; // increment/decrement depending upon insertion or deletion edge
      // update each sampler independently
      for (int i=0; i<num_samplers; i++) {
        int h=unique_hash_maps[i][v]; // TODO look at this
        for (int k=0; k<j; k++) if (h<=hash_lims[k]) update_s_sparse(v,e.value,num_rows,ps_s[i][k],phi_s[i][k],iota_s[i][k],tau_s[i][k]);
    }}

  }
  cout<<"DONE"<<endl;

  //  j_sample is shared
  int j_sample=log2(sparsity_estimate)-1; // -1 since 0 indexed

  set<vertex> sampled_neighbourhood, sampled_vertices;
  map<vertex,int> edge_count; int succ_count=0;

  for (int i=0; i<num_samplers; i++) {
    // extract from each sampler
    sampled_neighbourhood=recover_neighbourhood(num_cols,num_rows,phi_s[i][j_sample],iota_s[i][j_sample],tau_s[i][j_sample]);
    vertex sampled_vertex=recover_vertex(sampled_neighbourhood,s,unique_hash_maps[i]);

    if (sampled_vertex!=-1) {
      sampled_vertices.insert(sampled_vertex);
      succ_count+=1;
      if (edge_count.count(sampled_vertex)) edge_count[sampled_vertex]+=1; // update number of sample occurences
      else edge_count[sampled_vertex]=1;
    }
  }

  // report results
  cout<<endl<<"Recovered Vertices:";
  for (set<vertex>::iterator it=sampled_vertices.begin(); it!=sampled_vertices.end(); it++) cout<<*it<<",";
  cout<<endl;

  cout<<"SUCCESSES:"<<succ_count<<"/"<<num_samplers<<endl;
  cout<<"Unique:"<<sampled_vertices.size()<<endl;
  write_to_file("uniform_sample_test.csv",edge_count);

  // free space
  free_l0_sampler_counters(num_samplers,j,num_cols,num_rows,phi_s,iota_s,tau_s);
  free(phi_s); free(phi_s); free(tau_s); free(unique_hash_maps);

  for (int i=0; i<num_samplers; i++) free_2d_hash_params_array(ps_s[i],num_cols,num_rows);
  free(ps_s);
}

/*-----------*
 * L0 Stream *
 *-----------*/

// allocates space for the counters used for 1-sparse recovery
// long**** = [sampler][s-sparse][s-sparse col][s-sparse row]
void initialise_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long****& phi_s, long****& iota_s, long****& tau_s) {
  // prepare pointers
  phi_s =(long****) malloc(num_samplers*sizeof(long***));
  iota_s=(long****) malloc(num_samplers*sizeof(long***));
  tau_s =(long****) malloc(num_samplers*sizeof(long***));

  for (int i=0; i<num_samplers; i++) {
    if (i%100==0) cout<<i<<",";
    phi_s[i] =initalise_zero_3d_array(j,num_cols,num_rows); // sum of weights (sum ai)
    iota_s[i]=initalise_zero_3d_array(j,num_cols,num_rows); // weighted sum of weights (sum ai*i)
    tau_s[i] =initalise_zero_3d_array(j,num_cols,num_rows); // squared weighted sum of weights (sum ai*(i**2))
  }
}

// free space used for counters
void free_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long**** phi_s, long**** iota_s, long**** tau_s) {
  for (int i=0; i<num_samplers; i++) {
    free_3d_long_array(phi_s[i],j,num_cols,num_rows);
    free_3d_long_array(iota_s[i],j,num_cols,num_rows);
    free_3d_long_array(tau_s[i],j,num_cols,num_rows);
  }
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
void generate_random_hash(int n, int m, map<vertex,int>& hash) {
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
  free(arr); // NEW
}

// free space of 2d array of hash parameters
void free_2d_hash_params_array(hash_params** arr, int num_cols, int num_rows) {
  for (int c=0; c<num_cols; c++) {
    hash_params* ptr=arr[c];
    free(ptr);
  }
  free(arr);
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
