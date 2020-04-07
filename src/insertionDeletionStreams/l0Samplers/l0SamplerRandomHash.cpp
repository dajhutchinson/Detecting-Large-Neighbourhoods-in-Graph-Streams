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
  int s=1/delta; // sparsity to recover at
  int j=log2(num_vertices); // number of s-sparse recoveries to run
  cout<<"s="<<s<<",j="<<j<<endl;
  // hash function for whether to pass update to an s-spare recovery

  // check if edges sampled ~ uniformly
  map<vertex,int> edge_count;

  int succ_count=0;
  int num_runs=100000;
  for (int run=0; run<num_runs; run++) {
    map<int,int> unique_hash_map;
    generate_random_hash(num_vertices,pow(num_vertices,3),unique_hash_map);
    //cout<<run<<","<<"HASH FUNCTION CHOSEN"<<endl;

    /** Prepare j s-sparse recoveries **/
    int num_cols=2*s;
    int num_rows=log(s/gamma);
    //cout<<run<<","<<"NUMBER COLS & ROWS DECIDED ("<<num_cols<<"x"<<num_rows<<")"<<endl;

    // arrays for 1-sparse recovery
    long*** phi_s =initalise_zero_3d_array(j,num_cols,num_rows); // sum of weights (sum ai)
    long*** iota_s=initalise_zero_3d_array(j,num_cols,num_rows); // weighted sum of weights (sum ai*i)
    long*** tau_s =initalise_zero_3d_array(j,num_cols,num_rows); // squared weighted sum of weights (sum ai*(i**2))
    //cout<<run<<","<<"3D ARRAYS PREPARED"<<endl;

    //cout<<run<<","<<"NUM 1-SPARSE SAMPLERS "<<j*num_cols*num_rows<<endl;

    // choose hash function for each row
    hash_params** ps_s=initalise_2d_hash_params_array(j,num_rows); // each col is for one s-sparse recovery
    for (int i=0; i<j; i++) {
      for (int j=0; j<num_rows; j++) {
        ps_s[i][j]=generate_hash(num_cols);
      }
    }
    //cout<<run<<","<<"HASH FUNCTIONS CHOSEN"<<endl;

    // hash limits for each value of j
    int* hash_lims=new int[j];
    int n3=pow(num_vertices,3);
    for (int i=1; i<=j; i++) hash_lims[i]=n3/pow(2,i);
    //cout<<run<<","<<"HASH LIMITS FOUND"<<endl;

    // Process stream
    string line; edge e; int v;
    ifstream edge_stream(file_path);

    //cout<<run<<","<<"BEGINNING PROCESSING STREAM"<<endl;
    int sparsity_estimate=0; // r in survery paper algorithm 2
    while (getline(edge_stream,line)) {
      parse_edge(line,e);
      v=identify_endpoint(e,target); // determine if edge is connected to target, if so what is the other endpoint

      if (v!=-1) { // edge is connected to target vertex
        sparsity_estimate+=e.value; // increment/decrement depending upon insertion or deletion edge
        for (int i=0; i<j; i++) { // for each s-sparse recovery
          //int h=hash_function(v,ps); // calculate hash_value
          int h=unique_hash_map[v];
          if (h<=hash_lims[i]) { // update ith s-sparse recovery
            //cout<<run<<",PREPARING TO UPDATE"<<i<<"/"<<j<<endl;
            for (int r=0; r<num_rows; r++) { // decide which sampler in each row to update
              int c=hash_function(v,ps_s[i][r]); // col to update // TODO change this
              update_1_sparse_counters(v,e.value,r,c,phi_s[i],iota_s[i],tau_s[i]);
            }
            //cout<<run<<",UPDATED"<<i<<"/"<<j<<endl;
          }
        }
      }
    }
    //cout<<run<<","<<"STEAM PROCESSED"<<endl;

    // Gather sample from j_sample^th s-sparse recovery
    set<vertex> sampled_neighbourhood;

    //cout<<run<<","<<"SPARSITY ESTIMATE "<<sparsity_estimate<<endl;
    int j_sample=log2(sparsity_estimate)-1; // -1 since 0 indexed
    //cout<<"j_sample "<<j_sample<<" ("<<num_cols*num_rows<<")"<<endl;

    for (int c=0; c<num_cols; c++) {
      for (int r=0; r<num_rows; r++) {
        if (verify_1_sparse(phi_s[j_sample][c][r],iota_s[j_sample][c][r],tau_s[j_sample][c][r])) {
          sampled_neighbourhood.insert(iota_s[j_sample][c][r]/phi_s[j_sample][c][r]);
    }}}

    // NOTE the official algorithm is to only check the j_sample^th estimator but that has a low success rate
    // so instead I am do it for the first to succeed
    /*int i=0;
    while (sampled_neighbourhood.size()==0 && i<j) {
      for (int c=0; c<num_cols; c++) {
        for (int r=0; r<num_rows; r++) {
          if (verify_1_sparse(phi_s[i][c][r],iota_s[i][c][r],tau_s[i][c][r])) {
            sampled_neighbourhood.insert(iota_s[i][c][r]/phi_s[i][c][r]);
      }}}
      i++;
    }*/

    //cout<<run<<","<<"NEIGHBOURHOOD SIZE "<<sampled_neighbourhood.size()<<endl;


    if (sampled_neighbourhood.size()>s || sampled_neighbourhood.size()==0) {
      cout<<run<<","<<"FAIL"<<endl;
    } else {
      int min_hash=pow(num_vertices,3), min_val=-1;
      for (set<vertex>::iterator it=sampled_neighbourhood.begin(); it!=sampled_neighbourhood.end(); it++) {
        int h_i=unique_hash_map[*it];
        if (h_i<min_hash) {
          min_hash=h_i;
          min_val=*it;
        }
      }
      cout<<run<<","<<"SUCCESS ("<<min_val<<")"<<endl;

      if (edge_count.count(min_val)) edge_count[min_val]+=1; // update number of sample occurences
      else edge_count[min_val]=1;

      succ_count+=1;
    }
    free_3d_long_array(phi_s,j,num_cols,num_rows);
    free_3d_long_array(iota_s,j,num_cols,num_rows);
    free_3d_long_array(tau_s,j,num_cols,num_rows);
    free_2d_hash_params_array(ps_s,j,num_rows);
  }
  cout<<endl<<succ_count<<"/"<<num_runs<<endl;

  write_to_file("uniform_sample_test.csv", edge_count);
}

/*-------------------*
 * 1 SPARSE RECOVERY *
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
