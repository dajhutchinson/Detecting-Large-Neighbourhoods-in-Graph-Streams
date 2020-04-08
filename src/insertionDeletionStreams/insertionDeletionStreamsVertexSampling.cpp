/*
 * One pass c-approximation Streaming Algorithm for Neighbourhood Detection for Insertion-Deletion Graph Streams
 *
 *  gamma=.20 & delta=.10 over good balance between success rate and space used (~92% success rate)
 *
 *  TODO
 *    Test vertex sample size against success rate
 *    Test number of samplers per sampled vertex (Can be done with maths right??)
 */

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace std;

int BYTES; // space used atm
int RESEVOIR_BYTES; // space used by resevoirs atm
int DEGREE_BYTES; // space used to store degrees of vertices
int MAX_BYTES; // max space used at any time
int MAX_RESEVOIR_BYTES; // max space used by resevoirs
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

// Main algorithm
set<vertex> generate_vertex_sample(string file_path, int num_vertices, int sample_size);

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
void parse_vertex(string str, vertex& v);
int identify_endpoint(edge e,vertex target);
long*** initalise_zero_3d_array(int depth, int num_cols, int num_rows);
hash_params** initalise_2d_hash_params_array(int num_cols, int num_rows);
void free_3d_long_array(long*** arr, int depth, int num_cols, int num_rows);
void free_2d_hash_params_array(hash_params** arr, int num_cols, int num_rows);

 /*------*
 * BODY *
 *------*/

int main() {
  // details of graph to perform on
  string edge_file_path="../../data/artifical/test_deletion_double.edges";
  string vertex_file_path="../../data/artifical/test_deletion_double.vertices";
  int num_vertices=1000;

  // approximation parameters
  int d=515, c=31;

  // L0 sampling parameters
  double delta=0.2, gamma=0.1, success_rate=0.93; // success_rate ~= P(L0 sampler returning a vertex given delta & gamma)

  // calculate model feature
  int reasonable_degree=17; // reasonable expectation for the degree of a node
  int vertex_sample_size=sqrt(num_vertices); // size of vertex sample
  int samplers_per_l0=(1/success_rate)*((log(c*reasonable_degree-d)-log(c*reasonable_degree))/(log(reasonable_degree-1)-log(reasonable_degree))); // TODO play with these
  int total_samplers=vertex_sample_size*samplers_per_l0;

  cout<<"Vertex sample size:"<<vertex_sample_size<<endl<<"Samplers per vertex:"<<samplers_per_l0<<endl<<"Total Samplers:"<<total_samplers<<endl;
  cout<<"d/c="<<d/c<<endl;

  // generate vertex_sample
  set<vertex> vertex_sample=generate_vertex_sample(vertex_file_path,num_vertices,vertex_sample_size);
  cout<<"Vertex sample={";
  for (set<vertex>::iterator it=vertex_sample.begin(); it!=vertex_sample.end(); it++) cout<<*it<<",";
  cout<<"\b}"<<endl;

  // prepare samplers
  // sampler parameters
  int s=1/delta; // sparsity to recover at
  int j=log2(num_vertices); // number of s-sparse recoveries to run
  int num_cols=2*s;
  int num_rows=log(s/gamma);
  cout<<"Sparsity of s-sparse:"<<s<<endl<<"# s-sparse per L0:"<<j<<endl<<"# cols per s-sparse:"<<num_cols<<endl<<"# rows per s-sparse:"<<num_rows<<endl;

  // allocate space for counters
  long**** phi_s; long**** iota_s; long**** tau_s;
  cout<<"ALLOCATING COUNTERS"<<endl;
  initialise_l0_sampler_counters(total_samplers,j,num_cols,num_rows,phi_s,iota_s,tau_s);
  cout<<"\rDONE"<<endl;

  // generate unique_hash_map for each sampler
  cout<<"GENERATING UNIQUE HASH MAPS"<<endl;
  map<vertex,int>* unique_hash_maps=new map<vertex,int>[total_samplers];
  for (int i=0; i<total_samplers; i++) {
    if (i%100==0) cout<<"\r"<<i;
    map<vertex,int> hash_map;
    generate_random_hash(num_vertices,pow(num_vertices,1.5),hash_map);
    unique_hash_maps[i]=hash_map;
  }
  cout<<"\rDONE"<<endl;

  // generate hashs for s-sparse recovery
  hash_params*** ps_s=(hash_params***) malloc(total_samplers*sizeof(hash_params**));
  for (int i=0; i<total_samplers; i++) {
    ps_s[i]=initalise_2d_hash_params_array(j,num_rows); // each col is for one s-sparse recovery
    for (int k=0; k<j; k++) ps_s[i][k]=choose_hash_functions(num_cols,num_rows);
  }

  // hash limits for each value of j
  int* hash_lims=new int[j];
  int n15=pow(num_vertices,1.5);
  for (int i=1; i<=j; i++) hash_lims[i]=n15/pow(2,i);

  int* sparsity_estimates=new int[total_samplers];
  for (int i=0;i<total_samplers;i++) sparsity_estimates[i]=0;
  ifstream edge_stream(edge_file_path);

  cout<<"STREAM STARTING"<<endl;
  string line; edge e; vertex v; vertex target;
  set<vertex>::iterator it; // track which vertex from sample is target for given sampler
  int edge_counter=0;
  while (getline(edge_stream,line)) {
    edge_counter+=1;
    if (edge_counter%10000==0) cout<<"\r"<<edge_counter;

    parse_edge(line,e);
    it=vertex_sample.begin(); // whole sample for edge
    target=*it;
    for (int i=0; i<total_samplers; i++) { // update every l0 sampler
      // check if target should be changed
      if (i>0 && i%samplers_per_l0==0) it++; target=*it;
      v=identify_endpoint(e,target); // determine if edge is connected to target, if so what is the other endpoint

      if (v!=-1) { // edge is connected to current target
        sparsity_estimates[i]+=e.value;
        int h=unique_hash_maps[i][v]; // update certain s-sparse recoveries
        for (int k=0; k<j; k++) if (h<=hash_lims[k]) update_s_sparse(v,e.value,num_rows,ps_s[i][k],phi_s[i][k],iota_s[i][k],tau_s[i][k]);
      }
    }

  }
  cout<<"\rDONE                                             "<<endl; // spaces to "clear" line

  // recover neighbourhood for each
  // return first neighbourhood of size > d/c
  set<vertex> sampled_vertices, sampled_neighbourhood;
  it=vertex_sample.begin(); // track which vertex from sample is target for given sampler
  for (int i=0; i<total_samplers; i++) {
    if (i>0 && i%samplers_per_l0==0) {it++; target=*it; sampled_vertices.clear();} // restart sampled_vertices since have run through all L0 samplers which were associated to a single sampled vertex

    int j_sample=log2(sparsity_estimates[i])-1; // -1 since 0 indexed

    sampled_neighbourhood=recover_neighbourhood(num_cols,num_rows,phi_s[i][j_sample],iota_s[i][j_sample],tau_s[i][j_sample]);
    vertex sampled_vertex=recover_vertex(sampled_neighbourhood,s,unique_hash_maps[i]);

    if (sampled_vertex!=-1) {
      sampled_vertices.insert(sampled_vertex);
      if (sampled_vertices.size()>=d/c) {
        cout<<"SUCCESSES ("<<sampled_vertices.size()<<")"<<endl;
        cout<<"NEIGHBOURHOOD for "<<target<<"={";
        for (set<vertex>::iterator it=sampled_vertices.begin(); it!=sampled_vertices.end(); it++) cout<<*it<<",";
        cout<<"\b}"<<endl;
        break;
      }
    }

    if (i==total_samplers-1) cout<<"FAILED to find neighbourhood";
  }

  // free space
  free_l0_sampler_counters(total_samplers,j,num_cols,num_rows,phi_s,iota_s,tau_s);
  free(phi_s); free(phi_s); free(tau_s); free(unique_hash_maps);

  for (int i=0; i<total_samplers; i++) free_2d_hash_params_array(ps_s[i],num_cols,num_rows);
  free(ps_s);
}

/*----------------*
 * MAIN ALGORITHM *
 *----------------*/

// Generate sample of vertices
set<vertex> generate_vertex_sample(string file_path, int num_vertices, int sample_size) {
  // chose indices of vertices
  set<vertex> indices;
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
  uniform_int_distribution<int> distribution(0,num_vertices-1);
  while (indices.size()<sample_size) { // set only contains unique items
    int i=distribution(generator); // generate new value
    indices.insert(i);
  }

  // run through stream pic out indices
  ifstream vertex_stream(file_path);
  set<vertex> sample;
  string line; vertex v; int counter=0;
  set<vertex>::iterator it=indices.begin();

  while (getline(vertex_stream,line)) {
    if (*it==counter) { // sample vertex
      parse_vertex(line,v);
      sample.insert(v);
      it++;
      if (it==indices.end()) break;
    }
    counter++; // increment line count
  }

  return sample;
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
    if (i%100==0) cout<<"\r"<<i;
    try{
      phi_s[i] =initalise_zero_3d_array(j,num_cols,num_rows); // sum of weights (sum ai)
      iota_s[i]=initalise_zero_3d_array(j,num_cols,num_rows); // weighted sum of weights (sum ai*i)
      tau_s[i] =initalise_zero_3d_array(j,num_cols,num_rows); // squared weighted sum of weights (sum ai*(i**2))
    } catch (const exception &e) { cout<<"1,"<<e.what()<<endl; }
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

// parse vertex from line in file
void parse_vertex(string str, vertex& v) {
  string vertex_name="";
  for (char& c:str) {
    if (c==',') break; // name done
    vertex_name+=c;
  }
  v=stoi(vertex_name);
}

// returns endpoint which is not target (or -1 if target not on edge)
int identify_endpoint(edge e,vertex target) {
  if (e.fst==target) return e.snd;
  else if (e.snd==target) return e.fst;
  else return -1;
}

// allocate space of 3d array of long intergers, all values set of 0
long*** initalise_zero_3d_array(int depth, int num_cols, int num_rows) {
  long*** arr;
  try { arr = (long***) malloc(depth * sizeof(long**)); // allocate depth
  } catch (const exception &e) { cout<<"2,"<<e.what()<<endl; }
  for (int i=0; i<depth; i++) {
    try { arr[i] = (long**) malloc(num_cols * sizeof(long*)); // allocate cols
    } catch (const exception &e) { cout<<"3,"<<e.what()<<endl; }

    for (int j=0; j<num_cols; j++) {
      try { arr[i][j]=(long*) malloc(num_rows * sizeof(long)); // allocate rows
      } catch (const exception &e) { cout<<"4,"<<e.what()<<endl; }
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
