/*
 * One pass c-approximation Streaming Algorithm for Neighbourhood Detection for Insertion-Deletion Graph Streams (Using Edge Sampling)
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

uint64_t BYTES; // space used atm
uint64_t L0_HASH_BYTES; // space used atm
uint64_t GENERATING_L0_HASH_TIME; // space used atm
int P=1073741789; // >2^30

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

using vertex = int; // typemap vertex
using time_point=chrono::high_resolution_clock::time_point;

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

void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string edge_file_path, string vertex_file_path, string out_file);

// Main algorithm
void single_pass_insertion_deletion_stream(int c, int d, int num_vertices, string edge_file_path, string vertex_file_path, set<vertex>& neighbourhood, vertex& root);

// L0 Stream
void initialise_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long****& phi_s, long****& iota_s);
void free_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long**** phi_s, long**** iota_s);

// s-sparse
hash_params* choose_hash_functions(int num_cols, int num_rows);
void update_s_sparse(uint64_t endpoint, int edge_value, int num_rows, hash_params* ps_s, long** phi_s, long** iota_s);
set<uint64_t> recover_neighbourhood(int num_cols, int num_rows, long** phi_s, long** iota_s);
uint64_t recover_id(set<uint64_t> neighbourhood, int sparsity, uint64_t* hash_map);

// 1-sparse
bool verify_1_sparse(int phi,int iota);
void update_1_sparse_counters(uint64_t index,int delta,int row,int col,long** phi_s,long** iota_s);

// Hashing
hash_params generate_hash(int m);
int hash_function(int key, hash_params ps);
// n = num keys, m=possible hash values, hash=map from key to hash value
uint64_t* generate_random_hash(int n, uint64_t m);

// Utility
uint64_t edge_id(string str, int num_vertices);
edge unparse_edge_id(uint64_t id, int num_vertices, int edge_value);
bool edge_type(string str);
void parse_vertex(string str, vertex& v);
int identify_endpoint(edge e,vertex target);
long*** initalise_zero_3d_array(int depth, int num_cols, int num_rows);
hash_params** initalise_2d_hash_params_array(int num_cols, int num_rows);
void free_3d_long_array(long*** arr, int depth, int num_cols, int num_rows);
void free_2d_hash_params_array(hash_params** arr, int num_cols, int num_rows);
double variance(vector<uint64_t> vals);
uint64_t mean(vector<uint64_t> vals);

int main() {

  string edge_file_path, vertex_file_path, out_file;
  int num_vertices, d, reps;

  // details of graph to perform on
  //edge_file_path="../../data/gplus_deletion.edges"; vertex_file_path="../../data/gplus_deletion.vertices"; out_file="gplus_results.csv"; num_vertices=12417; d=4998; reps=10;
  //execute_test(2,20,1,reps,d,num_vertices,edge_file_path,vertex_file_path,out_file);

  // details of graph to perform on
  //edge_file_path="../../data/facebook_deletion.edges"; vertex_file_path="../../data/facebook_deletion.vertices"; num_vertices=747; d=267; reps=10;
  edge_file_path="../../data/facebook_small_deletion.edges"; vertex_file_path="../../data/facebook_deletion.vertices"; num_vertices=52; d=33; reps=10;
  out_file="edge_sampled_better_id.csv";
  execute_test(2,d,1,reps,d,num_vertices,edge_file_path,vertex_file_path,out_file);

  //set<vertex> neighbourhood; vertex root; // variables for returned values
  //int c=2;
  //single_pass_insertion_deletion_stream(c,d,num_vertices,edge_file_path,vertex_file_path,neighbourhood,root);
}

void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string edge_file_path, string vertex_file_path, string out_file) {
  ofstream outfile(out_file);
  outfile<<"name,"<<vertex_file_path<<endl<<"n,"<<n<<endl<<"d,"<<d<<endl<<"repetitions,"<<reps<<endl<<"delta,0.2"<<endl<<"gamma,0.3"<<endl<<"sample size,((num_vertices*d)/((double)c))*(1/((double)x)+1/((double)c))*2*log(num_vertices)"<<endl<<endl; // test details
  outfile<<"c,time (microseconds), generating l0 hash time (microseconds) ,mean max space (bytes), l0 hash space (bytes), variance time, variance hash time, variance max space, variance hash space,successes"<<endl; // headers
  set<vertex> neighbourhood; vertex root; // variables for returned values
  vector<uint64_t> times, total_space, hash_times, hash_space; // results of each run of c
  int successes;
  //for (int c=c_min;c<=c_max;c+=c_step) {
  for (int c=c_max;c>=c_min;c-=c_step) {
    successes=0;
    times.clear(); total_space.clear(); hash_times.clear(); hash_space.clear();
    for (int i=0;i<reps;i++) {
      cout<<"("<<i<<"/"<<reps<<") "<<c<<"/"<<c_max<<endl; // output to terminal

      // reset values
      BYTES=0; L0_HASH_BYTES=0; GENERATING_L0_HASH_TIME=0;
      neighbourhood.clear(); vertex* p=&root; p=nullptr;

      time_point before=chrono::high_resolution_clock::now(); // time before execution
      single_pass_insertion_deletion_stream(c,d,n,edge_file_path,vertex_file_path,neighbourhood,root);
      time_point after=chrono::high_resolution_clock::now(); // time after execution

      cout<<endl<<root<<endl;
      cout<<neighbourhood.size()<<"("<<d/c<<")"<<endl<<endl;

      if (neighbourhood.size()!=0) successes+=1;

      auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count(); // time passed
      cout<<"DURATION "<<duration/1000000<<"s"<<endl<<"SPACE "<<BYTES/pow(1024,2)<<"MBs"<<endl;
      cout<<"HASH DURATION "<<GENERATING_L0_HASH_TIME/1000000<<"s"<<endl<<"HASH SPACE "<<L0_HASH_BYTES/pow(1024,2)<<"MBs"<<endl<<endl;
      cout<<"STORING VALUES";
      times.push_back(duration); total_space.push_back(BYTES);
      cout<<" *";
      hash_times.push_back(GENERATING_L0_HASH_TIME); hash_space.push_back(L0_HASH_BYTES);
      cout<<"\rVALUES STORED               "<<endl;
    }
    cout<<"CALCULATING MEAN";
    uint64_t mean_duration   =mean(times);
    uint64_t mean_total_space=mean(total_space);
    uint64_t mean_hash_space =mean(hash_space);
    uint64_t mean_hash_time  =mean(hash_times);

    double variance_duration   =variance(times);
    double variance_total_space=variance(total_space);
    double variance_hash_space =variance(hash_space);
    double variance_hash_time  =variance(hash_times);
    cout<<"\rCALCULATED VARAIANCE"<<endl<<endl;
    outfile<<c<<","<<mean_duration<<","<<mean_hash_time<<","<<mean_total_space<<","<<mean_hash_space<<","<<variance_duration<<","<<variance_hash_time<<","<<variance_total_space<<","<<variance_hash_space<<","<<successes<<endl; // write values to file
  }
  outfile.close();
}

/*----------------*
 * MAIN ALGORITHM *
 *----------------*/

void single_pass_insertion_deletion_stream(int c, int d, int num_vertices, string edge_file_path, string vertex_file_path, set<vertex>& neighbourhood, vertex& root) {

  // L0 sampling parameters
  double delta=0.2, gamma=0.3; // success_rate ~= P(L0 sampler returning a vertex given delta & gamma)

  // calculate model feature
  int x=(num_vertices/(double)c>sqrt(num_vertices)) ? num_vertices/(double)c : sqrt(num_vertices);
  int total_samplers=((num_vertices*d)/((double)c))*(1/((double)x)+1/((double)c))*2*log(num_vertices);
  int possible_edges=ceil(.5*num_vertices*(num_vertices-1));
  cout<<"x:"<<x<<endl<<"Total Samplers:"<<total_samplers<<endl<<"d/c:"<<d/c<<endl;

  // prepare samplers
  // sampler parameters
  int s=1/delta; // sparsity to recover at
  int j=log2(possible_edges); // number of s-sparse recoveries to run
  int num_cols=2*s;
  int num_rows=log(s/gamma);

  cout<<"Sparsity of s-sparse:"<<s<<endl<<"# s-sparse per L0:"<<j<<endl<<"# cols per s-sparse:"<<num_cols<<endl<<"# rows per s-sparse:"<<num_rows<<endl;

  // allocate space for counters
  long**** phi_s; long**** iota_s;
  cout<<"ALLOCATING COUNTERS"<<endl;
  initialise_l0_sampler_counters(total_samplers,j,num_cols,num_rows,phi_s,iota_s);
  cout<<"\rDONE                                             "<<endl; // spaces to "clear" line

  // generate unique_hash_map for each sampler
  time_point before=chrono::high_resolution_clock::now(); // time before execution
  uint64_t** unique_hash_maps=(uint64_t**) malloc(total_samplers * sizeof(uint64_t*));
  uint64_t n3=pow(possible_edges,3);
  cout<<"possible_edges="<<possible_edges<<", n3="<<n3<<endl;
  for (int i=0; i<total_samplers; i++) {
    if (i%10==0) cout<<"\r"<<i;
    uint64_t* hash_map=generate_random_hash(possible_edges,n3);
    unique_hash_maps[i]=hash_map;
  }
  time_point after=chrono::high_resolution_clock::now(); // time before execution
  GENERATING_L0_HASH_TIME=chrono::duration_cast<chrono::microseconds>(after-before).count();

  // generate hashs for s-sparse recovery
  hash_params*** ps_s=(hash_params***) malloc(total_samplers*sizeof(hash_params**));
  for (int i=0; i<total_samplers; i++) {
    ps_s[i]=initalise_2d_hash_params_array(j,num_rows); // each col is for one s-sparse recovery
    for (int k=0; k<j; k++) ps_s[i][k]=choose_hash_functions(num_cols,num_rows);
  }

  // hash limits for each value of j
  uint64_t* hash_lims=new uint64_t[j];
  for (int i=1; i<=j; i++) hash_lims[i]=n3/pow(2,i);

  int sparsity_estimate=0;

  ifstream edge_stream(edge_file_path);

  string line; edge e;
  int edge_counter=0;
  int edge_value;
  uint64_t id;
  while (getline(edge_stream,line)) {
    edge_counter+=1;
    if (edge_counter%1000==0) cout<<"\r"<<edge_counter;

    id=edge_id(line,num_vertices);
    edge_value=edge_type(line) ? 1 : -1;
    sparsity_estimate+=edge_value;
    //cout<<line<<" "<<id<<" ("<<edge_value<<")"<<endl;

    for (int i=0; i<total_samplers; i++) { // pass to samplers
      uint64_t h=unique_hash_maps[i][id];
      for (int k=0; k<j; k++) if (h<=hash_lims[k]) {
        update_s_sparse(id,edge_value,num_rows,ps_s[i][k],phi_s[i][k],iota_s[i][k]);
      }
    }

  }
  cout<<"\rDONE                                             "<<endl; // spaces to "clear" line

  // recover from each sampler
  set<uint64_t> sampled_neighbourhood, sampled_ids;
  int successes=0;
  int j_sample=log2(sparsity_estimate)-1;
  cout<<"j_sample="<<j_sample<<endl;

  for (int i=0; i<total_samplers; i++) {
    cout<<"\r"<<i<<"/"<<total_samplers<<"   "<<phi_s[i][j_sample][0][0]<<","<<iota_s[i][j_sample][0][0];
    sampled_neighbourhood=recover_neighbourhood(num_cols,num_rows,phi_s[i][j_sample],iota_s[i][j_sample]);
    cout<<"*";
    uint64_t sampled_id=recover_id(sampled_neighbourhood,s,unique_hash_maps[i]);
    cout<<"*";
    if (sampled_id!=-1) {
      successes+=1;
      sampled_ids.insert(sampled_id);
    }
    cout<<"*";
  }
  cout<<"\rDONE                 "<<endl;
  cout<<sampled_ids.size()<<"/"<<successes<<"/"<<total_samplers<<endl;

  // check for sufficient neighbourhood
  map<vertex, vector<vertex>> neighbourhoods;
  for (set<uint64_t>::iterator it=sampled_ids.begin(); it!=sampled_ids.end(); it++) {

    edge e=unparse_edge_id(*it,num_vertices,1);

    if (neighbourhoods.find(e.fst)==neighbourhoods.end()) {
      vector<vertex> n;
      n.push_back(e.snd);
      neighbourhoods[e.fst]=n;
    } else neighbourhoods[e.fst].push_back(e.snd);

    if (neighbourhoods[e.fst].size()>=(int)d/c) { // sufficient neighbourhood found
      root=e.fst;
      for (vector<vertex>::iterator it=neighbourhoods[e.fst].begin(); it!=neighbourhoods[e.fst].end(); it++) neighbourhood.insert(*it);
      cout<<"NEIGHBOURHOOD OF "<<root<<" ("<<neighbourhood.size()<<")={";
      for (set<vertex>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) cout<<*it<<",";
      cout<<"\b}";
      return;
    }

    if (neighbourhoods.find(e.snd)==neighbourhoods.end()) {
      vector<vertex> n;
      n.push_back(e.fst);
      neighbourhoods[e.snd]=n;
    } else neighbourhoods[e.snd].push_back(e.snd);

    if (neighbourhoods[e.snd].size()>=(int)d/c) { // sufficient neighbourhood found
      root=e.snd;
      for (vector<vertex>::iterator it=neighbourhoods[e.snd].begin(); it!=neighbourhoods[e.snd].end(); it++) neighbourhood.insert(*it);
      cout<<"NEIGHBOURHOOD OF "<<root<<" ("<<neighbourhood.size()<<")={";
      for (set<vertex>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) cout<<*it<<",";
      cout<<"\b}";
      return;
    }

  }

  cout<<"\rFAILED to find neighbourhood";
  neighbourhood.clear();
  vertex* p=&root;
  p=nullptr;
  return;

}

/*-----------*
 * L0 Stream *
 *-----------*/

// allocates space for the counters used for 1-sparse recovery
// long**** = [sampler][s-sparse][s-sparse col][s-sparse row]
void initialise_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long****& phi_s, long****& iota_s) {
  // prepare pointers
  phi_s =(long****) malloc(num_samplers*sizeof(long***));
  iota_s=(long****) malloc(num_samplers*sizeof(long***));

  for (int i=0; i<num_samplers; i++) {
    if (i%100==0) cout<<"\r"<<i;
    try{
      phi_s[i] =initalise_zero_3d_array(j,num_cols,num_rows); // sum of weights (sum ai)
      iota_s[i]=initalise_zero_3d_array(j,num_cols,num_rows); // weighted sum of weights (sum ai*i)
    } catch (const exception &e) { cout<<"1,"<<e.what()<<endl; }
  }
}

// free space used for counters
void free_l0_sampler_counters(int num_samplers, int j, int num_cols, int num_rows, long**** phi_s, long**** iota_s) {
  for (int i=0; i<num_samplers; i++) {
    free_3d_long_array(phi_s[i],j,num_cols,num_rows);
    free_3d_long_array(iota_s[i],j,num_cols,num_rows);
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
void update_s_sparse(uint64_t endpoint, int edge_value, int num_rows, hash_params* ps_s, long** phi_s, long** iota_s) {
  for (int r=0; r<num_rows; r++) { // decide which sampler in each row to update
    int c=hash_function(endpoint,ps_s[r]); // col to update
    update_1_sparse_counters(endpoint,edge_value,r,c,phi_s,iota_s);
  }
}

// recover neighbourhood from s-sparse recovery counters
set<uint64_t> recover_neighbourhood(int num_cols, int num_rows, long** phi_s, long** iota_s) {
  set<uint64_t> neighbourhood;
  for (int c=0; c<num_cols; c++) {
    for (int r=0; r<num_rows; r++) {
      if (verify_1_sparse(phi_s[c][r],iota_s[c][r])) {
        neighbourhood.insert(iota_s[c][r]);
      }
  }}
  return neighbourhood;
}

// recover vertex from recovered neighbourhood
uint64_t recover_id(set<uint64_t> neighbourhood, int sparsity, uint64_t* hash_map) {
  if (neighbourhood.size()>sparsity || neighbourhood.size()==0) return -1; // s-sparse recovery failed
  else { // return vertex in neighbourhood with min hash value
    uint64_t min_hash=INT_MAX, min_val=-1;
    for (set<uint64_t>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) {
      uint64_t h_i=hash_map[*it];
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
void update_1_sparse_counters(uint64_t index,int delta,int row,int col,long** phi_s,long** iota_s) {
  try {
    //cout<<phi_s[col][row]<<","<<iota_s[col][row]<<","<<" ("<<delta<<","<<index<<") ";
    phi_s[col][row] +=delta;
    iota_s[col][row]+=delta*index;
    //cout<<phi_s[col][row]<<","<<iota_s[col][row]<<endl;
  } catch (exception exp) {
    cout<<"***********************************VERIFY 1 SPARSE COUNTERS FAILURE";
  }
}

// verify if array is 1_sparse
bool verify_1_sparse(int phi,int iota) {
  if (phi==1) return true;
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

// generate random hash with unique values for all keys
uint64_t* generate_random_hash(int n, uint64_t m) {
  vector<uint64_t> used; // record hash values which have been used
  uint64_t* hash_map=new uint64_t[n+1];
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
  uniform_int_distribution<uint64_t> distribution(0,m);
  for (int i=1; i<=n; i++) {
    // find unique hash_value
    uint64_t hash_value=distribution(generator);
    while (find(used.begin(), used.end(), hash_value)!=used.end()) hash_value=distribution(generator);

    hash_map[i]=hash_value;         // store in hash map
    used.push_back(hash_value); // record hash_value as used
  }
  return hash_map;
}


/*-----------*
 * UTILITIES *
 *-----------*/

// parse edge from string "v1 v2" or "(I/D) v1 v2"
uint64_t edge_id(string str, int num_vertices) {
  int spaces=count(str.begin(),str.end(),' ');

  if (spaces==2) { // insertion deletion edge
    str=str.substr(2,str.size());
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
    int f=stoi(fst);
    int s=stoi(snd);
    vertex max = (f>s) ? f-1 : s-1;
    vertex min = (f<s) ? f-1 : s-1;
    uint64_t id=(.5*(num_vertices)*(num_vertices-1))-(.5*(num_vertices-max)*(num_vertices-max-1))+min-max-1;

    return id;
  } catch (exception ex) {
    cout<<"EXCEPETION"<<endl;
    return -1;
  }
}

edge unparse_edge_id(uint64_t id, int num_vertices, int edge_value) {
  edge e;
  e.value=edge_value;

  // i = n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
  int min=num_vertices-2-floor((sqrt((-8*id)+(4*num_vertices*(num_vertices-1))-7)/2.0)-0.5);
  //j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
  int max=id+min+1-((1/2.0)*num_vertices*(num_vertices-1))+((1/2.0)*(num_vertices-min)*(num_vertices-min-1));

  e.fst=min+1;
  e.snd=max+1;

  return e;
}

// identify if insertion (true) or deletion (false)
bool edge_type(string str) {
  int spaces=count(str.begin(),str.end(),' ');
  if (spaces==2) { // insertion deletion edge
    if (str[0]=='D') return false;
    return true;
  }
  return true;
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
}

// free space of 2d array of hash parameters
void free_2d_hash_params_array(hash_params** arr, int num_cols, int num_rows) {
  for (int c=0; c<num_cols; c++) {
    hash_params* ptr=arr[c];
    free(ptr);
  }
}

// return variance of values in a vector
double variance(vector<uint64_t> vals) {
  if (vals.size()<=1) return 0;
  double var=0;
  double mean=accumulate(vals.begin(),vals.end(),0)/vals.size();

  for (vector<uint64_t>::iterator it=vals.begin(); it!=vals.end(); it++) var+=(*it-mean)*(*it-mean);
  var/=(vals.size()-1);

  return var;
}

uint64_t mean(vector<uint64_t> vals) {
  if (vals.size()==0) return 0;
  uint64_t sum=0;
  for (vector<uint64_t>::iterator it=vals.begin(); it!=vals.end(); it++) sum+=*it;

  return sum/vals.size();
}
