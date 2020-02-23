/*-----------------*
 * One pass c-approximation Streaming Algorithm for Neighbourhood Detection for Insertion-Only Graph Streams
 * TODO
 *  c has to be quite large??
 *-----------------*/

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <set>

using namespace std;
int P=1073741789; // >2^30

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

using vertex = string; // typemap vertex
using time_point=chrono::high_resolution_clock::time_point;

struct edge { // undirected edge
   vertex fst;
   vertex snd;
   bool insertion;
};

struct hash_params { // parameters for hash function
  unsigned long a;
  unsigned long b;
  unsigned long m;
};

/*-----------*
 * SIGNATURES *
 *------------*/

//  main algorithm
int vertex_sampling(string edge_file_path, string vertex_file_path, int num_vertices, int c, int d, set<vertex> &neighbourhood, vertex& root);
void sample_vertices(string vertex_file_path, int num_vertices, int sample_size, set<vertex>& vertex_sample, set<vertex>& vertex_not_in_sample);
void index_vertices(set<vertex> vertex_sample, set<vertex> vertex_not_in_sample, map<vertex,int>& vertex_to_index, map<int,vertex>& index_to_vertex);
map<vertex,int*> initialise_edge_vectors(set<vertex> vertex_sample, int num_vertices);
void construct_edge_vectors(string edge_file_path, set<vertex> vertex_sample, int num_vertices, map<vertex,int> vertex_to_index, map<vertex,int*>& edge_vectors);
int run_l0_samplers(set<vertex> vertex_sample, map<vertex,int*> edge_vectors, int num_vectors, map<int,vertex> index_to_vertex, int num_samplers, int neighbourhood_size, set<vertex>& neighbourhood, vertex&root);

// l0
int l0_sample(int* aj, int size);
int* sample(int* a, int size, int j, hash_params ps);
int recover(int* aj, int size);

// hashing
hash_params generate_hash(int m);
int hash_function(int key, hash_params ps);

// utility
void parse_edge(string str, edge& e);
int* copy_arr(int* arr, int size);

/*------*
 * BODY *
 *------*/

int main() {
  string edge_file_path="data/gplus_deletion.edges";
  string vertex_file_path="data/gplus.vertices";

  int n=12417; // number of vertices
  int c=200;
  int d=4998;
  /*
  string edge_file_path="data/gplus_large_deletion.edges", vertex_file_path="data/gplus.vertices";

  int n=102100, c=200, d=10000;
  */

  vertex root; set<vertex> neighbourhood;

  time_point before=chrono::high_resolution_clock::now(); // time before execution
  int result=vertex_sampling(edge_file_path,vertex_file_path,n,c,d,neighbourhood,root);
  time_point after=chrono::high_resolution_clock::now(); // time before execution

  if (result==0) { // success
    cout<<"SOLVED"<<endl<<"root="<<root<<endl<<"size neighbourhood="<<neighbourhood.size()<<endl<<"Neighbourhood={";
    for (set<vertex>::iterator it=neighbourhood.begin(); it!=neighbourhood.end(); it++) cout<<*it<<",";
    cout<<"}"<<endl;
  } else { // failure
    cout<<"FAIL"<<endl;
  }

  auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count(); // time passed
  cout<<"Duration="<<duration/1000<<"ms"<<endl;

  return 0;
}

/*----------------*
 * MAIN ALGORITHM *
 *----------------*/

int vertex_sampling(string edge_file_path, string vertex_file_path, int num_vertices, int c, int d, set<vertex> &neighbourhood, vertex& root) {
  int x=max((int)num_vertices/c,(int)sqrt(num_vertices));
  int sample_size=(int)(10*log(num_vertices)*x);

  cout<<"n="<<num_vertices<<", sample_size="<<sample_size<<endl<<"proportion="<<(float)sample_size/(float)num_vertices<<endl;

  set<vertex> vertex_sample, vertex_not_in_sample;
  map<vertex,int> vertex_to_index; map<int,vertex> index_to_vertex;

  sample_vertices(vertex_file_path,num_vertices,sample_size,vertex_sample,vertex_not_in_sample);
  cout<<vertex_sample.size()<<","<<vertex_not_in_sample.size()<<endl;

  index_vertices(vertex_sample,vertex_not_in_sample,vertex_to_index,index_to_vertex); // TODO stop requirement for double here
  map<vertex,int*> edge_vectors=initialise_edge_vectors(vertex_sample,num_vertices);
  construct_edge_vectors(edge_file_path,vertex_sample,num_vertices,vertex_to_index,edge_vectors);

  int neighbourhood_size=d/c;
  int num_samplers=10*(neighbourhood_size)*log(num_vertices);
  cout<<"d/c="<<neighbourhood_size<<endl<<"# samplers="<<num_samplers<<endl;
  int result=run_l0_samplers(vertex_sample,edge_vectors,num_vertices,index_to_vertex,num_samplers,neighbourhood_size,neighbourhood,root);
  return result;
}

// generate uniform random sample of vertices
void sample_vertices(string vertex_file_path, int num_vertices, int sample_size, set<vertex>& vertex_sample, set<vertex>& vertex_not_in_sample) {

  // choose which vertices to sample
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
  uniform_int_distribution<int> uniform_d(0,num_vertices-1);

  set<int> sample_indicies;
  while(sample_indicies.size()<sample_size) sample_indicies.insert(uniform_d(generator));

  // sample vertices from file
  string line; int i=0;
  set<int>::iterator it=sample_indicies.begin();
  ifstream stream(vertex_file_path);
  while (getline(stream,line)) {
    if(*it==i) { // sample vertex
      it++; // next index to find
      vertex_sample.insert(line);
    } else {
      vertex_not_in_sample.insert(line);
    }
    i++;
  }

}

// map vertices to an index used for edge vectors
void index_vertices(set<vertex> vertex_sample, set<vertex> vertex_not_in_sample, map<vertex,int>& vertex_to_index, map<int,vertex>& index_to_vertex) {
  int i=0;
  for (set<vertex>::iterator it=vertex_sample.begin();it!=vertex_sample.end(); it++) {
    vertex_to_index[*it]=i;
    index_to_vertex[i]=*it;
    i++;
  }

  for (set<vertex>::iterator it=vertex_not_in_sample.begin();it!=vertex_not_in_sample.end(); it++) {
    vertex_to_index[*it]=i;
    index_to_vertex[i]=*it;
    i++;
  }
}

// initalise edge vector for each vertex in sample
map<vertex,int*> initialise_edge_vectors(set<vertex> vertex_sample, int num_vertices) {
  map<vertex,int*> edge_vectors;
  for (set<vertex>::iterator it=vertex_sample.begin();it!=vertex_sample.end(); it++) {
    int* zeroes=new int[num_vertices]; // initalise zero vector
    for (int i=0; i<num_vertices; i++) zeroes[i]=0;
    edge_vectors[*it]=zeroes;
  }
  return edge_vectors;
}

// construct edge vector for each sample vertex
void construct_edge_vectors(string edge_file_path, set<vertex> vertex_sample, int num_vertices, map<vertex,int> vertex_to_index, map<vertex,int*>& edge_vectors) {
  ifstream stream(edge_file_path);
  string line; edge e;
  while (getline(stream,line)) {
    parse_edge(line,e);
    int fst_index=vertex_to_index[e.fst];
    int snd_index=vertex_to_index[e.snd];

    if (vertex_sample.find(e.fst)!=vertex_sample.end()) {
      edge_vectors[e.fst][snd_index]=e.insertion;
    }

    if (vertex_sample.find(e.snd)!=vertex_sample.end()) {
      edge_vectors[e.snd][fst_index]=e.insertion;
    }
  }
}

int run_l0_samplers(set<vertex> vertex_sample, map<vertex,int*> edge_vectors, int num_vectors, map<int,vertex> index_to_vertex, int num_samplers, int neighbourhood_size, set<vertex>& neighbourhood, vertex&root) {
  // run l0 sampler on each v in sample
  for (set<vertex>::iterator it=vertex_sample.begin(); it!=vertex_sample.end(); it++) {
    int* edge_vector=edge_vectors[*it];
    root=*it;
    neighbourhood.clear();
    for (int i=0; i<num_samplers; i++) {
      int sampled_index=l0_sample(edge_vector,num_vectors);
      if (sampled_index!=-1) neighbourhood.insert(index_to_vertex[sampled_index]);
      if (neighbourhood.size()>=neighbourhood_size) return 0; // success
    }
  }
  return -1;
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

/*-------------*
 * L0 SAMPLING *
 *-------------*/

// perform l0-sampling on a
int l0_sample(int* a, int size) {
  int j=1;
  int* a_j=copy_arr(a,size);

  while (true) {
    //cout<<"j="<<j<<endl;
    hash_params ps=generate_hash(pow(size,3)); // m=n^3 where n is |x|

    a_j=sample(a_j,size,j,ps); // sample from vector
    //for (int i=0; i<20; i++) cout << a_j[i]<<",";
    //cout<<endl;

    int r=recover(a_j,size); // try to recover an index
    if (r==-2) { // fail
      return -1;
    } else if (r!=-1) { // success
      return r;
    }

    j+=1; // try again
  }
}

// sample from vector a with condition j
int* sample(int* a, int size, int j, hash_params ps) {
  int* a_j= new int[size];
  for (int i=0; i<size; i+=1) {
    int h=hash_function(i,ps);
    if (h<=ps.m/pow(2,j)) a_j[i]=a[i];
    else a_j[i]=0;
  }
  return a_j;
}

// recover a non-zero index of a if a is 1-sparse
int recover(int* aj, int size) {
  int w1=0, w2=0; // w1=sum aj_i, w_2=sum i*aj_i

  for (int i=0; i<size; i++) {
      w1+=aj[i];
      w2+=i*aj[i];
  }

  if (w1==1) return w2; // 1-sparse
  if (w1==0) return -2; // zero string
  return -1; // not 1-sparse
}

/*-----------*
 * UTILITIES *
 *-----------*/

// parse insertion-deletion edge from string
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

// create copy of integer array
int* copy_arr(int* arr, int size) {
  int* new_arr=new int[size];
  for (int i=0; i<size; i++) new_arr[i]=arr[i];
  return new_arr;
}
