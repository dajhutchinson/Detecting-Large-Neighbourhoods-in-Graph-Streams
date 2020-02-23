/*-----------------*
 * One pass c-approximation Streaming Algorithm for Neighbourhood Detection for Insertion-Only Graph Streams
 * NOTE - I have cheated by making bidirectional maps between vectors & indicies but only counting one of them
 * TODO
 *  c has to be quite large?? st (|A'|<=|A|)
 *  For all c where n/c < sqrt{n} surely space required is ~ equal??
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
int BYTES;
int MAX_BYTES;

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
map<vertex,bool*> initialise_edge_vectors(set<vertex> vertex_sample, int num_vertices);
void construct_edge_vectors(string edge_file_path, set<vertex> vertex_sample, int num_vertices, map<vertex,int> vertex_to_index, map<vertex,bool*>& edge_vectors);
int run_l0_samplers(set<vertex> vertex_sample, map<vertex,bool*> edge_vectors, int num_vectors, map<int,vertex> index_to_vertex, int num_samplers, int neighbourhood_size, set<vertex>& neighbourhood, vertex&root);

// l0
int l0_sample(bool* aj, int size);
void sample(bool* a, int size, int j, hash_params ps);
int recover(bool* aj, int size);

// hashing
hash_params generate_hash(int m);
int hash_function(int key, hash_params ps);

// utility
void parse_edge(string str, edge& e);
bool* copy_arr(bool* arr, int size);

/*------*
 * BODY *
 *------*/

int main() {
  string edge_file_path="data/gplus_deletion.edges",vertex_file_path="data/gplus.vertices";
  int n=12417,c=200,d=4998;
  /*
  string edge_file_path="data/gplus_large_deletion.edges", vertex_file_path="data/gplus.vertices";
  int n=102100, c=200, d=10000;
  */

  cout<<"CCCCCCCCCCCCCCCCCCCCCCCCCCCC="<<c<<endl;
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
  if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
  cout<<"BYTES="<<BYTES/pow(1024,2)<<"mb, MAX_BYTES="<<MAX_BYTES/pow(1024,2)<<"mb"<<endl;

  return 0;
}

/*----------------*
 * MAIN ALGORITHM *
 *----------------*/

int vertex_sampling(string edge_file_path, string vertex_file_path, int num_vertices, int c, int d, set<vertex> &neighbourhood, vertex& root) {
  BYTES=0; MAX_BYTES=0;
  int x=max((int)num_vertices/c,(int)sqrt(num_vertices));
  int sample_size=(int)(10*log(num_vertices)*x);
  BYTES+=2*sizeof(int);

  cout<<"n="<<num_vertices<<", sample_size="<<sample_size<<endl<<"proportion="<<(float)sample_size/(float)num_vertices<<endl;

  set<vertex> vertex_sample, vertex_not_in_sample;
  map<vertex,int> vertex_to_index; map<int,vertex> index_to_vertex;
  BYTES+=2*sizeof(set<vertex>)+sizeof(map<int,vertex>)+sizeof(map<vertex,int>);

  sample_vertices(vertex_file_path,num_vertices,sample_size,vertex_sample,vertex_not_in_sample);
  cout<<vertex_sample.size()<<","<<vertex_not_in_sample.size()<<endl;

  index_vertices(vertex_sample,vertex_not_in_sample,vertex_to_index,index_to_vertex); // TODO stop requirement for double here
  map<vertex,bool*> edge_vectors=initialise_edge_vectors(vertex_sample,num_vertices);
  construct_edge_vectors(edge_file_path,vertex_sample,num_vertices,vertex_to_index,edge_vectors);

  int neighbourhood_size=d/c;
  int num_samplers=10*(neighbourhood_size)*log(num_vertices);
  cout<<"d/c="<<neighbourhood_size<<endl<<"# samplers="<<num_samplers<<endl;
  int result=run_l0_samplers(vertex_sample,edge_vectors,num_vertices,index_to_vertex,num_samplers,neighbourhood_size,neighbourhood,root);
  BYTES+=3*sizeof(int);
  return result;
}

// generate uniform random sample of vertices
void sample_vertices(string vertex_file_path, int num_vertices, int sample_size, set<vertex>& vertex_sample, set<vertex>& vertex_not_in_sample) {

  // choose which vertices to sample
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
  uniform_int_distribution<int> uniform_d(0,num_vertices-1);
  BYTES+=sizeof(default_random_engine)+sizeof(uniform_int_distribution<int>);

  set<int> sample_indicies;
  BYTES+=sizeof(set<int>);
  while(sample_indicies.size()<sample_size) sample_indicies.insert(uniform_d(generator));

  // sample vertices from file
  string line; int i=0;
  set<int>::iterator it=sample_indicies.begin();
  ifstream stream(vertex_file_path);
  while (getline(stream,line)) {
    if(*it==i) { // sample vertex
      it++; // next index to find
      vertex_sample.insert(line);
      BYTES+=sizeof(line);
    } else {
      vertex_not_in_sample.insert(line);
      BYTES+=sizeof(line);
    }
    i++;
  }

  BYTES-=sizeof(default_random_engine)+sizeof(uniform_int_distribution<int>);

}

// map vertices to an index used for edge vectors
void index_vertices(set<vertex> vertex_sample, set<vertex> vertex_not_in_sample, map<vertex,int>& vertex_to_index, map<int,vertex>& index_to_vertex) {
  int i=0;
  for (set<vertex>::iterator it=vertex_sample.begin();it!=vertex_sample.end(); it++) {
    vertex_to_index[*it]=i;
    index_to_vertex[i]=*it;
    BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*);
    i++;
  }

  for (set<vertex>::iterator it=vertex_not_in_sample.begin();it!=vertex_not_in_sample.end(); it++) {
    vertex_to_index[*it]=i;
    index_to_vertex[i]=*it;
    BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*);
    i++;
  }
}

// initalise edge vector for each vertex in sample
map<vertex,bool*> initialise_edge_vectors(set<vertex> vertex_sample, int num_vertices) {
  map<vertex,bool*> edge_vectors;
  for (set<vertex>::iterator it=vertex_sample.begin();it!=vertex_sample.end(); it++) {
    bool* zeroes=new bool[num_vertices]; // initalise zero vector
    for (int i=0; i<num_vertices; i++) zeroes[i]=false;
    edge_vectors[*it]=zeroes;
    BYTES+=sizeof(vertex)+num_vertices*sizeof(bool)+sizeof(void*);
  }
  return edge_vectors;
}

// construct edge vector for each sample vertex
void construct_edge_vectors(string edge_file_path, set<vertex> vertex_sample, int num_vertices, map<vertex,int> vertex_to_index, map<vertex,bool*>& edge_vectors) {
  ifstream stream(edge_file_path);
  string line; edge e;
  BYTES+=sizeof(string)+sizeof(edge);
  if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
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
  BYTES-=sizeof(string)+sizeof(edge);
}

// run l0 samplers on edge vertices to try and reconstruct a neighbourhood
int run_l0_samplers(set<vertex> vertex_sample, map<vertex,bool*> edge_vectors, int num_vectors, map<int,vertex> index_to_vertex, int num_samplers, int neighbourhood_size, set<vertex>& neighbourhood, vertex&root) {
  // run l0 sampler on each v in sample
  for (set<vertex>::iterator it=vertex_sample.begin(); it!=vertex_sample.end(); it++) {
    bool* edge_vector=edge_vectors[*it];
    root=*it;
    BYTES-=sizeof(vertex)*neighbourhood.size();
    neighbourhood.clear();
    for (int i=0; i<num_samplers; i++) {
      int sampled_index=l0_sample(edge_vector,num_vectors);
      if (sampled_index!=-1) {
        neighbourhood.insert(index_to_vertex[sampled_index]);
        BYTES+=sizeof(vertex);
        if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
      }
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
int l0_sample(bool* a, int size) {
  int j=1;
  bool* a_j=copy_arr(a,size);

  hash_params ps;

  BYTES+=sizeof(bool)*(size+1)+sizeof(hash_params);
  if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;

  while (true) {
    //cout<<"j="<<j<<endl;
    ps=generate_hash(pow(size,3)); // m=n^3 where n is |x|

    sample(a_j,size,j,ps); // sample from vector
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
  delete [] a_j;
  BYTES-=sizeof(bool)*(size+1)+sizeof(hash_params);
}

// sample from vector a with condition j
void sample(bool* a, int size, int j, hash_params ps) {
  for (int i=0; i<size; i+=1) {
    int h=hash_function(i,ps);
    if (h<=ps.m/pow(2,j)) a[i]=a[i];
    else a[i]=0;
  }
}

// recover a non-zero index of a if a is 1-sparse
int recover(bool* aj, int size) {
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
bool* copy_arr(bool* arr, int size) {
  bool* new_arr=new bool[size];
  for (int i=0; i<size; i++) new_arr[i]=arr[i];
  return new_arr;
}
