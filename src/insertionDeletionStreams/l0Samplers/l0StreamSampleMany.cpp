/*
 *  MULTIPLE L0-SAMPLERs APPLIED TO A STREAM (SINGLE PASS)
 *  uniformly samples a term from a stream
 *  --------------------------
 *  NOTES
 *   See l0StreamSampler.cpp
 *   The Hash function used hear (26/02) maps from [n] to [n^3]. Given n can be in the millions this breaks down due to bad choice of prime.
 *    maybe just map [n] to [n^2]. More collisions, but at least it works
 */

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <vector>

using namespace std;
int P=1073741789; // >2^30-1

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

// parameters for hash function
struct hash_params {
  unsigned long a;
  unsigned long b;
  unsigned long m;
};

/*------------*
 * SIGNATURES *
 *------------*/

// hashing
hash_params generate_hash(unsigned long m);
unsigned long hash_function(unsigned long key, hash_params ps);

/*------*
 * BODY *
 *------*/

int main() {
  string edge_file_path="..././data/facebook_deletion.edges";
  int n=66643; // number of edges
  int d=522; // max degree
  int c=200; // accuracy
  //int num_samplers=20*n*(int)(d/c)*log(n)*(((float)1/max((int)n/c,(int)sqrt(n)))+(float)(1/c));
  int num_samplers=100; // NOTE RUNS LIKE poorly FOR LARGE num_samplers (many due to first while loop which can't be optimised, i think?)
  cout<<num_samplers<<endl;

  set<string> output; hash_params ps;
  vector<string>* samples[num_samplers];
  for (int i=0; i<num_samplers; i++) samples[i]=new vector<string>;

  ifstream stream(edge_file_path);
  int i=0; int j=1;
  string line; unsigned long lim=(10*n)/pow(2,j);
  while (getline(stream,line)) { // single pass of stream
    for (int k=0; k<num_samplers; k++) {
      ps=generate_hash(10*n);
      unsigned long h=hash_function(i,ps);
      if (h<=lim) samples[k]->push_back(line);
    }
    i+=1;
  }

  for (int k=0; k<num_samplers; k++) {
    vector<string>& sample=*samples[k];
    int j=1;
    while(sample.size()>1) {
      ps=generate_hash(10*n);
      j+=1;
      vector<string> new_sample;
      lim=(10*n)/pow(2,j);
      for (int i=0; i<sample.size(); i++) {
        int h=hash_function(i,ps);
        if (h<=lim) new_sample.push_back(sample[i]);
      }
      sample=new_sample;
    }
    if (sample.size()==1) output.insert(sample[0]);
  }

  cout<<"Sample size="<<output.size()<<endl;
  for (set<string>::iterator it=output.begin(); it!=output.end(); it++) cout<<*it<<",";
  cout<<endl;

  return 0;
}

/*---------*
 * HASHING *
 *---------*/

// generate parameters to use in hash function
hash_params generate_hash(unsigned long m) {
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
unsigned long hash_function(unsigned long key, hash_params ps) {
  return (unsigned long)((ps.a*key+ps.b)%P)%ps.m;
}
