/*
 *  AN L0-SAMPLER FOR A STREAM
 *  uniformly samples a term from a stream
 *  --------------------------
 *  NOTES
 *    Probability of index i being sampled from vector |x| with x in R^n
 *      depends on how many non-zero terms are in x, and not on which terms
 *      are non-zero.
 *   For a stream, all entries are considered non-zero so the probability
 *      of index i being sampled depends on n (the length of the stream) only
 */

#include <chrono>
#include <iostream>
#include <fstream>
#include <random>
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
  string edge_file_path="../../data/facebook_deletion.edges";
  int n=66643; // number of edges in file
  vector<string> sample;
  hash_params ps=generate_hash(10*n);

  ifstream stream(edge_file_path);
  int j=22; // level //TODO this is a bodge due to issue with hash function (FIND BETTER HASH FUNCTION)
  int i=0; string line; unsigned long lim=(10*n)/pow(2,j);
  while (getline(stream,line)) { // single pass of stream
    unsigned long h=hash_function(i,ps);
    if (h<=lim) sample.push_back(line);
    i+=1;
  }

  while(sample.size()>1) {
    ps=generate_hash(10*n); // change hash function
    j+=1;
    vector<string> new_sample;
    lim=(10*n)/pow(2,j);
    for (int i=0; i<sample.size(); i++) {
      int h=hash_function(i,ps);
      if (h<=lim) new_sample.push_back(sample[i]);
    }
    sample=new_sample;
  }

  if (sample.size()==1) {
    cout<<"SUCCESS "<<sample[0]<<endl;
  } else {
    cout<<"FAIL"<<endl;
  }

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
