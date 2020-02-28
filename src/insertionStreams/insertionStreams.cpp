/*-----------------*
 * One pass c-approximation Streaming Algorithm for Neighbourhood Detection for Insertion-Only Graph Streams
 *-----------------*/

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <vector>

using namespace std;

int BYTES; // space used atm
int RESEVOIR_BYTES; // space used by resevoirs atm
int DEGREE_BYTES; // space used to store degrees of vertices
int MAX_BYTES; // max space used at any time
int MAX_RESEVOIR_BYTES; // max space used by resevoirs

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

using vertex = string; // typemap vertex
using time_point=chrono::high_resolution_clock::time_point;

struct edge { // undirected edge
  vertex fst;
  vertex snd;
};

/*-----------*
* SIGNATURES *
*------------*/

// size = size of resevoir
void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string file_name, string out_file); // for (int c=c_min;c<=c_max;c+=c_step). Reps is the number of times each c is tested, average is taken.
int single_pass_insertion_stream(int c, int d, int n,ifstream& stream, vector<vertex>& neighbourhood, vertex& root);
void update_resevoir(vertex n, int d1, int d2, int count, int size, vector<vertex>& resevoir, vector<edge>& edges);
void parse_edge(string str, edge& e);

/*-----*
* BODY *
*------*/

int main() {
  //int d=586, n=747, reps=100;
  //execute_test(2,20,1,reps,d,n,"../../data/facebook.edges","results/facebook_results.csv");
  //int d=5948, n=12417, reps=20;
  //execute_test(2,2,1,reps,d,n,"../../data/gplus.edges","results/gplus_results.csv");
  //int d=104947, n=102100, reps=10; // c=runs, d/c=d2, n=# vertices, NOTE - set d=max degree, n=number of vertices
  //execute_test(9,17,8,reps,d,n,"../../data/gplus_large.edges","results/gplus_large_results.csv");
  return 0;
}

// Runs algorithm multiple time, writing results to a csv file
void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string file_name, string out_file) {
  ofstream outfile(out_file);
  outfile<<"name,"<<file_name<<endl<<"n,"<<n<<endl<<"d,"<<d<<endl<<"repetitions,"<<reps<<endl<<endl; // test details
  outfile<<"c,time (microseconds),max space (bytes),resevoir space (bytes),degree space (bytes),mean edges checked,successes"<<endl; // headers
  vector<vertex> neighbourhood; vertex root; // variables for returned values
  vector<int> times, total_space, resevoir_space, degree_space; // results of each run of c
  int successes; long edges_checked=0;
  for (int c=c_min;c<=c_max;c+=c_step) {
    successes=0; edges_checked=0;
    times.clear(); total_space.clear(); resevoir_space.clear(); degree_space.clear();// reset for new run of c
    for (int i=0;i<reps;i++) {
      cout<<"("<<i<<"/"<<reps<<") "<<c<<"/"<<c_max<<endl; // output to terminal

      // reset values
      BYTES=0; RESEVOIR_BYTES=0; DEGREE_BYTES=0; MAX_BYTES=0; MAX_RESEVOIR_BYTES=0;
      neighbourhood.clear(); vertex* p=&root; p=nullptr;
      ifstream stream(file_name); // file to read

      time_point before=chrono::high_resolution_clock::now(); // time before execution
      edges_checked+=single_pass_insertion_stream(c,d,n,stream,neighbourhood,root);
      time_point after=chrono::high_resolution_clock::now(); // time after execution

      if (neighbourhood.size()!=0) successes+=1;

      stream.close();
      auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count(); // time passed
      if (RESEVOIR_BYTES>MAX_RESEVOIR_BYTES) MAX_RESEVOIR_BYTES=RESEVOIR_BYTES;
      times.push_back(duration); total_space.push_back(MAX_BYTES); resevoir_space.push_back(MAX_RESEVOIR_BYTES); degree_space.push_back(DEGREE_BYTES);
    }
    int average_duration = accumulate(times.begin(),times.end(),0)/reps;
    int average_total_space = accumulate(total_space.begin(),total_space.end(),0)/reps;
    int average_resevoir_space = accumulate(resevoir_space.begin(),resevoir_space.end(),0)/reps;
    int average_degree_space = accumulate(degree_space.begin(),degree_space.end(),0)/reps;
    outfile<<c<<","<<average_duration<<","<<average_total_space<<","<<average_resevoir_space<<","<<average_degree_space<<","<<(int) edges_checked/reps<<","<<successes<<endl; // write values to file
  }
  outfile.close();
}

// perform resevoir sampling
// returns number of edges which are read
int single_pass_insertion_stream(int c, int d, int n, ifstream& stream, vector<vertex>& neighbourhood, vertex& root) {
  int size=ceil(log10(n)*pow(n,(double)1/c));

  // initalise edge & vector sets for each parallel run
  vector<vertex>* resevoirs[c]; vector<edge>* edges[c];
  for (int i=0; i<c; i++) {
    resevoirs[i]=new vector<vertex>;
    edges[i]=new vector<edge>;
  }
  BYTES+=c*(sizeof(vector<vertex>)+sizeof(vector<edge>));
  RESEVOIR_BYTES+=c*(sizeof(vector<vertex>)+sizeof(vector<edge>));

  string line; edge e; map<vertex,int> degrees; // these are shared for each run
  int count[c];  // counts the number of vertexs >=d1
  BYTES+=sizeof(string)+sizeof(edge)+sizeof(map<vertex,int>)+sizeof(int);
  DEGREE_BYTES+=sizeof(map<vertex,int>);

  int edge_count=0;
  while (getline(stream,line)) { // While stream is not empty
    parse_edge(line,e);
    edge_count+=1;

    // increment degrees for each vertex
    if (degrees.count(e.fst)) {
      degrees[e.fst]+=1;
    } else {
      degrees[e.fst]=1;
      BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*); // sizeof(void*)=size of pointer
      DEGREE_BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*);
    }

    if (degrees.count(e.snd)) {
      degrees[e.snd]+=1;
    } else {
      degrees[e.snd]=1;
      BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*); // sizeof(void*)=size of pointer
      DEGREE_BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*);
    }

    for (int j=0; j<c; j++) { // perform parallel runs
      int d1=max(1,(j*d)/c), d2=d/c; // calculate degree bounds for run
      BYTES+=sizeof(int);
      // NB rest is standard degree-restricted sampling

      // Consider adding first vertex to the resevoir
      if (degrees[e.fst]==d1) {
        count[j]+=1; // increment number of d1 degree vertexs
        update_resevoir(e.fst,d1,d2,count[j],size,*resevoirs[j],*edges[j]); // possibly add value to resevoir
      }

      // Consider adding second vertex to the resevoir
      if (degrees[e.snd]==d1) {
        count[j]+=1; // increment number of d1 degree vertexs
        update_resevoir(e.snd,d1,d2,count[j],size,*resevoirs[j],*edges[j]); // possibly add value to resevoir
      }

      if (find(resevoirs[j]->begin(),resevoirs[j]->end(),e.fst)!=resevoirs[j]->end()) { // if first endpoint is in resevoir
        if (degrees[e.fst]<=d2+d1) {
          edges[j]->push_back(e);
          BYTES+=sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESEVOIR_BYTES+=sizeof(edge);
        }
        if (degrees[e.fst]==d2+d1) { // sufficient neighbourhood has been found, return it
          for (vector<edge>::iterator i=edges[j]->begin(); i!=edges[j]->end(); i++) { // construct neighbourhood to be returned
            if (i->fst==e.fst) neighbourhood.push_back(i->snd);
            else if (i->snd==e.fst) neighbourhood.push_back(i->fst);
          }
          BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESEVOIR_BYTES+=neighbourhood.size()*sizeof(edge);
          root=e.fst;
          return edge_count;
        }
      } else if (find(resevoirs[j]->begin(),resevoirs[j]->end(),e.snd)!=resevoirs[j]->end()) { // if second endpoint is in resevoir
        if (degrees[e.snd]<=d2+d1) {
          edges[j]->push_back(e);
          BYTES+=sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESEVOIR_BYTES+=sizeof(edge);
        }
        if (degrees[e.snd]==d2+d1) { // sufficient neighbourhood has been found, return it
          for (vector<edge>::iterator i=edges[j]->begin(); i!=edges[j]->end(); i++) { // construct neighbourhood to be returned
            if (i->fst==e.snd) neighbourhood.push_back(i->snd);
            else if (i->snd==e.snd) neighbourhood.push_back(i->fst);
          }
          BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESEVOIR_BYTES+=neighbourhood.size()*sizeof(edge);
          root=e.snd;
          return edge_count;
        }
      }
      BYTES-=sizeof(int); // deletion of d1
    }

  }

  // No sucessful runs
  neighbourhood.clear();
  vertex* p=&root;
  p=nullptr;

  return edge_count;
}

void update_resevoir(vertex n, int d1, int d2, int count, int size, vector<vertex>& resevoir, vector<edge>& edges) {
  if (resevoir.size()<size) { // resevoir is not full
    resevoir.push_back(n);
    BYTES+=sizeof(vertex); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
    RESEVOIR_BYTES+=2*sizeof(vertex);
  } else { // resevoir is full
    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
    bernoulli_distribution bernoulli_d((float)size/(float)count);
    BYTES+=sizeof(default_random_engine)+sizeof(bernoulli_distribution);
    if (bernoulli_d(generator)) { // if coin flip passes

      uniform_int_distribution<unsigned long> uniform_d(0,size-1); // decide which vertex to delete
      int to_delete=uniform_d(generator);

      vertex to_delete_val=resevoir[to_delete];
      BYTES+=2*sizeof(vertex);
      RESEVOIR_BYTES+=2*sizeof(vertex);
      if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;

      resevoir.erase(resevoir.begin()+to_delete);
      resevoir.push_back(n);
      // NB No space change

      vector<edge>::iterator i=edges.begin();
      while (i!=edges.end()) {
        if (i->fst==to_delete_val || i->snd==to_delete_val) {
          i=edges.erase(i);
          if (RESEVOIR_BYTES>MAX_RESEVOIR_BYTES) MAX_RESEVOIR_BYTES=RESEVOIR_BYTES;
          BYTES-=sizeof(edge);RESEVOIR_BYTES-=sizeof(edge);
        }
        else i++;
      }
      BYTES-=2*sizeof(int);
    }
    BYTES-=sizeof(default_random_engine)+sizeof(bernoulli_distribution);
  }
}

// parse ege from stream
void parse_edge(string str, edge& e) {
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
  //e.fst=stoi(fst);
  //e.snd=stoi(snd);
  e.fst=fst;e.snd=snd;
}
