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
#include <set>
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
void display_results(int c, int d, int n, string file_name);

// main algorithm
int single_pass_insertion_stream(int c, int d, int n,ifstream& stream, vector<vertex>& neighbourhood, vertex& root);

// resevoir sampling
void update_resevoir(vertex n, int d1, int d2, int count, int size, vector<vertex>& resevoir, vector<edge>& edges);

// utility
void parse_edge(string str, edge& e);
double variance(vector<int> vals);

/*-----*
* BODY *
*------*/

int main() {
  //int d=586, n=747, reps=100; int c=3;
  //display_results(c,d,n,"../../data/facebook.edges");
  //execute_test(2,100,1,reps,d,n,"../../data/facebook.edges","../../results/facebook_results.csv");
  int d=5948, n=12417, reps=10; // NOTE - # edges=1,179,613
  execute_test(3,20,1,reps,d,n,"../../data/gplus.edges","results_not_quit_early.csv");
  //int d=104947, n=102100, reps=10; // c=runs, d/c=d2, n=# vertices, NOTE - set d=max degree, n=number of vertices
  //execute_test(9,17,8,reps,d,n,"../../data/gplus_large.edges","../../results/gplus_large_results.csv");
  return 0;
}

void display_results(int c, int d, int n, string file_name) {
  vector<vertex> neighbourhood; vertex root; // variables for returned values
  ifstream stream(file_name);
  BYTES=0; RESEVOIR_BYTES=0; DEGREE_BYTES=0; MAX_BYTES=0; MAX_RESEVOIR_BYTES=0;

  time_point before=chrono::high_resolution_clock::now(); // time before execution
  int edges_checked=single_pass_insertion_stream(c,d,n,stream,neighbourhood,root);
  time_point after=chrono::high_resolution_clock::now(); // time after execution

  cout<<"Root Node - "<<root<<endl;
  cout<<"Neighbourhood Size - "<<neighbourhood.size()<<endl;
  cout<<"# edges checked - "<<edges_checked<<endl;
  cout<<"Time - "<<chrono::duration_cast<chrono::seconds>(after-before).count()<<" seconds"<<endl;

  if (RESEVOIR_BYTES>MAX_RESEVOIR_BYTES) MAX_RESEVOIR_BYTES=RESEVOIR_BYTES;

  cout<<"MAX RESEVOIR - "<<(float)MAX_RESEVOIR_BYTES/1048576<<" mb"<<endl;
  cout<<"MAX DEGREE - "<<(float)DEGREE_BYTES/1048576<<" mb"<<endl;
  //for (vector<vertex>::iterator i=neighbourhood.begin(); i!=neighbourhood.end(); i++) cout<<*i<<",";
}

// Runs algorithm multiple time, writing results to a csv file
void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string file_name, string out_file) {
  ofstream outfile(out_file);
  outfile<<"name,"<<file_name<<endl<<"n,"<<n<<endl<<"d,"<<d<<endl<<"repetitions,"<<reps<<endl<<endl; // test details
  outfile<<"c,time (microseconds),mean max space (bytes),mean resevoir space (bytes),mean degree space (bytes),mean edges checked,variance time, variance max space,variance resevoir space, variance degree space,varriance edges checked,successes"<<endl; // headers
  vector<vertex> neighbourhood; vertex root; // variables for returned values
  vector<int> times, total_space, resevoir_space, degree_space, edges_checked; // results of each run of c
  int successes;
  for (int c=c_min;c<=c_max;c+=c_step) {
    successes=0;
    times.clear(); total_space.clear(); resevoir_space.clear(); degree_space.clear(); edges_checked.clear();// reset for new run of c
    for (int i=0;i<reps;i++) {
      cout<<"("<<i<<"/"<<reps<<") "<<c<<"/"<<c_max<<endl; // output to terminal

      // reset values
      BYTES=0; RESEVOIR_BYTES=0; DEGREE_BYTES=0; MAX_BYTES=0; MAX_RESEVOIR_BYTES=0;
      neighbourhood.clear(); vertex* p=&root; p=nullptr;
      ifstream stream(file_name); // file to read

      time_point before=chrono::high_resolution_clock::now(); // time before execution
      edges_checked.push_back(single_pass_insertion_stream(c,d,n,stream,neighbourhood,root));
      time_point after=chrono::high_resolution_clock::now(); // time after execution

      cout<<root<<endl;

      stream.close();

      if (neighbourhood.size()!=0) successes+=1;

      auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count(); // time passed
      if (RESEVOIR_BYTES>MAX_RESEVOIR_BYTES) MAX_RESEVOIR_BYTES=RESEVOIR_BYTES;
      times.push_back(duration); total_space.push_back(MAX_BYTES); resevoir_space.push_back(MAX_RESEVOIR_BYTES); degree_space.push_back(DEGREE_BYTES);
    }
    int mean_duration      =accumulate(times.begin(),times.end(),0)/times.size();
    int mean_total_space   =accumulate(total_space.begin(),total_space.end(),0)/total_space.size();
    int mean_resevoir_space=accumulate(resevoir_space.begin(),resevoir_space.end(),0)/resevoir_space.size();
    int mean_degree_space  =accumulate(degree_space.begin(),degree_space.end(),0)/degree_space.size();
    int mean_edges_checked =accumulate(edges_checked.begin(),edges_checked.end(),0)/edges_checked.size();

    double variance_duration      =variance(times);
    double variance_total_space   =variance(total_space);
    double variance_resevoir_space=variance(resevoir_space);
    double variance_degree_space  =variance(degree_space);
    double variance_edges_checked =variance(edges_checked);

    outfile<<c<<","<<mean_duration<<","<<mean_total_space<<","<<mean_resevoir_space<<","<<mean_degree_space<<","<<mean_edges_checked<<","<<variance_duration<<","<<variance_total_space<<","<<variance_resevoir_space<<","<<variance_degree_space<<","<<variance_edges_checked<<","<<successes<<endl; // write values to file
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
    if (edge_count%10000==0) cout<<"\r"<<edge_count;

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
      } else if (find(resevoirs[j]->begin(),resevoirs[j]->end(),e.snd)!=resevoirs[j]->end()) { // if second endpoint is in resevoir
        if (degrees[e.snd]<=d2+d1) {
          edges[j]->push_back(e);
          BYTES+=sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESEVOIR_BYTES+=sizeof(edge);
        }
      }
      BYTES-=sizeof(int); // deletion of d1
    }

  }
  cout<<"\rDONE                         "<<endl;

  set<vertex> successful_in_run; // set so each root is chosen uniformly
  set<pair<int,vertex> > successful_overall;
  BYTES+=sizeof(set<vertex>)+sizeof(set<pair<int,vertex> >); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
  RESEVOIR_BYTES+=sizeof(set<vertex>)+sizeof(set<pair<int,vertex> >);

  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time

  for (int j=0; j<c; j++) { // find all successful runs
    int d1=max(1,(j*d)/c), d2=d/c; // calculate degree bounds for run
    BYTES-=sizeof(vertex)*successful_in_run.size();
    RESEVOIR_BYTES-=sizeof(vertex)*successful_in_run.size();
    successful_in_run.clear();

    // find all successful runs for this resevoir sampler
    for (vector<vertex>::iterator it=resevoirs[j]->begin(); it!=resevoirs[j]->end(); it++) {
      if (degrees[*it]>=d1+d2) successful_in_run.insert(*it);
    }
    BYTES+=sizeof(vertex)*successful_in_run.size(); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
    RESEVOIR_BYTES+=sizeof(vertex)*successful_in_run.size();

    // uniformly choose one at random
    uniform_int_distribution<int> d(0,successful_in_run.size()-1);
    int chosen_index=d(generator);
    set<vertex>::iterator it=successful_in_run.begin();
    advance(it,chosen_index);
    vertex chosen_vertex=*it;
    pair<int,vertex> p(j,chosen_vertex);
    successful_overall.insert(p);
    BYTES+=sizeof(int)+sizeof(vertex)+sizeof(pair<int,vertex>); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
    RESEVOIR_BYTES+=sizeof(vertex)*successful_in_run.size();
  }

  if (successful_overall.size()!=0) { // successful run found
    // choose uniformly at random one of the successful runs
    uniform_int_distribution<int> overall_d(0,successful_overall.size()-1);
    int chosen_index=overall_d(generator);

    set<pair<int,vertex> >::iterator it=successful_overall.begin();
    advance(it,chosen_index);
    int chosen_run=it->first;
    vertex chosen_vertex=it->second;

    for (vector<edge>::iterator i=edges[chosen_run]->begin(); i!=edges[chosen_run]->end(); i++) { // construct neighbourhood to be returned
      if (i->fst==chosen_vertex) neighbourhood.push_back(i->snd);
      else if (i->snd==chosen_vertex) neighbourhood.push_back(i->fst);
    }

    BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
    RESEVOIR_BYTES+=neighbourhood.size()*sizeof(edge);
    root=e.snd;
    return edge_count;
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

// return variance of values in a vector
double variance(vector<int> vals) {
  double var=0;
  double mean=accumulate(vals.begin(),vals.end(),0)/vals.size();

  for (vector<int>::iterator it=vals.begin(); it!=vals.end(); it++) var+=(*it-mean)*(*it-mean);
  var/=(vals.size()-1);

  return var;
}
