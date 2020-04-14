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
int RESERVOIR_BYTES; // space used by reservoirs atm
int DEGREE_BYTES; // space used to store degrees of vertices
int MAX_BYTES; // max space used at any time
int MAX_RESERVOIR_BYTES; // max space used by reservoirs

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

// size = size of reservoir
void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string file_name, string out_file); // for (int c=c_min;c<=c_max;c+=c_step). Reps is the number of times each c is tested, average is taken.
void display_results(int c, int d, int n, string file_name);

// main algorithm
int single_pass_insertion_stream(int c, int d, int n,ifstream& stream, vector<vertex>& neighbourhood, vertex& root);

// reservoir sampling
void update_reservoir(vertex n, int d1, int d2, int count, int size, vector<vertex>& reservoir, vector<edge>& edges);

// utility
void parse_edge(string str, edge& e);
double variance(vector<int> vals);

/*-----*
* BODY *
*------*/

int main() {
  //int d=586, n=747, reps=100; int c=3;
  //display_results(c,d,n,"../../data/facebook.edges");

  string out_file_path="results_quit_early.csv";
  int n,d,reps; string edge_file_path;

  // c=runs, d/c=d2, n=# vertices, NOTE - set d=max degree, n=number of vertices
  //n=12417; d=5948; reps=10; edge_file_path="../../data/gplus.edges"; // NOTE - # edges=1,179,613
  //n=102100; d=104947; reps=10; edge_file_path="../../data/gplus_large.edges"; // NOTE - # edges=30,238,035
  //n=52; d=35; reps=10; edge_file_path="../../data/facebook_small.edges"; // NOTE - # edges=292
  //n=747; d=586; reps=10; edge_file_path="../../data/facebook.edges"; // NOTE - # edges=60,050
  //n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000star.edges"; // NOTE - # edges=999
  //n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000complete.edges"; // NOTE - # edges=499,500
  //execute_test(3,20,1,reps,d,n,edge_file_path,out_file_path);

  n=102100; d=104947; reps=10; edge_file_path="../../data/gplus_large.edges"; // NOTE - # edges=30,238,035
  out_file_path="results_quit_early_gplus_large.csv";
  execute_test(3,20,1,reps,d,n,edge_file_path,out_file_path);

  return 0;
}

void display_results(int c, int d, int n, string file_name) {
  vector<vertex> neighbourhood; vertex root; // variables for returned values
  ifstream stream(file_name);
  BYTES=0; RESERVOIR_BYTES=0; DEGREE_BYTES=0; MAX_BYTES=0; MAX_RESERVOIR_BYTES=0;

  time_point before=chrono::high_resolution_clock::now(); // time before execution
  int edges_checked=single_pass_insertion_stream(c,d,n,stream,neighbourhood,root);
  time_point after=chrono::high_resolution_clock::now(); // time after execution

  cout<<"Root Node - "<<root<<endl;
  cout<<"Neighbourhood Size - "<<neighbourhood.size()<<endl;
  cout<<"# edges checked - "<<edges_checked<<endl;
  cout<<"Time - "<<chrono::duration_cast<chrono::seconds>(after-before).count()<<" seconds"<<endl;

  if (RESERVOIR_BYTES>MAX_RESERVOIR_BYTES) MAX_RESERVOIR_BYTES=RESERVOIR_BYTES;

  cout<<"MAX RESERVOIR - "<<(float)MAX_RESERVOIR_BYTES/1048576<<" mb"<<endl;
  cout<<"MAX DEGREE - "<<(float)DEGREE_BYTES/1048576<<" mb"<<endl;
  //for (vector<vertex>::iterator i=neighbourhood.begin(); i!=neighbourhood.end(); i++) cout<<*i<<",";
}

// Runs algorithm multiple time, writing results to a csv file
void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string file_name, string out_file) {
  ofstream outfile(out_file);
  outfile<<"name,"<<file_name<<endl<<"n,"<<n<<endl<<"d,"<<d<<endl<<"repetitions,"<<reps<<endl<<endl; // test details
  outfile<<"c,time (microseconds),mean max space (bytes),mean reservoir space (bytes),mean degree space (bytes),mean edges checked,variance time, variance max space,variance reservoir space, variance degree space,varriance edges checked,successes"<<endl; // headers
  vector<vertex> neighbourhood; vertex root; // variables for returned values
  vector<int> times, total_space, reservoir_space, degree_space, edges_checked; // results of each run of c
  int successes;
  for (int c=c_min;c<=c_max;c+=c_step) {
    successes=0;
    times.clear(); total_space.clear(); reservoir_space.clear(); degree_space.clear(); edges_checked.clear();// reset for new run of c
    for (int i=0;i<reps;i++) {
      cout<<"("<<i<<"/"<<reps<<") "<<c<<"/"<<c_max<<endl; // output to terminal

      // reset values
      BYTES=0; RESERVOIR_BYTES=0; DEGREE_BYTES=0; MAX_BYTES=0; MAX_RESERVOIR_BYTES=0;
      neighbourhood.clear(); vertex* p=&root; p=nullptr;
      ifstream stream(file_name); // file to read

      time_point before=chrono::high_resolution_clock::now(); // time before execution
      edges_checked.push_back(single_pass_insertion_stream(c,d,n,stream,neighbourhood,root));
      time_point after=chrono::high_resolution_clock::now(); // time after execution

      cout<<root<<endl;

      stream.close();

      if (neighbourhood.size()!=0) successes+=1;

      auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count(); // time passed
      if (RESERVOIR_BYTES>MAX_RESERVOIR_BYTES) MAX_RESERVOIR_BYTES=RESERVOIR_BYTES;
      times.push_back(duration); total_space.push_back(MAX_BYTES); reservoir_space.push_back(MAX_RESERVOIR_BYTES); degree_space.push_back(DEGREE_BYTES);
    }
    int mean_duration      =accumulate(times.begin(),times.end(),0)/times.size();
    int mean_total_space   =accumulate(total_space.begin(),total_space.end(),0)/total_space.size();
    int mean_reservoir_space=accumulate(reservoir_space.begin(),reservoir_space.end(),0)/reservoir_space.size();
    int mean_degree_space  =accumulate(degree_space.begin(),degree_space.end(),0)/degree_space.size();
    int mean_edges_checked =accumulate(edges_checked.begin(),edges_checked.end(),0)/edges_checked.size();

    double variance_duration      =variance(times);
    double variance_total_space   =variance(total_space);
    double variance_reservoir_space=variance(reservoir_space);
    double variance_degree_space  =variance(degree_space);
    double variance_edges_checked =variance(edges_checked);

    outfile<<c<<","<<mean_duration<<","<<mean_total_space<<","<<mean_reservoir_space<<","<<mean_degree_space<<","<<mean_edges_checked<<","<<variance_duration<<","<<variance_total_space<<","<<variance_reservoir_space<<","<<variance_degree_space<<","<<variance_edges_checked<<","<<successes<<endl; // write values to file
  }
  outfile.close();
}

// perform reservoir sampling
// returns number of edges which are read
int single_pass_insertion_stream(int c, int d, int n, ifstream& stream, vector<vertex>& neighbourhood, vertex& root) {
  int size=ceil(log10(n)*pow(n,(double)1/c));

  // initalise edge & vector sets for each parallel run
  vector<vertex>* reservoirs[c]; vector<edge>* edges[c];
  for (int i=0; i<c; i++) {
    reservoirs[i]=new vector<vertex>;
    edges[i]=new vector<edge>;
  }
  BYTES+=c*(sizeof(vector<vertex>)+sizeof(vector<edge>));
  RESERVOIR_BYTES+=c*(sizeof(vector<vertex>)+sizeof(vector<edge>));

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

      // Consider adding first vertex to the reservoir
      if (degrees[e.fst]==d1) {
        count[j]+=1; // increment number of d1 degree vertexs
        update_reservoir(e.fst,d1,d2,count[j],size,*reservoirs[j],*edges[j]); // possibly add value to reservoir
      }

      // Consider adding second vertex to the reservoir
      if (degrees[e.snd]==d1) {
        count[j]+=1; // increment number of d1 degree vertexs
        update_reservoir(e.snd,d1,d2,count[j],size,*reservoirs[j],*edges[j]); // possibly add value to reservoir
      }

      if (find(reservoirs[j]->begin(),reservoirs[j]->end(),e.fst)!=reservoirs[j]->end()) { // if first endpoint is in reservoir
        if (degrees[e.fst]<=d2+d1) {
          edges[j]->push_back(e);
          BYTES+=sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=sizeof(edge);
        }
        if (degrees[e.fst]==d2+d1) { // sufficient neighbourhood has been found, return it
          for (vector<edge>::iterator i=edges[j]->begin(); i!=edges[j]->end(); i++) { // construct neighbourhood to be returned
            if (i->fst==e.fst) neighbourhood.push_back(i->snd);
            else if (i->snd==e.fst) neighbourhood.push_back(i->fst);
          }
          BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=neighbourhood.size()*sizeof(edge);
          root=e.fst;
          return edge_count;
        }
      } else if (find(reservoirs[j]->begin(),reservoirs[j]->end(),e.snd)!=reservoirs[j]->end()) { // if second endpoint is in reservoir
        if (degrees[e.snd]<=d2+d1) {
          edges[j]->push_back(e);
          BYTES+=sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=sizeof(edge);
        }
        if (degrees[e.snd]==d2+d1) { // sufficient neighbourhood has been found, return it
          for (vector<edge>::iterator i=edges[j]->begin(); i!=edges[j]->end(); i++) { // construct neighbourhood to be returned
            if (i->fst==e.snd) neighbourhood.push_back(i->snd);
            else if (i->snd==e.snd) neighbourhood.push_back(i->fst);
          }
          BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=neighbourhood.size()*sizeof(edge);
          root=e.snd;
          return edge_count;
        }
      }
      BYTES-=sizeof(int); // deletion of d1
    }

  }
  cout<<"\rDONE                         "<<endl;

  // No sucessful runs
  neighbourhood.clear();
  vertex* p=&root;
  p=nullptr;

  return edge_count;
}

void update_reservoir(vertex n, int d1, int d2, int count, int size, vector<vertex>& reservoir, vector<edge>& edges) {
  if (reservoir.size()<size) { // reservoir is not full
    reservoir.push_back(n);
    BYTES+=sizeof(vertex); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
    RESERVOIR_BYTES+=2*sizeof(vertex);
  } else { // reservoir is full
    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
    bernoulli_distribution bernoulli_d((float)size/(float)count);
    BYTES+=sizeof(default_random_engine)+sizeof(bernoulli_distribution);
    if (bernoulli_d(generator)) { // if coin flip passes
      uniform_int_distribution<unsigned long> uniform_d(0,size-1); // decide which vertex to delete
      int to_delete=uniform_d(generator);

      vertex to_delete_val=reservoir[to_delete];
      BYTES+=2*sizeof(vertex);
      RESERVOIR_BYTES+=2*sizeof(vertex);
      if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;

      reservoir.erase(reservoir.begin()+to_delete);
      reservoir.push_back(n);
      // NB No space change

      // removes edges adjacent to vertex to be deleted (which are not adjacent to another vertex in the reservoir)
      vector<edge>::iterator i=edges.begin();
      while (i!=edges.end()) {

        if (i->fst==to_delete_val) {
          if (find(reservoir.begin(),reservoir.end(),i->snd)==reservoir.end()) { // edge not adjacent to another value in the reservoir
            i=edges.erase(i); // next edge
            if (RESERVOIR_BYTES>MAX_RESERVOIR_BYTES) MAX_RESERVOIR_BYTES=RESERVOIR_BYTES;
            BYTES-=sizeof(edge);RESERVOIR_BYTES-=sizeof(edge);
          } else i++; // next edge

        } else if (i->snd==to_delete_val) {
          if (find(reservoir.begin(),reservoir.end(),i->fst)==reservoir.end()) { // edge not adjacent to another value in the reservoir
            i=edges.erase(i); // next edge
            if (RESERVOIR_BYTES>MAX_RESERVOIR_BYTES) MAX_RESERVOIR_BYTES=RESERVOIR_BYTES;
            BYTES-=sizeof(edge);RESERVOIR_BYTES-=sizeof(edge);
          } else i++; // next edge

        } else i++; // next edge

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
