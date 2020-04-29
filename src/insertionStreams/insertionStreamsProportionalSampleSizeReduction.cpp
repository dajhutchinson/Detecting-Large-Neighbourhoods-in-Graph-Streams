/*-----------------*
 * One pass c-approximation Streaming Algorithm for Neighbourhood Detection for Insertion-Only Graph Streams
 *
 * This uses the implementation which terminates early AND has a shared edge set AND uses a reduced number of samplers.
 * Here we can test reducing the sample size by a certain proportion of its proposed value
 *-----------------*/

#include <algorithm>
#include <bitset>
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
void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string file_name, string out_file, double p_min, double p_max, double p_step); // for (int c=c_min;c<=c_max;c+=c_step). Reps is the number of times each c is tested, average is taken.
void display_results(int c, int d, int n, string file_name);

// main algorithm
// sample size = p*proposed size
int single_pass_insertion_stream(int c, int d, int n,ifstream& stream, vector<vertex>& neighbourhood, vertex& root, double prop);

// reservoir sampling
void update_reservoir(vertex n,int c,int res_num, int d1, int d2, int count, int size, vector<vertex>& reservoir, vector<edge>& edges, map<vertex,bool*>& num_reservoirs);

// utility
void parse_edge(string str, edge& e);
double variance(vector<int> vals);

/*-----*
* BODY *
*------*/

int main() {
  //int d=586, n=747, reps=100; int c=3;
  //display_results(c,d,n,"../../data/facebook.edges");

  string out_file_path="results_proportional_sample_size_reduction.csv";
  int n,d,reps; string edge_file_path;

  // c=runs, d/c=d2, n=# vertices, NOTE - set d=max degree, n=number of vertices
  //n=12417; d=5948; reps=10; edge_file_path="../../data/gplus.edges"; // NOTE - # edges=1,179,613
  //n=102100; d=104947; reps=10; edge_file_path="../../data/gplus_large.edges"; // NOTE - # edges=30,238,035
  //n=52; d=35; reps=10; edge_file_path="../../data/facebook_small.edges"; // NOTE - # edges=292
  //n=747; d=586; reps=10; edge_file_path="../../data/facebook.edges"; // NOTE - # edges=60,050
  //n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000star.edges"; // NOTE - # edges=999
  //n=1000; d=999; reps=10; edge_file_path="../../data/artifical/1000complete.edges"; // NOTE - # edges=499,500
  //execute_test(3,20,1,reps,d,n,edge_file_path,out_file_path);

  n=747; d=586; reps=10; edge_file_path="../../data/facebook.edges"; // NOTE - # edges=60,050
  out_file_path="results_proportional_sample_size_reduction_facebook.csv";
  execute_test(3,20,1,reps,d,n,edge_file_path,out_file_path,0.01,2,0.1);

  n=12417; d=5948; reps=10; edge_file_path="../../data/gplus.edges"; // NOTE - # edges=1,179,613
  out_file_path="results_proportional_sample_size_reduction_gplus.csv";
  execute_test(3,20,1,reps,d,n,edge_file_path,out_file_path,0.01,2,0.1);


  return 0;
}

// Runs algorithm multiple time, writing results to a csv file
void execute_test(int c_min, int c_max, int c_step, int reps, int d, int n, string file_name, string out_file, double p_min, double p_max, double p_step) {
  ofstream outfile(out_file);
  outfile<<"name,"<<file_name<<endl<<"n,"<<n<<endl<<"d,"<<d<<endl<<"repetitions,"<<reps<<endl<<endl; // test details
  outfile<<"c,p,time (microseconds),mean max space (bytes),mean reservoir space (bytes),mean degree space (bytes),mean edges checked,variance time, variance max space,variance reservoir space, variance degree space,varriance edges checked,successes"<<endl; // headers
  vector<vertex> neighbourhood; vertex root; // variables for returned values
  vector<int> times, total_space, reservoir_space, degree_space, edges_checked; // results of each run of c
  int successes;
  for (int c=c_min;c<=c_max;c+=c_step) {

    double prop;
    for (prop=p_max; prop>=p_min; prop-=p_step) {

      times.clear(); total_space.clear(); reservoir_space.clear(); degree_space.clear(); edges_checked.clear();// reset for new run of c
      successes=0;
      for (int i=0;i<reps;i++) {
        cout<<"("<<i<<"/"<<reps<<") "<<c<<"/"<<c_max<<" <"<<p_min<<"-"<<prop<<"-"<<p_max<<">"<<endl; // output to terminal

        // reset values
        BYTES=0; RESERVOIR_BYTES=0; DEGREE_BYTES=0; MAX_BYTES=0; MAX_RESERVOIR_BYTES=0;
        neighbourhood.clear(); vertex* p=&root; p=nullptr;
        ifstream stream(file_name); // file to read

        time_point before=chrono::high_resolution_clock::now(); // time before execution
        edges_checked.push_back(single_pass_insertion_stream(c,d,n,stream,neighbourhood,root,prop));
        time_point after=chrono::high_resolution_clock::now(); // time after execution

        cout<<root<<endl;
        cout<<neighbourhood.size()<<"("<<d/c<<")"<<endl;

        stream.close();

        if (neighbourhood.size()!=0) successes+=1;

        auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count(); // time passed
        cout<<duration/1000000<<"s"<<endl<<endl;

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

      outfile<<c<<","<<prop<<","<<mean_duration<<","<<mean_total_space<<","<<mean_reservoir_space<<","<<mean_degree_space<<","<<mean_edges_checked<<","<<variance_duration<<","<<variance_total_space<<","<<variance_reservoir_space<<","<<variance_degree_space<<","<<variance_edges_checked<<","<<successes<<endl; // write values to file

    }

  }
  outfile.close();
}

// perform reservoir sampling
// returns number of edges which are read
int single_pass_insertion_stream(int c, int d, int n,ifstream& stream, vector<vertex>& neighbourhood, vertex& root, double prop) {
  int size=ceil(prop*ceil(log10(n)*pow(n,(double)1/c)));
  int num_samplers=(2>log(n)/5) ? 2 : log(n)/5;
  num_samplers=c<num_samplers ? c : num_samplers;
  cout<<"num_samplers="<<num_samplers<<endl;

  // initalise edge & vector sets for each parallel run
  vector<vertex>* reservoirs[num_samplers]; vector<edge> edges;
  for (int i=0; i<num_samplers; i++) reservoirs[i]=new vector<vertex>;
  BYTES+=num_samplers*(sizeof(vector<vertex>))+sizeof(vector<edge>);
  RESERVOIR_BYTES+=num_samplers*(sizeof(vector<vertex>))+sizeof(vector<edge>);

  string line; edge e; map<vertex,int> degrees; // these are shared for each run
  map<vertex,bool*> num_reservoirs;
  int count[num_samplers];  // counts the number of vertexs >=d1
  BYTES+=sizeof(string)+sizeof(edge)+2*sizeof(map<vertex,int>)+sizeof(int);
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

    // update reservoirs
    for (int j=0; j<num_samplers; j++) { // perform parallel runs
      int d1=max(1,(j*d)/c), d2=d/c; // calculate degree bounds for run
      BYTES+=sizeof(int);
      // NB rest is standard degree-restricted sampling

      // Consider adding first vertex to the reservoir
      if (degrees[e.fst]==d1) {
        count[j]+=1; // increment number of d1 degree vertexs
        update_reservoir(e.fst,num_samplers,j,d1,d2,count[j],size,*reservoirs[j],edges,num_reservoirs); // possibly add value to reservoir
      }

      // Consider adding second vertex to the reservoir
      if (degrees[e.snd]==d1) {
        count[j]+=1; // increment number of d1 degree vertexs
        update_reservoir(e.snd,num_samplers,j,d1,d2,count[j],size,*reservoirs[j],edges,num_reservoirs); // possibly add value to reservoir
      }

      BYTES-=sizeof(int); // deletion of d1
    }

    // update edge set
    // TODO think (kind of seems like I am cheating here as it essentially adds the edge if it has an endpoint in any reservoir regardless of the bounds of that reservoir, but I'm pretty sure this is acceptable as the edge would always be added if it is in a reservoir)
    // Only variation here is that
    bool inserted=false;
    if (num_reservoirs.find(e.fst)!=num_reservoirs.end()) { // if first endpoint is in reservoir
      for (int j=0; j<num_samplers; j++) {
        // NOTE inserted==false always at this point
        int d1=max(1,(j*d)/c), d2=d/c; // calculate degree bounds for run
        if (num_reservoirs[e.fst][j]==true && degrees[e.fst]<=d2+d1) {
          edges.push_back(e);
          BYTES+=sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=sizeof(edge);
          inserted=true;
        }
        if (num_reservoirs[e.fst][j]==true && degrees[e.fst]>=d2+d1) { // sufficient neighbourhood has been found, return it
          cout<<endl<<"*"<<degrees[e.fst]<<"/"<<edges.size()<<endl;
          for (vector<edge>::iterator i=edges.begin(); i!=edges.end(); i++) { // construct neighbourhood to be returned
            if (i->fst==e.fst) neighbourhood.push_back(i->snd);
            else if (i->snd==e.fst) neighbourhood.push_back(i->fst);
          }
          BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=neighbourhood.size()*sizeof(edge);
          root=e.fst;
          return edge_count;
        }
        if (inserted) break; // don't check any more reservoirs if it has been inserted (termination won't succeed anyway & stops duplication)
      }
    } else if (num_reservoirs.find(e.snd)!=num_reservoirs.end()) { // if second endpoint is in a reservoir
      for (int j=0; j<num_samplers; j++) {
        // NOTE inserted==false always at this point
        int d1=max(1,(j*d)/c), d2=d/c; // calculate degree bounds for run
        if (num_reservoirs[e.snd][j]==true && degrees[e.snd]<=d2+d1) {
          edges.push_back(e);
          BYTES+=sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=sizeof(edge);
          inserted=true;
        }
        if (num_reservoirs[e.snd][j]==true && degrees[e.snd]>=d2+d1) { // sufficient neighbourhood has been found, return it
          cout<<endl<<"*"<<degrees[e.snd]<<"/"<<edges.size()<<endl;
          for (vector<edge>::iterator i=edges.begin(); i!=edges.end(); i++) { // construct neighbourhood to be returned
            if (i->fst==e.snd) neighbourhood.push_back(i->snd);
            else if (i->snd==e.snd) neighbourhood.push_back(i->fst);
          }
          BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
          RESERVOIR_BYTES+=neighbourhood.size()*sizeof(edge);
          root=e.snd;
          return edge_count;
        }
        if (inserted) break; // don't check any more reservoirs if it has been inserted (termination won't succeed anyway & stops duplication)
      }
    }

  }
  cout<<"\rDONE                         "<<endl;

  // No sucessful runs
  neighbourhood.clear();
  vertex* p=&root;
  p=nullptr;

  return edge_count;
}

void update_reservoir(vertex n,int c,int res_num, int d1, int d2, int count, int size, vector<vertex>& reservoir, vector<edge>& edges, map<vertex,bool*>& num_reservoirs) {
  if (reservoir.size()<size) { // reservoir is not full
    reservoir.push_back(n);
    BYTES+=sizeof(vertex); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
    RESERVOIR_BYTES+=2*sizeof(vertex);
    if (num_reservoirs.find(n)==num_reservoirs.end()) {
      //cout<<"NEW ENTRY";
      bool* arr_p=new bool[c];
      for (int j=0; j<c; j++) arr_p[j]=false;
      arr_p[res_num]=true;
      num_reservoirs[n]=arr_p;
      BYTES+=sizeof(vertex)+sizeof(bool*)+sizeof(bool)*c; if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
      RESERVOIR_BYTES+=sizeof(vertex)+sizeof(bool*)+sizeof(bool)*c;
    } else num_reservoirs[n][res_num]+=1;
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
      if (RESERVOIR_BYTES>MAX_RESERVOIR_BYTES) MAX_RESERVOIR_BYTES=RESERVOIR_BYTES;
      if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;

      // update reservoirs
      reservoir.erase(reservoir.begin()+to_delete);
      reservoir.push_back(n);

      // delete if empty
      num_reservoirs[to_delete_val][res_num]=false;
      int count=0;
      for (int j=0; j<c; j++) count+=int(num_reservoirs[to_delete_val][j]);
      if (count==0) {
        num_reservoirs.erase(to_delete_val);
        BYTES-=sizeof(vertex)+sizeof(bool*)+sizeof(bool)*c;
        RESERVOIR_BYTES-=sizeof(vertex)+sizeof(bool*)+sizeof(bool)*c;
      }

      // add new
      if (num_reservoirs.find(n)==num_reservoirs.end()) {
        bool* arr_p=new bool[c];
        for (int j=0; j<c; j++) arr_p[j]=false;
        arr_p[res_num]=true;
        num_reservoirs[n]=arr_p;
        BYTES+=sizeof(vertex)+sizeof(bool*)+sizeof(bool)*c; if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
        RESERVOIR_BYTES+=sizeof(vertex)+sizeof(bool*)+sizeof(bool)*c;
      } else num_reservoirs[n][res_num]=1;

      // NB No space change
      // removes edges adjacent to vertex to be deleted (which are not adjacent to another vertex in the reservoir)
      vector<edge>::iterator i=edges.begin();
      while (i!=edges.end()) {
        // count number of resevoirs endpoints appear in
        bool needed_fst=false, needed_snd=false;
        if (num_reservoirs.find(i->fst)!=num_reservoirs.end()) for (int j=0; j<c; j++) needed_fst=(needed_fst || num_reservoirs[i->fst][j]);
        if (num_reservoirs.find(i->snd)!=num_reservoirs.end()) for (int j=0; j<c; j++) needed_snd=(needed_snd || num_reservoirs[i->snd][j]);

        if (needed_fst==false && needed_snd==false) { // edge not adjacent to another value in any other reservoir
          i=edges.erase(i); // next edge
          BYTES-=sizeof(edge);RESERVOIR_BYTES-=sizeof(edge);
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
