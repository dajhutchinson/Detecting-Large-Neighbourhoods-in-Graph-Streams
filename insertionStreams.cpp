/*-----------------*
 * TODO
 * Record time take & space used
 *-----------------*/

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace std;

int BYTES=0; // space used atm
int MAX_BYTES=0; // max space used at any time

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

using vertex = string; // typemap vertex

struct edge { // undirected edge
  vertex fst;
  vertex snd;
};

typedef chrono::high_resolution_clock::time_point time_point;

/*-----------*
* SIGNATURES *
*------------*/

// size = size of resevoir
void single_pass_insertion_stream(int c, int d, int size, ifstream& stream, vector<vertex>& neighbourhood, vertex& root);
void update_resevoir(vertex n, int d1, int d2, int count, int size, vector<vertex>& resevoir, vector<edge>& edges);
void parse_edge(string str, edge& e);

/*-----*
* BODY *
*------*/

int main() {
  int c=2, d=40, s=2; // c=runs, d/c=d2, s=size
  ifstream stream("data/gplus.edges"); // file to read
  vector<vertex> neighbourhood; vertex root; // variables for returned values

  time_point before=chrono::high_resolution_clock::now(); // time before execution
  single_pass_insertion_stream(c,d,s,stream,neighbourhood,root); // NB does not tell you whose neighbourhood it is
  time_point after=chrono::high_resolution_clock::now(); // time after execution

  // Print out returned neighbourhood, if one exists
  vertex* p=&root;
  if (p==nullptr) cout<<"NO SUCCESSES"<<endl;
  else {
    cout<<"Neighbourhood for <"<<root<<">"<<endl;
    cout<<"<";
    for (vector<vertex>::iterator i=neighbourhood.begin(); i!=neighbourhood.end(); i++) cout<<*i<<",";
    cout<<">"<<endl;
  }
  stream.close();

  auto duration = chrono::duration_cast<chrono::microseconds>(after-before).count();
  cout<<"Execution Time - "<<duration<<" microseconds ("<<duration/1000<<" milliseconds)"<<endl;
  cout<<"Max space - "<<MAX_BYTES<<" bytes ("<<MAX_BYTES/1024<<" kb)"<<endl;
  return 0;
}

// perform resevoir sampling
void single_pass_insertion_stream(int c, int d, int size, ifstream& stream, vector<vertex>& neighbourhood, vertex& root) {

  // initalise edge & vector sets for each parallel run
  vector<vertex>* resevoirs[c]; vector<edge>* edges[c];
  for (int i=0; i<c; i++) {
    resevoirs[i]=new vector<vertex>;
    edges[i]=new vector<edge>;
  }
  BYTES+=c*(sizeof(vector<vertex>)+sizeof(vector<edge>));

  string line; edge e; map<vertex,int> degrees; // these are shared for each run
  int count[c];  // counts the number of vertexs >=d1
  BYTES+=sizeof(string)+sizeof(edge)+sizeof(map<vertex,int>)+sizeof(int);

  while (getline(stream,line)) { // While stream is not empty
    parse_edge(line,e);

    // increment degrees for each vertex
    if (degrees.count(e.fst)) {
      degrees[e.fst]+=1;
    } else {
      degrees[e.fst]=1;
      BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*); // sizeof(void*)=size of pointer
    }

    if (degrees.count(e.snd)) {
      degrees[e.snd]+=1;
    } else {
      degrees[e.snd]=1;
      BYTES+=sizeof(vertex)+sizeof(int)+sizeof(void*); // sizeof(void*)=size of pointer
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

      // If one the endpoints is in the resevoir add the edge to collection of edges
      if ((find(resevoirs[j]->begin(),resevoirs[j]->end(),e.fst)!=resevoirs[j]->end() && degrees[e.fst]<d2+d1)
        || (find(resevoirs[j]->begin(),resevoirs[j]->end(),e.snd)!=resevoirs[j]->end() && degrees[e.snd]<d2+d1)) {
        edges[j]->push_back(e);
        BYTES+=sizeof(edge);
      }

      if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
      BYTES-=sizeof(int);
    }

  }

  // find successful neighbourhoods for each run
  vector<tuple<int,vertex> > success; // (run, vertex_id)
  for (int j=0; j<c; j++) { // cheack each run
    int d1=max(1,(j*d)/c), d2=d/c; // calculate required degree bounds
    for (vector<vertex>::iterator i=resevoirs[j]->begin(); i!=resevoirs[j]->end(); i++) {
      if (degrees[*i]>=d1+d2) success.push_back(make_tuple(j,*i)); // record which vertexs had sufficient degree in each run
    }
  }
  BYTES+=sizeof(vector<tuple<int,vertex> >)+success.size()*sizeof(tuple<int,vertex>);
  if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;

  if (success.size()==0) { // No sucessful runs
    neighbourhood.clear();
    vertex* p=&root;
    p=nullptr;
    return;
  } else { // at least one successful run
    int x=rand()%success.size(); // randomly choose a successful (run, vertex_id) pair
    BYTES+=sizeof(int);
    root=get<1>(success[x]);
    vector<edge>* edge_set=edges[get<0>(success[x])];
    for (vector<edge>::iterator i=edge_set->begin(); i!=edge_set->end(); i++) { // construct neighbourhood to be returned
      if (i->fst==root) neighbourhood.push_back(i->snd);
      else if (i->snd==root) neighbourhood.push_back(i->fst);
    }
    BYTES+=neighbourhood.size()*sizeof(edge); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
    BYTES-=sizeof(int);
    return;
  }

}

void update_resevoir(vertex n, int d1, int d2, int count, int size, vector<vertex>& resevoir, vector<edge>& edges) {
  if (resevoir.size()<size) { // resevoir is not full
    resevoir.push_back(n);
    BYTES+=sizeof(vertex); if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;
  } else { // resevoir is full
    default_random_engine generator;
    bernoulli_distribution d((float)size/(float)count);
    BYTES+=sizeof(default_random_engine)+sizeof(bernoulli_distribution);
    if (d(generator)) { // if coin flip passes
      int to_delete=rand()%size;
      vertex to_delete_val=resevoir[to_delete];
      BYTES+=2*sizeof(vertex);
      if (BYTES>MAX_BYTES) MAX_BYTES=BYTES;

      resevoir.erase(resevoir.begin()+to_delete);
      resevoir.push_back(n);
      // NB No space change

      vector<edge>::iterator i=edges.begin();
      while (i!=edges.end()) {
        if (i->fst==to_delete_val || i->snd==to_delete_val) {
          i=edges.erase(i);
          BYTES-=sizeof(edge);
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