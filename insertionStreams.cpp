#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>

using namespace std;

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

struct edge { // undirected edge
  int fst;
  int snd;
};

/*-----------*
* SIGNATURES *
*------------*/

// size = size of resevoir
void single_pass_insertion_stream(int c, int d, int size, ifstream& stream, vector<int> &neighbourhood);
void update_resevoir(int node, int d1, int d2, int count, int size, vector<int>& resevoir, vector<edge>& edges);
void parse_edge(string str, edge& e);

/*-----*
* BODY *
*------*/

int main() {
  int c=10, d=400, s=2; // c=runs, d/c=d2, s=size
  ifstream stream("data/facebook.edges");
  vector<int> n;
  single_pass_insertion_stream(c,d,s,stream,n); // NB does not tell you whose neighbourhood it is

  // Print out returned neighbourhood, if one exists
  if (n.size()==0) {
    cout<<"NO SUCCESSES"<<endl;
  } else {
    for (vector<int>::iterator i=n.begin(); i!=n.end(); i++) cout<<*i<<endl;
  }
  stream.close();
  return 0;
}

// perform resevoir sampling
void single_pass_insertion_stream(int c, int d, int size, ifstream& stream, vector<int> &neighbourhood) {

  // initalise edge & vector sets for each parallel run
  vector<int>* resevoirs[c]; vector<edge>* edges[c];
  for (int i=0; i<c; i++) {
    resevoirs[i]=new vector<int>;
    edges[i]=new vector<edge>;
  }

  string line; edge e; map<int,int> degrees; // these are shared for each run
  int count[c];  // number of nodes >=d1

  while (getline(stream,line)) { // While stream is not empty
    parse_edge(line,e);

    // increment degrees for each node
    if (degrees.count(e.fst)) {
      degrees[e.fst]+=1;
    } else {
      degrees[e.fst]=1;
    }

    if (degrees.count(e.snd)) {
      degrees[e.snd]+=1;
    } else {
      degrees[e.snd]=1;
    }

    for (int j=0; j<c; j++) { // perform parallel runs
      int d1=max(1,(j*d)/c), d2=d/c; // calculate degree bounds for run
      // NB rest is standard degree-restricted sampling

      // Consider adding first node to the resevoir
      if (degrees[e.fst]==d1) {
        count[j]+=1; // increment number of d1 degree nodes
        update_resevoir(e.fst,d1,d2,count[j],size,*resevoirs[j],*edges[j]); // possibly add value to resevoir
      }

      // Consider adding second node to the resevoir
      if (degrees[e.snd]==d1) {
        count[j]+=1; // increment number of d1 degree nodes
        update_resevoir(e.snd,d1,d2,count[j],size,*resevoirs[j],*edges[j]); // possibly add value to resevoir
      }

      // If one the endpoints is in the resevoir add the edge to collection of edges
      if ((find(resevoirs[j]->begin(),resevoirs[j]->end(),e.fst)!=resevoirs[j]->end() && degrees[e.fst]<d2+d1)
        || (find(resevoirs[j]->begin(),resevoirs[j]->end(),e.snd)!=resevoirs[j]->end() && degrees[e.snd]<d2+d1)) {
        edges[j]->push_back(e);
      }

    }

  }

  // find successful neighbourhoods for each run
  vector<tuple<int, int> > success; // (run, node_id)
  for (int j=0; j<c; j++) { // cheack each run
    int d1=max(1,(j*d)/c), d2=d/c; // calculate required degree bounds
    for (vector<int>::iterator i=resevoirs[j]->begin(); i!=resevoirs[j]->end(); i++) {
      if (degrees[*i]>=d1+d2) success.push_back(make_tuple(j,*i)); // record which nodes had sufficient degree in each run
    }
  }

  if (success.size()==0) { // No sucessful runs
    neighbourhood.clear();
    return
  } else { // at least one successful run
    int x=rand()%success.size(); // randomly choose a successful (run, node_id) pair
    int node=get<1>(success[x]);
    vector<edge>* edge_set=edges[get<0>(success[x])];
    for (vector<edge>::iterator i=edge_set->begin(); i!=edge_set->end(); i++) { // construct neighbourhood to be returned
      if (i->fst==node) neighbourhood.push_back(i->snd);
      else if (i->snd==node) neighbourhood.push_back(i->fst);
    }
    return
  }

}

void update_resevoir(int node, int d1, int d2, int count, int size, vector<int>& resevoir, vector<edge>& edges) {
  if (resevoir.size()<size) { // resevoir is not full
    resevoir.push_back(node);
  } else { // resevoir is full
    default_random_engine generator;
    bernoulli_distribution d((float)size/(float)count);
    if (d(generator)) { // if coin flip passes
      int to_delete=rand()%size;
      int to_delete_val=resevoir[to_delete];
      resevoir.erase(resevoir.begin()+to_delete);
      resevoir.push_back(node);

      vector<edge>::iterator i=edges.begin();
      while (i!=edges.end()) {
        if (i->fst==to_delete_val || i->snd==to_delete_val) i=edges.erase(i);
        else i++;
      }
    }
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
  e.fst=stoi(fst);
  e.snd=stoi(snd);
}
