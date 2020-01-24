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
void deg_res_sampling(int d1, int d2, int size, ifstream& stream, vector<int> &neighbourhood);
void update_resevoir(int node, int d1, int d2, int count, int size, vector<int>& resevoir, vector<edge>& edges);
void parse_edge(string str, edge& e);

/*-----*
* BODY *
*------*/

int main() {
  ifstream stream("data/facebook.edges");
  vector<int> n;
  deg_res_sampling(1,20,2,stream,n); // NB does not tell you whose neighbourhood it is

  for (vector<int>::iterator i=n.begin(); i!=n.end(); i++) cout<<*i<<endl;

  return 0;
}

// perform resevoir sampling
void deg_res_sampling(int d1, int d2, int size, ifstream& stream, vector<int> &neighbourhood) {
  string line; edge e; map<int,int> degrees; vector<edge> edges; vector<int> resevoir;
  int count;  // number of nodes >=d1

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

    // Consider adding first node to the resevoir
    if (degrees[e.fst]==d1) {
      count+=1; // increment number of d1 degree nodes
      update_resevoir(e.fst,d1,d2,count,size,resevoir,edges); // possibly add value to resevoir
    }

    // Consider adding second node to the resevoir
    if (degrees[e.snd]==d1) {
      count+=1; // increment number of d1 degree nodes
      update_resevoir(e.snd,d1,d2,count,size,resevoir,edges); // possibly add value to resevoir
    }

    // If one the endpoints is in the resevoir add the edge to collection of edges
    if ((find(resevoir.begin(),resevoir.end(),e.fst)!=resevoir.end() && degrees[e.fst]<d2+d1)
      || (find(resevoir.begin(),resevoir.end(),e.snd)!=resevoir.end() && degrees[e.snd]<d2+d1)) {
      edges.push_back(e);
    }

  }

  while (true) {
    if (resevoir.size()==0) { // unsuccessful
      neighbourhood.clear();
      return;
    }

    int x=rand()%resevoir.size(); // randomly choose a node
    int node=resevoir[x];
    if (degrees[node]>=d2+d1) { // check it has a sufficiently large neighbourhood
      for (vector<edge>::iterator i=edges.begin(); i!=edges.end(); i++) { // construct & return the neighbourhood
        if (i->fst==node) neighbourhood.push_back(i->snd);
        if (i->snd==node) neighbourhood.push_back(i->fst);
      }
      return;
    } else { // pick a different node
      resevoir.erase(resevoir.begin()+x);
    }
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
