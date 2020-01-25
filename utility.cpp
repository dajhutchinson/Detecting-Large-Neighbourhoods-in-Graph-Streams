/*
 *  Takes an insertion only stream and creates an insertion-deletion stream
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>

using namespace std;

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

struct edge { // undirected edge
  int fst;
  int snd;
};


/*------------*
 * SIGNATURES *
 *------------*/

void parse_edge(string str, edge& e);

/*------*
 * BODY *
 *------*/

int main() {
  string file_name="data/facebook";
  ofstream outfile(file_name+"_deletion.edges");
  ifstream stream(file_name+".edges");

  vector<edge> edges; // edges in system
  string line; edge e;

  float p=0.1; // Probability to delete an edge // VARY THIS
  default_random_engine generator;
  bernoulli_distribution d(p);

  while (getline(stream,line)) {
    parse_edge(line,e);
    edges.push_back(e);
    outfile<<"I "<<e.fst<<" "<<e.snd<<endl;
    cout<<"I "<<e.fst<<" "<<e.snd<<endl;
    while (d(generator) && edges.size()>0) {
      int to_delete=rand()%edges.size();
      edge edge_to_delete=edges[to_delete];
      outfile<<"D "<<edge_to_delete.fst<<" "<<edge_to_delete.snd<<endl;
      cout<<"D "<<edge_to_delete.fst<<" "<<edge_to_delete.snd<<endl;
      edges.erase(edges.begin()+to_delete);
    }
  }

  stream.close();
  outfile.close();
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
