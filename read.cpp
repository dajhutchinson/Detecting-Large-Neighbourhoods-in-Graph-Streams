/*
 *  Read & print each edge from a file
 */

#include <iostream>
#include <fstream>
#include <string>

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
  ifstream stream("data/facebook_small.edges"); // File to read from
  string line; // current line being read
  string fst,snd;
  edge e;

  while (getline(stream,line)) { // Read each line
    parse_edge(line,e);
    cout<<e.fst<<" & "<<e.snd<<endl; // print
  }

  return 0;
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
