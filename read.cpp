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
  edge e;

  while (getline(stream,line)) { // Read each line
    parse_edge(line,e);
    cout<<e.fst<<" & "<<e.snd<<endl; // print
  }

  stream.close();

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

void parse_insertion_deletion_edge(string str, edge& e) {
  string fst="",snd="";
  int count=0;
  bool insert=NULL; // True for insertion, false for deletion

  for (char& c:str) {
    if (c==' ')count+=1;
    else if (count==0 && c=='I') insert=true;
    else if (count==0 && c=='D') insert=false;
    else if (count==1) fst+=c;
    else if (count==2) snd+=c;
  }

  e.fst=stoi(fst);
  e.snd=stoi(snd);
}
