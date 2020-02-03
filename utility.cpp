/*
 *  Takes an insertion only stream and creates an insertion-deletion stream
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <random>
#include <set>
#include <map>

using namespace std;

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

struct edge { // undirected edge
  string fst;
  string snd;
  bool insertion;
};


/*------------*
 * SIGNATURES *
 *------------*/

void parse_edge(string str, edge& e);
void generate_insertion_deletion(string file_name);
void list_vertices(string file_name);
void merge_files(string* file_names, string output_name);
void greatest_degree(string str, string& vertex, int& degree);
void count_edges(string file_name, int& count);

/*------*
 * BODY *
 *------*/

int main() {
  /*
  generate_insertion_deletion("data/gplus");
  cout<<"Deletions created"<<endl;
  string greatest; int deg;
  greatest_degree("data/gplus_deletion.edges",greatest,deg);
  cout<<greatest<<" "<<deg<<endl;

  list_vertices("data/gplus_deletion");
  int count;
  count_edges("data/gplus_deletion.edges",count);
  cout<<count<<endl;
  */
  return 0;
 }

// merge files into one by appending to the end of each other
void merge_files(string* file_names, string output_name) {
  ofstream outfile(output_name);
  string line;
  for (int i=0; i<sizeof(file_names); i++) { // each input file
    ifstream stream(file_names[i]); // output file
    while (getline(stream,line)) outfile<<line<<endl; // write lines to end of outputfile
    stream.close();
  }
  outfile.close();
}

// Takes an insertion only stream & generates an insertion-deletion stream
void generate_insertion_deletion(string file_name) {
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

    while (d(generator) && edges.size()>0) {
      int to_delete=rand()%edges.size();
      edge edge_to_delete=edges[to_delete];
      outfile<<"D "<<edge_to_delete.fst<<" "<<edge_to_delete.snd<<endl;
      edges.erase(edges.begin()+to_delete);
    }
  }

  stream.close();
  outfile.close();
}

void list_vertices(string file_name) {
  ofstream outfile(file_name+".vertices");
  ifstream stream(file_name+".edges");

  set<string> vertices; string line; edge e;

  while(getline(stream,line)) {
    parse_edge(line,e);

    if (vertices.find(e.fst)==vertices.end()) {
      vertices.insert(e.fst);
      outfile<<e.fst<<endl;
    }

    if (vertices.find(e.snd)==vertices.end()) {
      vertices.insert(e.snd);
      outfile<<e.snd<<endl;
    }
  }
}

// parse ege from stream
void parse_edge(string str, edge& e) {
  int spaces=count(str.begin(),str.end(),' ');

  if (spaces==2) { // insertion deletion edge
    if (str[0]=='D') e.insertion=false;
    else e.insertion=true;
    str=str.substr(2,str.size());
  } else if (spaces==1) { // insertion edge
    e.insertion=true;
  }

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
  e.fst=fst;
  e.snd=snd;
}

// return list of vertices with greatest degree & the degree
void greatest_degree(string file_name, string& vertex, int& degree) {
  degree=0; vertex=""; // initalise

  ifstream stream(file_name);
  map<string,int> degrees; string line; edge e;

  while (getline(stream,line)) {
    parse_edge(line,e);

    if (degrees.count(e.fst)) {
      if (e.insertion) degrees[e.fst]+=1;
      else degrees[e.fst]-=1;
    } else {
      degrees[e.fst]=1;
    }

    if (degrees.count(e.snd)) {
      if (e.insertion) degrees[e.snd]+=1;
      else degrees[e.snd]-=1;
    } else {
      degrees[e.snd]=1;
    }
  }

  for (map<string,int>::iterator i=degrees.begin(); i!=degrees.end(); i++) {
    if (i->second>degree) {
      degree=i->second;
      vertex=i->first;
    } else if (i->second==degree) {
      vertex.append(",").append(i->first);
    }
  }
}

// return number of edges in graph
void count_edges(string file_name, int& count) {
  ifstream stream(file_name);
  string line; count=0; edge e;
  while (getline(stream,line)) {
    parse_edge(line,e);
    if (e.insertion) count+=1;
    else count-=1;
  }
}