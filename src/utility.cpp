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
#include <dirent.h>
#include <sys/types.h>

using namespace std;

/*-----------------*
 * DATA STRUCTURES *
 *-----------------*/

using vertex=string;

struct edge { // undirected edge
  vertex fst;
  vertex snd;
  bool insertion;
};


/*------------*
 * SIGNATURES *
 *------------*/

void relabel_vertices(string input_path, string output_path);
void parse_edge(string str, edge& e);
void generate_insertion_deletion(string file_name);
void list_vertices(string file_name);
void merge_files(string* file_names, string output_name);
void greatest_degree(string str, string& vertex, int& degree);
void count_final_edges(string file_name, int& count);
void count_total_edges(string file_name, int& count);
void merge_directory(const char* path, string output_name);

/*------*
 * BODY *
 *------*/

int main() {
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

void merge_directory(const char* path, string output_name) {
  DIR *dir = opendir(path);
  if (dir == NULL) return; // directory doesnt exist

  ofstream outfile(output_name); dirent *entry;
  int i=0;

  while ((entry = readdir(dir)) != NULL) { // while there are files in directory
    string file_name=entry->d_name;
    if (file_name.size()>5 && file_name.substr(file_name.size()-5)=="edges") { // is an edge file
      cout<<"STARTING <"<<i<<">-"<<file_name<<endl;
      string line;
      ifstream stream(path+file_name); // output file
      while (getline(stream,line)) outfile<<line<<endl; // write lines to end of outputfile
      stream.close();
      cout<<"ENDING-"<<file_name<<endl;
      i+=1;
    }
  }
  outfile.close();
  closedir(dir);
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
  int i=0;

  while (getline(stream,line)) {
    i+=1;
    if (i%10000==0) cout<<i<<", ";
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
  cout<<endl;
  stream.close();
  outfile.close();
}

// produces list of vertices & their degree
// works for insetion & insertion-deletion streams
void list_vertices(string file_name) {
  ifstream stream(file_name+".edges");
  ofstream outfile(file_name+".vertices");

  map<string,int> degrees; string line; edge e;
  int i=0;

  while (getline(stream,line)) {
    i+=1;
    if (i%100000==0) cout<<"\r"<<i; // update on how many edges have been checked
    parse_edge(line,e);

    if (degrees.count(e.fst)) { // vertex already in graph
      if (e.insertion) degrees[e.fst]+=1; // insertion edge
      else degrees[e.fst]-=1; // deletion edge
    } else { // vertex not in graph
      degrees[e.fst]=1; // cannot delete an edge which is not in the graph
    }

    if (degrees.count(e.snd)) { // vertex already in graph
      if (e.insertion) degrees[e.snd]+=1; // insertion edge
      else degrees[e.snd]-=1; // deletion edge
    } else { // vertex not in graph
      degrees[e.snd]=1; // cannot delete an edge which is not in the graph
    }
  }
  cout<<endl;
  // write to vertices file
  for (map<string,int>::iterator i=degrees.begin(); i!=degrees.end(); i++) {
    outfile<<i->first<<","<<i->second<<endl; // name,degree
  }
}

// relabel vertices of ID stream so they are all in [0,n)
void relabel_vertices(string input_path, string output_path) {
  ifstream stream(input_path);
  ofstream outfile(output_path);

  map<vertex,int> new_labels; string line; edge e;
  map<vertex,int>::iterator it;
  int count=0, i=0;

  while (getline(stream,line)) {
    i++;
    if (i%100000==0) cout<<"\r"<<i; // update on how many edges have been checked
    parse_edge(line,e);

    int new_fst, new_snd;

    // get new label for fst
    it=new_labels.find(e.fst);
    if (it==new_labels.end()) { // not currently mapped
      new_labels[e.fst]=count;
      new_fst=count;
      count++;
    } else new_fst=it->second;

    // get new label for snd
    it=new_labels.find(e.snd);
    if (it==new_labels.end()) { // not currently mapped
      new_labels[e.snd]=count;
      new_snd=count;
      count++;
    } else new_snd=it->second;

    outfile<<line[0]<<" "<<new_fst<<" "<<new_snd<<endl;
  }
  cout<<"\rDONE                 ";
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
  int i=0;

  while (getline(stream,line)) {
    i+=1;
    if (i%100000==0) cout<<i<<",";
    parse_edge(line,e);

    if (degrees.count(e.fst)) {
      if (e.insertion) degrees[e.fst]+=1;
      else degrees[e.fst]-=1;
    } else {
      // You cannot delete an edge which does not exist
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
void count_final_edges(string file_name, int& count) {
  ifstream stream(file_name);
  string line; count=0; edge e;
  while (getline(stream,line)) {
    parse_edge(line,e);
    if (e.insertion) count+=1;
    else count-=1;
  }
}

void count_total_edges(string file_name, int& count) {
  ifstream stream(file_name);
  count=0; string line;
  while (getline(stream,line)) count++;
}
