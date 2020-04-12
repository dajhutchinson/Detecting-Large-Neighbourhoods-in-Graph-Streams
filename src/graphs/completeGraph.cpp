/*-----------------*
 * Creates an insertion graph stream for a complete graph with a given number of nodes
 *
 * Complete graphs contain an edge between each vertex
 * In a Star Graph of order N, every node has degree N-1
 *
 * Creates
 *    _.edges = Insertion edge stream (formatted as "v1 v2")
 *    _.vertices = list of vertices and their degrees (formatted as CSV file)
 *-----------------*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <chrono>

using namespace std;

// signatures
void generate_complete_insertion_only(string file_name, int num_vertices);
void generate_complete_insertion_deletion(string file_name, int num_vertices, double proportion);

int main(int argc, char* argv[]) {
  cout<<argc<<endl;
  if (argc!=4 && argc!=5) {
    cout<<"ERROR: ./cg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;
  } else if (argc==4) {
    string file_name=argv[1];
    string type=argv[2];
    int num_vertices=atoi(argv[3]); // Graph degree

    if (type=="I") generate_complete_insertion_only(file_name,num_vertices);
    else if (type=="ID") generate_complete_insertion_deletion(file_name,num_vertices,.01);
    else cout<<"INVALID TYPE - ./cg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;
  } else if (argc==4) {
    string file_name=argv[1];
    string type=argv[2];
    int num_vertices=atoi(argv[3]); // Graph degree
    double p=stod(argv[4]); // Graph degree

    if (p<0 || p>=1) cout<<"INVALID PROBABILITY - ./cg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;

    if (type=="ID") generate_complete_insertion_deletion(file_name,num_vertices,p);
    else cout<<"INVALID TYPE - ./cg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;
  }
  return 0;
}

/*-------*
 *  BODY *
 *-------*/

void generate_complete_insertion_only(string file_name, int num_vertices) {
  // Prepare files for output
  ofstream edge_file(file_name+".edges");
  ofstream vertex_file(file_name+".vertices");

  for (int v=1; v<num_vertices+1; v++) {
    vertex_file<<v<<","<<num_vertices-1<<endl;
    for (int u=v+1; u<num_vertices+1; u++) edge_file<<v<<" "<<u<<endl;
  }

  edge_file.close();
  vertex_file.close();
}

void generate_complete_insertion_deletion(string file_name, int num_vertices, double proportion) {
  ofstream edge_file(file_name+".edges");
  ofstream vertex_file(file_name+".vertices");

  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time
  bernoulli_distribution bernoulli_d(proportion);

  for (int v=1; v<num_vertices+1; v++) {
    vertex_file<<v<<","<<num_vertices-1<<endl;
    for (int u=v+1; u<num_vertices+1; u++) {
      edge_file<<"I "<<v<<" "<<u<<endl;
      if (bernoulli_d(generator)) { // delete & insert
        edge_file<<"D "<<v<<" "<<u<<endl;
        edge_file<<"I "<<v<<" "<<u<<endl;
      }
    }
  }

  edge_file.close();
  vertex_file.close();
}
