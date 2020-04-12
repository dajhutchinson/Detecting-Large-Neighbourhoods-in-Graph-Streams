/*-----------------*
 * Creates an insertion graph stream for a star graph with a given number of nodes
 *
 * Star graphs have one node which is connected to every other node (the other nodes are only conneted to this node)
 * In a Star Graph of order N, one node has degree N-1 and the remaining 1 have degree 1.
 *
 * Creates
 *    _.edges = Insertion edge stream (formatted as "v1 v2")
 *    _.vertices = list of vertices and their degrees (formatted as CSV file)
 *-----------------*/

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace std;

// signatures
void generate_star_insertion_only(string file_name, int num_vertices);
void generate_star_insertion_deletion(string file_name, int num_vertices, double proportion);

/*-------*
 *  BODY *
 *-------*/

int main(int argc, char* argv[]) {

  if (argc!=4 && argc!=5) {
    cout<<"ERROR: ./sg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;
  } else if (argc==4) {
    string file_name=argv[1];
    string type=argv[2];
    int num_vertices=atoi(argv[3]); // Graph degree

    if (type=="ID") generate_star_insertion_deletion(file_name,num_vertices,0.01);
    else if (type=="I") generate_star_insertion_only(file_name,num_vertices);
    else cout<<"INVALID TYPE - ./sg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;

  } else if (argc==5) {
    string file_name=argv[1];
    int num_vertices=atoi(argv[2]); // Graph degree
    string type=argv[3];
    double p=stod(argv[4]);

    if (p<0 || p>=1) cout<<"INVALID PROBABILITY - ./cg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;

    if (type=="ID") generate_star_insertion_deletion(file_name,num_vertices,p);
    else cout<<"INVALID TYPE - ./sg [file_name] [I/ID] [graph_order] (deletion_probabilty)"<<endl;
  }
  return 0;
}

/*-------*
*  BODY *
*-------*/

void generate_star_insertion_only(string file_name, int num_vertices) {
  // Prepare files for output
  ofstream edge_file(file_name+".edges");
  ofstream vertex_file(file_name+".vertices");

  // Choose node to be centre (Uniformly)
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed generator with current time
  uniform_int_distribution<int> uniform_d(1,num_vertices);
  int centre=uniform_d(generator);
  cout<<"Centre:"<<centre<<endl;

  // create edge between centre and every other node
  for (int v=1; v<num_vertices+1; v++) {
    if (v!=centre) {
      edge_file<<centre<<" "<<v<<endl;
      vertex_file<<v<<",1"<<endl;
    } else { // no edge
      vertex_file<<v<<","<<num_vertices-1<<endl;
    }
  }

  edge_file.close();
  vertex_file.close();
}

void generate_star_insertion_deletion(string file_name, int num_vertices, double proportion) {
  // Prepare files for output
  ofstream edge_file(file_name+".edges");
  ofstream vertex_file(file_name+".vertices");

  // Choose node to be centre (Uniformly)
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed generator with current time
  uniform_int_distribution<int> uniform_d(1,num_vertices);
  bernoulli_distribution bernoulli_d(proportion);

  int centre=uniform_d(generator);
  cout<<"Centre:"<<centre<<endl;

  // create edge between centre and every other node
  for (int v=1; v<num_vertices+1; v++) {
    if (v!=centre) {
      edge_file<<"I "<<centre<<" "<<v<<endl;
      vertex_file<<v<<",1"<<endl;

      if (bernoulli_d(generator)) {
        edge_file<<"D "<<centre<<" "<<v<<endl;
        edge_file<<"I "<<centre<<" "<<v<<endl;
      }

    } else { // no edge
      vertex_file<<v<<","<<num_vertices-1<<endl;
    }
  }

  edge_file.close();
  vertex_file.close();
}
