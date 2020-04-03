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

/*-------*
 *  BODY *
 *-------*/

int main(int argc, char* argv[]) {

  if (argc!=3) {
    cout<<"ERROR: ./sg [file_name] [graph_order]"<<endl;
  } else {
    string file_name=argv[1];
    int N=atoi(argv[2]); // Graph degree

    // Prepare files for output
    ofstream edge_file(file_name+".edges");
    ofstream vertex_file(file_name+".vertices");

    // Choose node to be centre (Uniformly)
    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed generator with current time
    uniform_int_distribution<int> uniform_d(1,N);
    int centre=uniform_d(generator);
    cout<<"Centre:"<<centre<<endl;

    // create edge between centre and every other node
    for (int v=1; v<N+1; v++) {
      if (v!=centre) {
        edge_file<<centre<<" "<<v<<endl;
        vertex_file<<v<<",1"<<endl;
      } else { // no edge
        vertex_file<<v<<","<<N-1<<endl;
      }
    }

    edge_file.close();
    vertex_file.close();

  }
  return 0;
}
