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


    for (int v=1; v<N+1; v++) {
      vertex_file<<v<<","<<N-1<<endl;
      for (int u=v+1; u<N+1; u++) edge_file<<v<<" "<<u<<endl;
    }

    edge_file.close();
    vertex_file.close();

  }

  return 0;
}
