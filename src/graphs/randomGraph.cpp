/*
 * Creates an insertion graph stream for a graph with a defined number of edges & vertices
 *
 * Models creating the upper right triangle of a matrix which represents the edges in the graph
 *
 * Creates
 *    _.edges = Insertion edge stream (formatted as "v1 v2")
 *    _.vertices = list of vertices and their degrees (formatted as CSV file)
 */

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace std;

/*-------*
 *  BODY *
 *-------*/

int main(int argc, char* argv[]) {

  if (argc!=3 && argc!=4) {
    cout<<"ERROR: ./sg [file_name] [graph_order] (prob_edge)"<<endl;
  } else {
    string file_name=argv[1];
    int N=atoi(argv[2]); // Graph degree

    // Prepare files for output
    ofstream edge_file(file_name+".edges");
    ofstream vertex_file(file_name+".vertices");

    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time

    map<int,int> degrees;
    for (int i=1; i<N; i++) degrees[i]=0; // record degree of every node

    // bernoulli_d(generator);
    for (int i=1; i<N+1; i++) {
      double P;
      if (argc==4) { // specified probability
        P=atof(argv[3]); // Probs of making an edge
      } else { // randomly choose probability
        uniform_int_distribution<int> uniform_d(1,100);
        P=(double)uniform_d(generator)/100;
      }
      bernoulli_distribution bernoulli_d(P);

      for (int j=i; j<N+1; j++) {
        if (j!=i && bernoulli_d(generator)==1) {
          degrees[i]+=1; degrees[j]+=1; // update degrees
          edge_file<<i<<" "<<j<<endl;
      }}
      vertex_file<<i<<","<<degrees[i]<<endl;
    }

    edge_file.close();
    vertex_file.close();
  }

  return 0;
}
