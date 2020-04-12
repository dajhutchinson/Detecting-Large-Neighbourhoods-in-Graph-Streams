/*
 * Creates insertion-only or insertion-deletion graph streams
 *
 * USING - ./rg [file_name] [I/ID] [graph_order] (prob_edge)
 *  file_name   = name of .edges & .vertices file to output (can be a path)
 *  I/ID        = Insertion-only or Insertion-Deletion
 *  graph_order = number of vertices in graph
 *  prob_edge   = probability of an edge (insertion-only & optional)
 *
 * INSERTION ONLY
 *  If prob_edge is _not_ specified then for each vertex U[0,1] is sampled to define its probability of having an edge to another arbitrary vertex
 *  Otherwise, all vertices have the same probability (prob_edge).
 *  This probability is used in a Bernoulli distribution to decide whether to add an edge.
 *
 * INSERTION-DELETION
 *  Each vertex samples from Poi(0.8), this value is then used as mean for another Poisson distribution.
 *  The second distribution is sampled to decide the number of edges between this vertex and another one.
 *  Edges added as Insertion then deletion until the required number of edges is produced
 *
 * Returns
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

/*-------------*
 *  SIGNATURES *
 *-------------*/

// p=NULL to randomly generate
void insertion_only(ofstream& edge_file, ofstream& vertex_file, int num_vertices, double p, default_random_engine generator);
void insertion_deletion(ofstream& edge_file, ofstream& vertex_file, int num_vertices, default_random_engine generator);
/*-------*
 *  BODY *
 *-------*/

int main(int argc, char* argv[]) {

  if (argc!=4 && argc!=5) {
    cout<<"ERROR: rg [file_name] [I/ID] [graph_order] (prob_edge)"<<endl; // I=Insertion, ID=Insertion-Deletion
  } else {
    string file_name=argv[1];
    string type=argv[2];
    int num_vertices=atoi(argv[3]); // Graph degree

    double p=-1.0;
    if (argc==5) {
      p=stod(argv[4]); // specified
      if (p<0 || p>1) {
        cout<<"ERROR: prob_edge must be in [0,1]\n./rg [file_name] [I/ID] [graph_order] (prob_edge)"<<endl;
        return -1;
      }
    }

    // Prepare files for output
    if (type=="I" || type=="ID") {
      ofstream edge_file(file_name+".edges");
      ofstream vertex_file(file_name+".vertices");

      default_random_engine generator;
      generator.seed(chrono::system_clock::now().time_since_epoch().count()); // seed with current time

      if (type=="I") insertion_only(edge_file,vertex_file,num_vertices,p,generator);
      else if (type=="ID") insertion_deletion(edge_file,vertex_file,num_vertices,generator);

      edge_file.close();
      vertex_file.close();
    } else { // type not recognised
      cout<<"ERROR: must state type of stream 'I' or 'ID'.\n./rg [file_name] [I/ID] [graph_order] (prob_edge)"<<endl;
      return -1;
    }
  }

  return 0;
}

/*
 *
 */

// generates insertion-only graph stream
// p==-1 if we want to randomly generate for each vertex
void insertion_only(ofstream& edge_file, ofstream& vertex_file, int num_vertices, double p, default_random_engine generator) {

  map<int,int> degrees;
  for (int i=1; i<num_vertices; i++) degrees[i]=0; // record degree of every node

  for (int i=1; i<num_vertices+1; i++) {
    double prob;
    if (p!=-1) { // specified probability
      prob=p;
    } else { // randomly choose probability
      uniform_int_distribution<int> uniform_d(1,100);
      prob=(double)uniform_d(generator)/100;
    }
    bernoulli_distribution bernoulli_d(prob);

    for (int j=i+1; j<num_vertices+1; j++) {
      if (bernoulli_d(generator)==1) {
        degrees[i]+=1; degrees[j]+=1; // update degrees
        edge_file<<i<<" "<<j<<endl;
    }}
    vertex_file<<i<<","<<degrees[i]<<endl;
  }

}

// Generates an insertion_deletion stream
void insertion_deletion(ofstream& edge_file, ofstream& vertex_file, int num_vertices, default_random_engine generator) {

  map<int,int> degrees;
  for (int i=1; i<num_vertices; i++) degrees[i]=0; // record degree of every node

  poisson_distribution<> poisson_d(.5); // ~Poisson(n)

  for (int i=1; i<num_vertices+1; i++) { // upper right of triangle (above main diagonal)
    int mean_edges=poisson_d(generator);
    poisson_distribution<> edge_d(mean_edges);
    for (int j=i+1; j<num_vertices+1; j++) {
      int num_edges=edge_d(generator);
      degrees[i]+=num_edges%2; // update degrees
      degrees[j]+=num_edges%2;

      int cur_edge=0; // ensure I D I D I D
      while (num_edges>0) { // add edges
        if (cur_edge%2==0) edge_file<<"I "<<i<<" "<<j<<endl; // insertion edge
        else edge_file<<"D "<<i<<" "<<j<<endl; // deletion edge
        num_edges-=1; cur_edge+=1;
      }
    }
    vertex_file<<i<<","<<degrees[i]<<endl;
  }

}
