/*
 *  Read & print each edge from a file
 */
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

struct edge {
  int fst;
  int snd;
};

int main() {

  ifstream stream("data/facebook_small.edges"); // File to read from
  string line; // current line being read
  string fst,snd;
  edge e;

  while (getline(stream,line)) { // Read each line
    fst=""; snd=""; // initialise temp variables
    bool after=false; // record whether after seperator (ie if we are reading first or second id)

    for (char& c:line) { // read each character
      if (c==' ') { // seperator
        after=true;
      } else if (after) { // second id
        snd+=c;
      } else { // first id
        fst+=c;
      }
    }
    // Create edge
    e.fst=stoi(fst);
    e.snd=stoi(snd);

    cout<<e.fst<<" & "<<e.snd<<endl; // print
  }

  return 0;
}
