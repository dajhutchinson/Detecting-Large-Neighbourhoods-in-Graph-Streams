/*
 *  Read & print each edge from a file
 */
#include <iostream>
#include <fstream>

using namespace std;

int main() {

  ifstream stream("data/facebook_small.edges"); // File to read from
  string line; // current line being read
  string fst,snd;

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
    cout<<fst<<" & "<<snd<<endl; // print
  }

  return 0;
}
