/*
 *  An l0-sample of a vector, x, is a uniform sample from the non-zero elements of x
 *  NOTE
 *  hash function should be chosen from "k-wise independent family of hash functions,"???
 *  Always recovers same index (WHY???) bad hash function
 */

#include <ctime>
#include <iostream>
#include <random>

using namespace std;
int P=1073741789; // >2^30

 /*-----------------*
  * DATA STRUCTURES *
  *-----------------*/

// parameters for hash function
 struct hash_params {
   unsigned long a;
   unsigned long b;
   unsigned long m;
 };

 /*-----------*
 * SIGNATURES *
 *------------*/

// hashing
hash_params generate_hash(int m);
int hash_function(int key, hash_params ps);

// l0
int* sample(int* a, int size, int j, hash_params ps);
int recover(int* aj, int size);

/*-----*
* BODY *
*------*/

int main() {
  int a[]={0,0,0,1,0,1,1,0,1,0,0,0,1,0,0,1,1,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,0,1,1,0,0,1,1,1,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,1,1,
           0,1,1,0,1,1,1,1,0,1,0,0,0,1,1,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,1,
           0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,0,1,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,
           1,1,0,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,0,0,0,0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0,1,1,1,0,
           1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0,1,0,1,1,1,1,0,1,0,1,1,1,0,0,1,1,0,0,1,0,1,1,0,
           1,1,1,0,0,1,0,0,1,0,1,1,0,1,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,0,1,0,0,1,1,0,1,0,
           0,0,0,1,0,1,1,0,1,0,0,0,1,0,0,1,1,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,0,1,1,0,0,1,1,1,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,1,1,
           0,1,1,0,1,1,1,1,0,1,0,0,0,1,1,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,1,
           0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,0,1,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,
           1,1,0,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,0,0,0,0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0,1,1,1,0,
           1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0,1,0,1,1,1,1,0,1,0,1,1,1,0,0,1,1,0,0,1,0,1,1,0,
           1,1,1,0,0,1,0,0,1,0,1,1,0,1,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,0,1,0,0,1,1,0,1,0}; // 402 1s, 367 0s
  int size=sizeof(a)/sizeof(int);

  hash_params ps=generate_hash(pow(size,3)); // m=n^3 where n is |x|
  cout<<"{a="<<ps.a<<",b="<<ps.b<<",m="<<ps.m<<"}"<<endl;

  int j=1;
  while (true) {
    //cout<<"j="<<j<<endl;

    int* a_j=sample(a,size,j,ps); // sample from vector
    //for (int i=0; i<20; i++) cout << a_j[i]<<",";
    //cout<<endl;

    int r=recover(a_j,size); // try to recover an index
    if (r==-2) { // fail
      cout<<"FAIL-a_"<<j<<" is zero string"<<endl;
      break;
    } else if (r!=-1) { // success
      cout<<"j="<<j<<",r="<<r<<",a[r]="<<a[r]<<endl;
      break;
    } else j+=1;// try again
  }

  return 0;
}

// generate parameters to use in hash function
hash_params generate_hash(int m) {
  default_random_engine generator{static_cast<long unsigned int>(time(0))};
  uniform_int_distribution<unsigned long> distribution(0,P-1);
  hash_params ps={
    distribution(generator),
    distribution(generator),
    (unsigned long)m
  };
  return ps;
}

// hash a key
int hash_function(int key, hash_params ps) {
  return ((ps.a*key+ps.b)%P)%ps.m;
}

// sample from vector a with condition j
int* sample(int* a, int size, int j, hash_params ps) {
  int* a_j= new int[size];
  for (int i=0; i<size; i+=1) {
    int h=hash_function(i,ps);
    if (h<=ps.m/pow(2,j)) a_j[i]=a[i];
    else a_j[i]=0;
  }
  return a_j;
}

// recover a non-zero index of a if a is 1-sparse
int recover(int* aj, int size) {
  int w1=0, w2=0; // w1=sum aj_i, w_2=sum i*aj_i

  for (int i=0; i<size; i++) {
      w1+=aj[i];
      w2+=i*aj[i];
  }
  cout<<"w1="<<w1<<",w2="<<w2<<endl;
  if (w1==1) return w2; // 1-sparse
  if (w1==0) return -2; // zero string
  return -1; // not 1-sparse
}
