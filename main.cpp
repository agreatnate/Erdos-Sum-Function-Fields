#include <iostream>
#include <gmpxx.h>  //Used for large integer/rational classes
#include <mpfr.h>   //Used for correct-rounding high-precision floats
#include "erdosSum.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;





int main (int argc, char *argv[]){
    if(argc != 5){
      cout<<"This program will compute upper and lower bounds for the value of the sum of 1/(deg(f)*q^(deg(f))"<<endl;
      cout<<"\t ranging over all polynomials f with coeffecients in the finite field F_q having a fixed number (j)"<<endl<<"\t of irreducible factors."<<endl;
      cout<<"Syntax:"<<endl;
      cout<<'\t'<<"erdosSum q k N p"<<endl<<endl;
      cout<<'\t'<<"q = size of the finite field"<<endl;
      cout<<'\t'<<"k - upper and lower bounds for the sum will be computed for each j<= k"<<endl;
      cout<<'\t'<<"N - The exact value of the sum will be computed over all polynomials whose irreducible factors"<<endl<<"\t\thave degree at most N, and then the tails will be estimated for the contribtion from "<<endl<<"\t\tdegrees greater than N (N must be even)"<<endl;
      cout<<'\t'<<"p = precision in bits to use when performing the floating point arithmetic."<<endl;
      return 0;
    }
    int q=atoi(argv[1]); //size of Finite Field to consider
    int k=atoi(argv[2]); //maximum degree to consider
    int N=atoi(argv[3]); //Bound size of irreducibles to precompute
    int p=atoi(argv[4]); //precision in bytes
    if (N%2==1){
       cerr<<"Only implemented for even values of N."<<endl;
       return 1;
    }

//    cout<<N<<'\t'<<k<<'\t'<<q<<'\t'<<MordellSum(N,k,q)<<endl;
    erdosKsum(q,k,N,p);

    return 0;
}
