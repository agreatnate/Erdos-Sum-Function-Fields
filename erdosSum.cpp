#include <sstream>
#include <iostream>
#include <string>
#include <gmpxx.h>  //Used for large integer/rational classes
#include <mpfr.h>   //Used for correct-rounding high-precision floats
#include "MordellSum.h"
#include "numberTheory.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;




void PrintErdSumValsWithRounding(int k, int prec, mpfr_t *erdSumsU, mpfr_t *erdSumsL){
    for(int i =1;i<=k;i++){          //Print out values of partial sums
        mpfr_t diff;
        mpfr_init2(diff,prec);
        mpfr_sub(diff,erdSumsU[i-1],erdSumsL[i-1],MPFR_RNDU);
        mpfr_log10(diff,diff,GMP_RNDU);
        long precision=(-1)*mpfr_get_si(diff,GMP_RNDU);
        string formats = "%i %."+std::to_string(precision+1)+"RDf %."+std::to_string(precision+1)+"RUf \n";
        mpfr_printf(formats.c_str(),i,erdSumsL[i-1],erdSumsU[i-1]);  //This makes sure rounding is done correctly when printing
//      cout<<i<< ' '<<erdSumsL[i-1]<<' '<<erdSumsU[i-1]<<endl;      //This is what is effectively being done here
    }
}




void ComputeExactPartialSums(int q, int k, int N, mpz_class **counts, mpfr_t *erdSumsU, mpfr_t *erdSumsL){
    //Computes the exact value of the partial Erdos Sum over F_q[x] for all polynomials with at most k factors, all with degree at most N 
    //This is the quantity S_{k,N,q} in the paper.

    for(int i=1;i<=N;i++){       //Range over all degrees up to N
       mpz_class pii(pi_q(q,i));  //Irreducibles of degree i
       mpz_class piij[k];        //List of polynomials having exactly j+1 irreducible factors of degree i (and no others)
       for(int j=1;j<=k;j++){    //First we consider the polynomials which only have irreducible factors of degree i
           piij[j-1]==0;     //Initialize
           mpz_class temp; 
           temp = pii+j-1;
           mpz_bin_ui(piij[j-1].get_mpz_t(), temp.get_mpz_t(), j);  //Number of polynomials with j irreducible factors is (pii+j-1 choose j)
           if (j<k)  counts[j-1][i*j-1]+=piij[j-1];  //Don't need to refer back to this value if j=k.  Otherwise, save it to use later.
           mpz_class qpow;
           mpz_ui_pow_ui (qpow.get_mpz_t(), q, i*j);  //qpow = q^(i*j)
           mpq_class tmpq(piij[j-1],qpow*i*j);        //tmpq = piij[j-1]/(q^(i*j)*i*j)
           mpfr_add_q(erdSumsU[j-1],erdSumsU[j-1],tmpq.get_mpq_t(),MPFR_RNDU);
           mpfr_add_q(erdSumsL[j-1],erdSumsL[j-1],tmpq.get_mpq_t(),MPFR_RNDD);
       }

       //Now consider all polynomials whose largest irreducible factor has degree i
       for(int b=k;b>1;b--){                        // Number of factors of polynomials being created
          for(int j=1;j<b;j++){                     // Number of factors of degree i being multiplied
              for(int a=(b-j);a<=(i-1)*(b-j);a++){  // Degree of the (i-1)-smooth component of the polynomial
                 if (b<k) counts[b-1][a+i*j-1]+=counts[b-j-1][a-1]*piij[j-1]; //We Won't need to refer back to this value later if b=k.
                 mpz_class qpow;                    //qpow = q^(a+i*j)
                 mpz_ui_pow_ui (qpow.get_mpz_t(), q, a+i*j);
                 mpq_class tmpq(counts[b-j-1][a-1]*piij[j-1],qpow*(a+i*j));
                 mpfr_add_q(erdSumsU[b-1],erdSumsU[b-1],tmpq.get_mpq_t(),MPFR_RNDU);
                 mpfr_add_q(erdSumsL[b-1],erdSumsL[b-1],tmpq.get_mpq_t(),MPFR_RNDD);
              }
          }
       }
    }
}



void AddBoundsForTails(int q, int k, int N, int prec, mpz_class **counts, mpfr_t *erdSumsU, mpfr_t *erdSumsL,mpfr_t ***MScacheU,mpfr_t ***MScacheL){

    mpz_class qpowB2;
    mpz_ui_pow_ui (qpowB2.get_mpz_t(), q, N/2); //qpowB2 = q^(N/2)  used frequently in the bounds below.

    mpfr_t temp;
    mpfr_init2(temp,prec);

  //This adds the term R_{k,N,q} from the paper to both the upper and lower bounds 
    //The first for loop considers the contribution when at least one of the factors has degree <= N
    for(int j=1;j<k;j++){                // Number of factors in the smooth component
        for(int a=j;a<=N*j;a++){         // Total size of the N-smooth component
            mpz_class qpow;
            mpz_ui_pow_ui (qpow.get_mpz_t(), q, a);
            for(int b=1;b<=k-j;b++){     //Number of factors of degree >N
                mpz_class factb;
                mpz_fac_ui(factb.get_mpz_t(),b);               //b factorial
                mpq_class tmpq(counts[j-1][a-1],qpow*factb);
                //First compute upper bound
                MordellSum(temp,N+1,b,a,true,prec,MScacheU,MScacheL);
                mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDU);
                mpfr_add(erdSumsU[b+j-1],erdSumsU[b+j-1],temp,MPFR_RNDU);
                //Redo Calculation for Lower Bound
                MordellSum(temp,N+1,b,a,false,prec,MScacheU,MScacheL);
                mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
                mpfr_add(erdSumsL[b+j-1],erdSumsL[b+j-1],temp,MPFR_RNDD);
            }
        }
    }
    //This for loop adds the case when all of the polynomials have degree >N
    for(int i=1;i<=k;i++){  //Add final tails
        mpz_class facti;
        mpz_fac_ui(facti.get_mpz_t(),i);               //i factorial
        MordellSum(temp,N+1,i,0,true,prec,MScacheU,MScacheL);
        mpfr_div_z(temp,temp,facti.get_mpz_t(),MPFR_RNDU);
        mpfr_add(erdSumsU[i-1],erdSumsU[i-1],temp,MPFR_RNDU);
        MordellSum(temp,N+1,i,0,false,prec,MScacheU,MScacheL);
        mpq_class tmpq(qpowB2*(q-1)-i*q,(q-1)*qpowB2*facti);
        mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
        mpfr_add(erdSumsL[i-1],erdSumsL[i-1],temp,MPFR_RNDD);
    }


    //Now subtract extra terms from the lower bound
    for(int j=1;j<k;j++){                // Number of factors in the smooth component
        for(int a=j;a<=N*j;a++){         // Total size of the N-smooth component
            mpz_class qpow;
            mpz_ui_pow_ui (qpow.get_mpz_t(), q, a);
            for(int b=1;b<=k-j;b++){     //Number of factors of degree >N
                mpz_class factb;
                mpz_fac_ui(factb.get_mpz_t(),b);               //k factorial
                MordellSum(temp,N+1,b,a,false,prec,MScacheU,MScacheL);
                mpfr_div_z(temp,temp,factb.get_mpz_t(),MPFR_RNDD);
                mpq_class tmpq(counts[j-1][a-1]*b*q,qpow*(q-1)*qpowB2);
                mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
                mpfr_sub(erdSumsL[b+j-1],erdSumsL[b+j-1],temp,MPFR_RNDD);
            }
        }
    }

    //Finally add extra terms to upper bound (To deal with non-square-free polynomials)
    for(int i=1;i<=k;i++){  //Add final tails to upper bounds
        mpz_class qpow;
        mpz_ui_pow_ui (qpow.get_mpz_t(), q, N);
        mpq_class tmpq(2,N*qpow);
        mpfr_add_q(erdSumsU[i-1],erdSumsU[i-1],tmpq.get_mpq_t(),MPFR_RNDU);
    }
}




void erdosKsum(int q, int k, int N, int prec) {
    //Computes upper and lower bounds for the sum of 1/(deg(f) * q^deg(f) ) ranging over all polynomials f in F_q
    //with exactly j irreducible factors for each j<=k.  The bounds are computed by first computing the exact contribution
    //coming from polynomials, all of whose irreducible factors are at most N, and then bounding the size of the tails


    mpfr_t ***MScacheU;  //These arrays store cached upper and lower bounds for the Mordell Sums
    mpfr_t ***MScacheL;

    initializeMScache(k,N,prec,MScacheU,MScacheL);


    mpfr_t erdSumsU[k];  //Arrays storing computed upper and lower bounds for the sums with j irreducible factors
    mpfr_t erdSumsL[k];
    mpz_class **counts = new mpz_class*[k-1];  //Array to store counts of polynomials with given factorizations (See below)

    for (int i =0;i<k;i++){   //Initialize erdSum Arrays to 0
      mpfr_init2(erdSumsU[i],prec);
      mpfr_set_ui(erdSumsU[i],0,MPFR_RNDU);
      mpfr_init2(erdSumsL[i],prec);
      mpfr_set_ui(erdSumsL[i],0,MPFR_RNDD);
    }

    for (int i =0;i<k-1;i++){
      counts[i] = new mpz_class[N*(i+1)];
      for(int j=0;j<N*(i+1);j++){
        counts[i][j]=0;          //This will store the number of irreducible polynomials with [i+1] irreducible factors of degree [j+1]
      }
    }

    ComputeExactPartialSums(q,k,N,counts,erdSumsU,erdSumsL);

    cout<<"Computed lower and upper bounds for the partial sum:"<<endl;
    PrintErdSumValsWithRounding(k, prec, erdSumsU, erdSumsL);

    //Finished creating the exact (partial) sum, now taking into account the bounds for the tails

    AddBoundsForTails(q,k,N,prec,counts,erdSumsU,erdSumsL,MScacheU,MScacheL);

    cout<<"Computed upper and lower bounds for the entire sum:"<<endl;
    PrintErdSumValsWithRounding(k, prec, erdSumsU, erdSumsL);
}




