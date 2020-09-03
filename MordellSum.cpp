#include "MordellSum.h"
#include <gmpxx.h>  //Used for large integer/rational classes
//#include <mpfr.h>   //Used for correct-rounding high-precision floats

#define IS_ODD(n) (n % 2)

//MordellSum Constructor
MordellSum::MordellSum(int r, int N, int prec) : myPrecision(prec){
  InitializeMScache(r, N);
}


void MordellSum::InitializeMScache(const int r, const int N){
                                                 //Initialize tables of entries of Cached values of Mordell Sums
                                                 //These are stored in MSCachU and MSCacheL respectively
                                                 //Each entry is initialized with myPrecision of bytes and an initial value of 0
    MScacheU = new mpfr_t**[r];
    MScacheL = new mpfr_t**[r];
    for (int i =0;i<r-1;i++){
      MScacheU[i] = new mpfr_t*[N+1];
      MScacheL[i] = new mpfr_t*[N+1];
      for(int j=0;j<N+1;j++){
        int l=(r-(i+1))*(N+1);
        MScacheU[i][j]= new mpfr_t[l+1];
        MScacheL[i][j]= new mpfr_t[l+1];
        for(int k=0;k<l+1;k++){
          mpfr_init2(MScacheU[i][j][k],myPrecision);
          mpfr_set_ui(MScacheU[i][j][k],0,MPFR_RNDU);
          mpfr_init2(MScacheL[i][j][k],myPrecision);
          mpfr_set_ui(MScacheL[i][j][k],0,MPFR_RNDD);
        }
      }
    }
    MScacheU[r-1] = new mpfr_t*[N+1];
    MScacheL[r-1] = new mpfr_t*[N+1];
    int i=r-1;
    for(int j=0;j<N+1;j++){
      MScacheU[i][j]= new mpfr_t[1];
      MScacheL[i][j]= new mpfr_t[1];
      mpfr_init2(MScacheU[i][j][0],myPrecision);
      mpfr_set_ui(MScacheU[i][j][0],0,MPFR_RNDU);
      mpfr_init2(MScacheL[i][j][0],myPrecision);
      mpfr_set_ui(MScacheL[i][j][0],0,MPFR_RNDD);
    }
}









void MordellSum::ComputeSum(mpfr_t retv, const int n, const int k, const int c, UpOrLowerBound upOrLower){

    //Computes a bound (Upper/Lower depending on upOrLower) for
    //the value of sum_{n<=i_1,i_2,ldots i_k} 1/(i_1i_2 \cdots i_k(i_1+i_2+\cdots i_k +c))
    //This is the function M(k,n,c) appearing in the paper

    //Returned value is stored in retv

    mpfr_rnd_t  RoundDir;
    if(upOrLower == Upper) {
       RoundDir = MPFR_RNDU;
    }
    else {
       RoundDir = MPFR_RNDD;
    }
    mpfr_init2(retv,myPrecision);

    mpz_class factk;                           //k factorial
    mpz_fac_ui(factk.get_mpz_t(),k);


    if (k==0){                     //Trivial case
        mpq_class frac(1,c);
        mpfr_set_q(retv,frac.get_mpq_t(),RoundDir);
        return;
    }

    if (mpfr_cmp_ui (MScacheU[k-1][n-1][c],0)>0) {             //If this is a nonzero value then it has already been computed
                                                               //In that case we just return the cached value
        if(upOrLower == Upper){
           mpfr_set(retv,MScacheU[k-1][n-1][c],MPFR_RNDU);
           return;
        }
        else{
           mpfr_set(retv,MScacheL[k-1][n-1][c],MPFR_RNDD);
           return;
        }
    }

    // Otherwise, the value needs to be computed
    if (c==0 && n==1){                                //The special case when c=0 and n=1 results in an
                                                      //infinite sum which must be treated separately
         mpfr_t z;
         mpfr_init2(z,myPrecision);

         //First do the calculation, rounding up
         mpfr_zeta_ui(z,k+1,MPFR_RNDU);
         mpfr_mul_z (z, z, factk.get_mpz_t(), MPFR_RNDU); //z=k!*zeta(k+1) (rounded up)
         mpfr_set(MScacheU[k-1][n-1][0],z,MPFR_RNDU);

         //Redo the calculation, rounding down
         mpfr_zeta_ui(z,k+1, MPFR_RNDD);
         mpfr_mul_z (z, z, factk.get_mpz_t(), MPFR_RNDD); //z=k!*zeta(k+1) (rounded down)
         mpfr_set(MScacheL[k-1][n-1][0],z,MPFR_RNDD);
    }


    //if c>0 and n=1, Mordell's theorem results in a finite sum, and hence a rational number which we compute exactly.
    else if (n==1) {
        mpq_class s(0);                                           //running sum will be stored here

        //s = sum_{0<=i<=c} (-1)^i * (c-1 choose i) * 1/(i+1)^(k+1)
        for(int i=0;i<c;i++){
           mpz_class tmpz;
           tmpz = (i+1);
           mpz_pow_ui (tmpz.get_mpz_t(), tmpz.get_mpz_t(), k+1);  //tmpz= (i+1)^(k+1)
           mpq_class tmpq(1,tmpz);                                //tmpq = 1/(i+1)^(k+1)
           mpz_class binomz;
           int sign=1;
           if(IS_ODD(i)) sign=-1;
           mpz_bin_uiui (binomz.get_mpz_t(), c-1, i);             //binomz = (c-1 choose i)
           s+=sign*binomz*tmpq;
        }
        s*=factk;  //Multiply by k!
        mpfr_set_q(MScacheU[k-1][n-1][c],s.get_mpq_t(),GMP_RNDU);  //Save off computed values for use later
        mpfr_set_q(MScacheL[k-1][n-1][c],s.get_mpq_t(),GMP_RNDD);
    }    else{   //otherwise, if n > 1 we compute the sum by recursively using previous values of the sum for smaller values of n

        // M(k,n,c) = \sum_{i=1}^k (k choose i) (-1)^i M(k-i,n-1,a+i(n-1)) / (n-1)^i
        mpfr_t su;  //upper bound of sum
        mpfr_init2(su,myPrecision);
        mpfr_set_ui(su,0,MPFR_RNDU);
        mpfr_t sl;  //lower bound of sum
        mpfr_init2(sl,myPrecision);
        mpfr_set_ui(sl,0,MPFR_RNDD);
        for(int i=0;i<=k;i++){
           mpz_class tmpz(n-1);
           mpz_pow_ui (tmpz.get_mpz_t(), tmpz.get_mpz_t(), i);
           mpq_class tmpq(1,tmpz);
           mpz_class binomz;
           mpz_bin_uiui (binomz.get_mpz_t(), k, i);
           tmpq*=binomz;
           mpfr_t temp;  //Used to store return value from MordellSum
           mpfr_init2(temp,myPrecision);
           if(IS_ODD(i)){
              //sign=-1;
              ComputeSum(temp,n-1,k-i,c+(n-1)*i,Lower);       //First Round Down
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
              mpfr_sub(su,su,temp,MPFR_RNDU);                 //+=(-1)*binomz*tmpq*MordellSum(n-1,k-i,c+(n-1)*i,false,myPrecision);
              ComputeSum(temp,n-1,k-i,c+(n-1)*i,Upper);       //Now Round Up
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDU);
              mpfr_sub(sl,sl,temp,MPFR_RNDD);                 //+=(-1)*binomz*tmpq*MordellSum(n-1,k-i,c+(n-1)*i,false,myPrecision);
           }
           else{
              //sign=1;
              ComputeSum(temp,n-1,k-i,c+(n-1)*i,Upper);      //First Round Up
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDU);
              mpfr_add(su,su,temp,MPFR_RNDU);                //+=(-1)*binomz*tmpq*MordellSum(n-1,k-i,c+(n-1)*i,false,myPrecision);
              ComputeSum(temp,n-1,k-i,c+(n-1)*i,Lower);      //Now Round Down
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
              mpfr_add(sl,sl,temp,GMP_RNDD);                 //+=(-1)*binomz*tmpq*MordellSum(n-1,k-i,c+(n-1)*i,false,myPrecision);
           }
        }
        mpfr_set(MScacheU[k-1][n-1][c],su,MPFR_RNDU);  //Save off computed values for use later
        mpfr_set(MScacheL[k-1][n-1][c],sl,MPFR_RNDD);
    }

    //Return requested value
    if(upOrLower = Upper){  //Asked to return an upper bound
       mpfr_set(retv,MScacheU[k-1][n-1][c],MPFR_RNDU);
       return;
    }
    else{    //Asked to return a lower bound
       mpfr_set(retv,MScacheL[k-1][n-1][c],MPFR_RNDD);
       return;
    }
}




