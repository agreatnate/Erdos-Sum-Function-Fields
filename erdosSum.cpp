#include <sstream>
#include <iostream>
#include <algorithm>
#include <complex>
#include <string>
#include <gmpxx.h>
#include <mpfr.h>

using namespace std;

mpfr_t ***MScacheU;  //These arrays store cached upper and lower bounds for the Mordell Sums
mpfr_t ***MScacheL;

void MordellSum(mpfr_t retv, int n, int r,int c,bool up,int prec);

// Returns value of mobius(n)
int mobius(int n){
    if (n==1) return 1;
    int p = 0;
    // Handling 2 separately 
    if (n%2 == 0) {
        n = n/2; 
        p++; 
        if (n % 2 == 0) 
           return 0; 
    } 
    // Check for all other prime factors 
    for (int i = 3; i <= n; i = i+2)     {
        // If i divides n 
        if (n%i == 0)  { 
            n = n/i; 
            p++; 
            // If i^2 also divides N 
            if (n % i == 0) 
               return 0; 
        }
    }
    return (p % 2 == 0)? 1 : -1; 
}

mpz_class piq(int q, int i) {  
   //Returns the exact number of irreducible polynomials of degree i over F_q
   mpz_class res(0);
   for (int j=1;j<=sqrt(i);j++){
      if (i%j==0){
         mpz_class r;
         mpz_ui_pow_ui (r.get_mpz_t(), q, j);
         res+=r*mobius(i/j);
         if (i/j!=j){
             mpz_ui_pow_ui (r.get_mpz_t(), q, i/j);
             res+=r*mobius(j);
         }
      }
   }
   return res/i;
}



void LowBdRem(mpfr_t retv, int B,int k,int C,int prec){    //Lower bound for sum_{B<=i_1<=i_2<ldots l_k}1/(i_1i_2 \cdots l_k(l_1+l_2+\cdots l_k +C))
                                                    //Note that the integers are NOT necessarily distinct in this sum
     //return value is stored in retv
     mpfr_init2(retv,prec);
     mpz_class factk;
     mpz_fac_ui(factk.get_mpz_t(),k);               //k factorial
     mpfr_t v;
     mpfr_init2(v,prec);
     MordellSum(v,B,k,C,false,prec);
     mpfr_div_z(retv,v,factk.get_mpz_t(),MPFR_RNDD);
}


void UpBdRem(mpfr_t retv,int B,int k,int C,int prec){ //Upper bound for sum_{B<=i_1 < i_2 < ldots l_k}1/(i_1i_2 \cdots l_k(l_1+l_2+\cdots l_k +C))
                                                //Note that the integers are distinct in this sum
     //return value is stored in retv
     mpfr_init2(retv,prec);
     mpz_class factk;
     mpz_fac_ui(factk.get_mpz_t(),k);               //k factorial
     mpfr_t v;
     mpfr_init2(v,prec);
     MordellSum(v,B,k,C,true,prec);
     mpfr_div_z(retv,v,factk.get_mpz_t(),MPFR_RNDU);
}



void erdKsum(int q, int k, int B, int prec) {
     //Computes upper and lower bounds for the sum of 1/(deg(f) * q^deg(f) ) ranging over all polynomials f in F_q
     //with exactly j irreducible factors for each j<=k.  The bounds are computed by first computing the exact contribution
     //coming from polynomials, all of whose irreducible factors are at most B, and then estimating the size of the tails

    mpfr_t erdSumsU[k];  //Arrays storing computed upper and lower bounds for the sums with j irreducible factors
    mpfr_t erdSumsL[k];

    for (int i =0;i<k;i++){
      mpfr_init2(erdSumsU[i],prec);
      mpfr_set_ui(erdSumsU[i],0,MPFR_RNDU);
      mpfr_init2(erdSumsL[i],prec);
      mpfr_set_ui(erdSumsL[i],0,MPFR_RNDD);
    }

    mpz_class **counts = new mpz_class*[k-1];
    for (int i =0;i<k-1;i++){
      counts[i] = new mpz_class[B*(i+1)];
      for(int j=0;j<B*(i+1);j++){
        counts[i][j]=0;          //Number of irreducible polynomials with [i+1] irreducible factors of degree [j+1]
      }
    }


    for(int i=1;i<=B;i++){
       mpz_class pii(piq(q,i));  //Irreducibles of degree i
       mpz_class piij[k];        //List of polynomials having exactly j+1 irreducible factors of degree i (and no others)
       for(int j=1;j<=k;j++){    //Consider the polynomials which only have irreducible factors of degree i
         piij[j-1]==0;
         mpz_class temp;
         temp =(pii+j-1);
         mpz_bin_ui(piij[j-1].get_mpz_t(), temp.get_mpz_t(), j);
         if (j<k)  counts[j-1][i*j-1]+=piij[j-1];  //Don't need to refer back to this value if j=k.
         mpz_class qpow;
         mpz_ui_pow_ui (qpow.get_mpz_t(), q, i*j);
         mpq_class tmpq(piij[j-1],qpow*i*j);
         mpfr_add_q(erdSumsU[j-1],erdSumsU[j-1],tmpq.get_mpq_t(),MPFR_RNDU);
         mpfr_add_q(erdSumsL[j-1],erdSumsL[j-1],tmpq.get_mpq_t(),MPFR_RNDD);
       }

       //Now consider polynomials whose largest irreducible factor has degree i
       for(int b=k;b>1;b--){// in range(k,1,-1): #Number of factors of things being created
          for(int j=1;j<b;j++){// in [1..b-1]: #number of copies of things of size i being multiplied
              for(int a=(b-j);a<=(i-1)*(b-j);a++){//]: #degree of the i-smooth component
                 if (b<k) counts[b-1][a+i*j-1]+=counts[b-j-1][a-1]*piij[j-1]; //Don't need to refer back to this value if b=k.
                 mpz_class qpow;
                 mpz_ui_pow_ui (qpow.get_mpz_t(), q, a+i*j);
                 mpq_class tmpq(counts[b-j-1][a-1]*piij[j-1],qpow*(a+i*j));
                 mpfr_add_q(erdSumsU[b-1],erdSumsU[b-1],tmpq.get_mpq_t(),MPFR_RNDU);
                 mpfr_add_q(erdSumsL[b-1],erdSumsL[b-1],tmpq.get_mpq_t(),MPFR_RNDD);
              }
          }
       }
    }
    cout<<"Computed lower and upper bounds for the partial sum:"<<endl;
    for(int i =1;i<=k;i++){
        mpfr_t diff;
        mpfr_init2(diff,prec);
        mpfr_sub(diff,erdSumsU[i-1],erdSumsL[i-1],MPFR_RNDD);
        mpfr_log10(diff,diff,GMP_RNDU);
        long precision=(-1)*mpfr_get_si(diff,GMP_RNDU);
        string formats = "%i %."+std::to_string(precision+1)+"RDf %."+std::to_string(precision+1)+"RUf \n";
        mpfr_printf(formats.c_str(),i,erdSumsL[i-1],erdSumsU[i-1]);
//      cout<<i<< ' '<<erdSumsL[i-1]<<' '<<erdSumsU[i-1]<<endl;
    }


    //Finished creating the exact (partial) sum, now taking into account the tails

    mpz_class qpowB2;
    mpz_ui_pow_ui (qpowB2.get_mpz_t(), q, B/2); //qpowB2 = q^(B/2)  used frequently in the bounds below.

    for(int j=1;j<k;j++){// in [1..k]:  //Number of factors in the smooth component
        for(int a=j;a<=B*j;a++){// in [1..B*j]: //Total size of the B-smooth component
            mpz_class qpow;
            mpz_ui_pow_ui (qpow.get_mpz_t(), q, a);
            for(int b=1;b<=k-j;b++){// in [1..k-j]: //Number of factors of degree >B
                mpfr_t temp;
                mpfr_init2(temp,prec);
                UpBdRem(temp,B+1,b,a,prec);
                mpq_class tmpq(counts[j-1][a-1],qpow);
                mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDU);
                mpfr_add(erdSumsU[b+j-1],erdSumsU[b+j-1],temp,MPFR_RNDU);
                LowBdRem(temp,B+1,b,a,prec);
                tmpq = mpq_class(counts[j-1][a-1]*(qpowB2*(q-1)-b*q),qpow*(q-1)*qpowB2);
                mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
                mpfr_add(erdSumsL[b+j-1],erdSumsL[b+j-1],temp,MPFR_RNDD);
            }
        }
    }
    cout<<"Computed upper and lower bounds for the entire sum:"<<endl;

    for(int i=1;i<=k;i++){
        mpz_class qpow;
        mpz_ui_pow_ui (qpow.get_mpz_t(), q, B);
        mpfr_t temp;
        mpfr_init2(temp,prec);
        UpBdRem(temp,B+1,i,0,prec);
        mpfr_add(erdSumsU[i-1],erdSumsU[i-1],temp,MPFR_RNDU);
        mpq_class tmpq(2,B*qpow*(q-1));
        mpfr_add_q(erdSumsU[i-1],erdSumsU[i-1],tmpq.get_mpq_t(),MPFR_RNDU);
        LowBdRem(temp,B+1,i,0,prec);
        tmpq = mpq_class(qpowB2*(q-1)-i*q,(q-1)*qpowB2);
        mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
        mpfr_add(erdSumsL[i-1],erdSumsL[i-1],temp,MPFR_RNDD);
        mpfr_t diff;
        mpfr_init2(diff,prec);
        mpfr_sub(diff,erdSumsU[i-1],erdSumsL[i-1],MPFR_RNDD);
        mpfr_log10(diff,diff,GMP_RNDU);
        long precision=(-1)*mpfr_get_si(diff,GMP_RNDU);
        string formats = "%i %."+std::to_string(precision+1)+"RDf %."+std::to_string(precision+1)+"RUf \n";
        mpfr_printf(formats.c_str(),i,erdSumsL[i-1],erdSumsU[i-1]);
    }
}

void initMScache(int r, int B, int p){  //Initialize tables of entries of Cached values of Mordell Sums
    MScacheU = new mpfr_t**[r];
    MScacheL = new mpfr_t**[r];
    for (int i =0;i<r-1;i++){
      MScacheU[i] = new mpfr_t*[B+1];
      MScacheL[i] = new mpfr_t*[B+1];
      for(int j=0;j<B+1;j++){
        int l=(r-(i+1))*(B+1);
        MScacheU[i][j]= new mpfr_t[l+1];
        MScacheL[i][j]= new mpfr_t[l+1];
        for(int k=0;k<l+1;k++){
          mpfr_init2(MScacheU[i][j][k],p);
          mpfr_set_ui(MScacheU[i][j][k],0,MPFR_RNDU);
          mpfr_init2(MScacheL[i][j][k],p);
          mpfr_set_ui(MScacheL[i][j][k],0,MPFR_RNDD);
        }
      }
    }
    MScacheU[r-1] = new mpfr_t*[B+1];
    MScacheL[r-1] = new mpfr_t*[B+1];
    int i=r-1;
    for(int j=0;j<B+1;j++){
      MScacheU[i][j]= new mpfr_t[1];
      MScacheL[i][j]= new mpfr_t[1];
      mpfr_init2(MScacheU[i][j][0],p);
      mpfr_set_ui(MScacheU[i][j][0],0,MPFR_RNDU);
      mpfr_init2(MScacheL[i][j][0],p);
      mpfr_set_ui(MScacheL[i][j][0],0,MPFR_RNDD);
    }
}

void MordellSum(mpfr_t retv, int n,int r,int c,bool up,int prec){
    //Computes a bound (Upper/Lower) for sum_{n<=i_1,i_2,ldots l_r, distinct??} 1/(i_1i_2 \cdots l_r(l_1+l_2+\cdots l_r +c))
    mpfr_rnd_t  RoundDir;
    if(up) RoundDir = MPFR_RNDU;
    else RoundDir = MPFR_RNDD;
    mpfr_init2(retv,prec);
    if (r==0){
        mpq_class frac(1,c);
        mpfr_set_q(retv,frac.get_mpq_t(),RoundDir);
        return;
    }
    if (c==0){                                        //Special case when c=0, get infinite sum 
       if (n==1){
           mpfr_t z;
           mpfr_init2(z,prec);
           mpz_class factr;                           //r factorial
           mpz_fac_ui(factr.get_mpz_t(),r);
           mpfr_zeta_ui(z,r+1,MPFR_RNDU);
           mpfr_mul_z (z, z, factr.get_mpz_t(), MPFR_RNDU); //z=r!*zeta(r+1) (rounded up)
           mpfr_set(MScacheU[r-1][n-1][0],z,MPFR_RNDU);
           //Redo the calculation, rounding down
           mpfr_zeta_ui(z,r+1, MPFR_RNDD);
           mpfr_mul_z (z, z, factr.get_mpz_t(), MPFR_RNDD); //z=r!*zeta(r+1) (rounded down)
           mpfr_set(MScacheL[r-1][n-1][0],z,MPFR_RNDD);
           if(up){
              mpfr_set(retv,MScacheU[r-1][n-1][0],MPFR_RNDU);
              return;
           }
           else{
              mpfr_set(retv,MScacheL[r-1][n-1][0],MPFR_RNDD);
              return;
           }
       }
    }
    if (mpfr_cmp_ui (MScacheU[r-1][n-1][c],0)>0) {                      //Already Computed
        if(up){
           mpfr_set(retv,MScacheU[r-1][n-1][c],MPFR_RNDU);
           return;
        }
        else{
           mpfr_set(retv,MScacheL[r-1][n-1][c],MPFR_RNDD);
           return;
        }
    }

    if (n==1) {
        mpq_class s(0);
        for(int i=0;i<c;i++){
           mpz_class tmpz;
           tmpz = (i+1);
           mpz_pow_ui (tmpz.get_mpz_t(), tmpz.get_mpz_t(), r+1);  //tmpz= (i+1)^(r+1)
           mpq_class tmpq(1,tmpz);                                //tmpq = 1/(i+1)^(r+1)
           mpz_class binomz;
           int sign=1;
           if(i%2==1) sign=-1;
           mpz_bin_uiui (binomz.get_mpz_t(), c-1, i);             //binomz = (c-1 choose i)
           s+=sign*binomz*tmpq;
        }
                                                                  //s = sum_{0<=i<=c} (-1)^i * (c-1 choose i) * 1/(i+1)^(r+1)
        //v=factorial(r)*sum([(-1)^i/(i+1)^(r+1)*binomial(c-1,i) for i in [0..c-1]])
        mpz_class factr;
        mpz_fac_ui(factr.get_mpz_t(),r);
        s*=factr;
//        cout<<'!'<<n<<'\t'<<r<<'\t'<<c<<endl;
        mpfr_set_q(MScacheU[r-1][n-1][c],s.get_mpq_t(),GMP_RNDU);
        mpfr_set_q(MScacheL[r-1][n-1][c],s.get_mpq_t(),GMP_RNDD);
        if(up){
           mpfr_set(retv,MScacheU[r-1][n-1][c],MPFR_RNDU);
           return;
        }
        else{
           mpfr_set(retv,MScacheL[r-1][n-1][c],MPFR_RNDD);
           return;
        }
    }else{
        mpfr_t su;  //upper bound of sum
        mpfr_init2(su,prec);
        mpfr_set_ui(su,0,MPFR_RNDU);
        mpfr_t sl;  //lower bound of sum
        mpfr_init2(sl,prec);
        mpfr_set_ui(sl,0,MPFR_RNDD);
        for(int j=0;j<=r;j++){
           mpz_class tmpz(n-1);
           mpz_pow_ui (tmpz.get_mpz_t(), tmpz.get_mpz_t(), j);
           mpq_class tmpq(1,tmpz);
           mpz_class binomz;
//           int sign=1;
           mpz_bin_uiui (binomz.get_mpz_t(), r, j);
           tmpq*=binomz;
           mpfr_t temp;
           mpfr_init2(temp,prec);
           if(j%2==1){
              //sign=-1;
              MordellSum(temp,n-1,r-j,c+(n-1)*j,false,prec);
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
              mpfr_sub(su,su,temp,MPFR_RNDU);//+=(-1)*binomz*tmpq*MordellSum(n-1,r-j,c+(n-1)*j,false,prec);
              MordellSum(temp,n-1,r-j,c+(n-1)*j,true,prec);
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDU);
              mpfr_sub(sl,sl,temp,MPFR_RNDD);//+=(-1)*binomz*tmpq*MordellSum(n-1,r-j,c+(n-1)*j,false,prec);
           }
           else{
              //sign=1;
              MordellSum(temp,n-1,r-j,c+(n-1)*j,true,prec);
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDU);
              mpfr_add(su,su,temp,MPFR_RNDU);//+=(-1)*binomz*tmpq*MordellSum(n-1,r-j,c+(n-1)*j,false,prec);
              MordellSum(temp,n-1,r-j,c+(n-1)*j,false,prec);
              mpfr_mul_q(temp,temp,tmpq.get_mpq_t(),MPFR_RNDD);
              mpfr_add(sl,sl,temp,GMP_RNDD);//+=(-1)*binomz*tmpq*MordellSum(n-1,r-j,c+(n-1)*j,false,prec);
           }
        }
        mpfr_set(MScacheU[r-1][n-1][c],su,MPFR_RNDU);
        mpfr_set(MScacheL[r-1][n-1][c],su,MPFR_RNDD);
        if(up){
           mpfr_set(retv,su,MPFR_RNDU);
           return;
        }
        else{
           mpfr_set(retv,sl,MPFR_RNDD);
           return;
        }
    }
}


int main (int argc, char *argv[]){
    if(argc != 5){
      cout<<"This program will compute upper and lower bounds for the value of the sum of 1/(deg(f)*q^(deg(f))"<<endl;
      cout<<"\t ranging over all polynomials f with coeffecients in the finite field F_q having a fixed number (j)"<<endl<<"\t of irreducible factors."<<endl;
      cout<<"Syntax:"<<endl;
      cout<<'\t'<<"erdosSum q k B p"<<endl<<endl;
      cout<<'\t'<<"q = size of the finite field"<<endl;
      cout<<'\t'<<"k - upper and lower bounds for the sum will be computed for each j<= k"<<endl;
      cout<<'\t'<<"B - The exact value of the sum will be computed over all polynomials whose irreducible factors"<<endl<<"\t\thave degree at most B, and then the tails will be estimated for the contribtion from "<<endl<<"\t\tdegrees greater than B"<<endl;
      cout<<'\t'<<"p = precision in bits to use when performing the floating point arithmetic."<<endl;
      return 0;
    }
    cout.precision(50);
    cout<<fixed;
    int q=atoi(argv[1]); //size of Finite Field to consider
    int k=atoi(argv[2]); //maximum degree to consider
    int B=atoi(argv[3]); //Bound size of irreducibles to precompute
    int p=atoi(argv[4]); //precision in bytes
    initMScache(k,B,p);
//    cout<<B<<'\t'<<k<<'\t'<<q<<'\t'<<MordellSum(B,k,q)<<endl;
    erdKsum(q,k,B,p);

    return 0;
}
