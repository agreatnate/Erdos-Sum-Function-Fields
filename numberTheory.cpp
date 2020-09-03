#include <gmpxx.h>  //Used for large integer/rational classes
#include <mpfr.h>   //Used for correct-rounding high-precision floats
#include "math.h"

#define IS_EVEN(n) !(n % 2)

// Returns value of mobius(n)
int mobius(int n){
    if (n==1) return 1;
    int p = 0;   //Counts the number of distinct prime factors of n
    // Treat 2 separately
    if (IS_EVEN(n)) {
        n = n/2;
        p++;
        if (IS_EVEN(n))
           return 0;
    }
    // Check for all other prime factors
    for (int i = 3; i <= n; i = i+2)     {
        // If i divides n
        if (n % i == 0)  {
            n = n/i;
            p++;
            // If i^2 also divides N
            if (n % i == 0)
               return 0;
        }
    }
    if(IS_EVEN(p)) return 1;
    else return -1;
}


mpz_class pi_q(const int q, int i) {
   //Returns the exact number of irreducible polynomials of degree i over F_q
   mpz_class res(0);
   double sqrtI= sqrt(i);
   for (int j=1;j<=sqrtI;j++){
      if (i % j==0){
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



