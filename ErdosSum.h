
#include <gmpxx.h>  //Used for large integer/rational classes
#include <mpfr.h>   //Used for correct-rounding high-precision floats


class ErdosSum {

public:

  ErdosSum(int q, int k, int N, int prec);

  void ComputeSum();

protected:
  void PrintErdSumValsWithRounding();

  void ComputeExactPartialSums();

  void AddBoundsForTails();


  int myQ;          // size of the set of allowed coefficients
  int myK;          // number of factors
  int myN;          // maximum degree polynomials to count exactly
  int myPrecision;  // Precision to use in sum computation
  mpfr_t *erdSumsU; // Arrays which will be used to store the computed upper/lower bounds for the Erdos Sums
  mpfr_t *erdSumsL;
  mpz_class **counts; //Array which stores the counts of the exact number of polynomials in F_q[x] having degree i+1 and j+1 irreducible factors

};
