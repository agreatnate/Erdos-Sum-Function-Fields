
#include <mpfr.h>   //Used for correct-rounding high-precision floats

  enum UpOrLowerBound { Upper, Lower };

class MordellSum {

public:

  MordellSum(int r, int n, int prec);

  void ComputeSum(mpfr_t retv, const int n, const int r, const int c, UpOrLowerBound upOrLow);

protected:
  // Initialize my cache arrays
  void InitializeMScache(const int r, const int N);

  // 3 dimensional cache of upper and lower bound values
  mpfr_t ***MScacheU;
  mpfr_t ***MScacheL;

  int myPrecision;

};
