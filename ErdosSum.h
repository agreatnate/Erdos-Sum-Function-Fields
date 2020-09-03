
void PrintErdSumValsWithRounding(int k, int prec, mpfr_t *erdSumsU, mpfr_t *erdSumsL);

void ComputeExactPartialSums(int q, int k, int N, mpz_class **counts, mpfr_t *erdSumsU, mpfr_t *erdSumsL);

void AddBoundsForTails(int q, int k, int N, int prec, mpz_class **counts, mpfr_t *erdSumsU, mpfr_t *erdSumsL,mpfr_t ***MScacheU,mpfr_t ***MScacheL);

void erdosKsum(int q, int k, int N, int prec);
