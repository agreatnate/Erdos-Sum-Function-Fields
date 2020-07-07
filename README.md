# Erdos-Sum-Function-Fields
Compute rigorous upper and lower bounds for a sum over polynomials with a fixed number of irreducible factors over a finite field

Code used to compute the numerical values in the appendix of https://arxiv.org/abs/2007.02301 using the algorithm described in Section 6 of that paper.

Uses MPFR for arbitrary precision floating point arithmetic with correct rounding.

Requires:
 - mpfr  (https://www.mpfr.org/)
 - gmp (https://gmplib.org/)

Run make to compile.

Usage:
 erdosSum q k B p

        q = size of the finite field
        k - upper and lower bounds for the sum will be computed for each j<= k
        B - The exact value of the sum will be computed over all polynomials whose irreducible factors
                have degree at most B, and then the tails will be estimated for the contribtion from 
                degrees greater than B
        p = precision in bits to use when performing the floating point arithmetic.

Outputs upper and lower bounds for each of the computed sums, with one more digit of precision than the difference between the computed upper and lower bounds.   
