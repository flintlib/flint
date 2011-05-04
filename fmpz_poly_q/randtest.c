#include "fmpz_poly_q.h"

void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state,
                          long len1, mp_bitcnt_t bits1, 
                          long len2, mp_bitcnt_t bits2)
{
    len2  = FLINT_MAX(len2, 1);
    bits2 = FLINT_MAX(bits2, 1);

    fmpz_poly_randtest(poly->num, state, len1, bits1);
    fmpz_poly_randtest_not_zero(poly->den, state, len2, bits2);
    fmpz_poly_q_canonicalise(poly);
}

void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, 
                                   long len1, mp_bitcnt_t bits1, 
                                   long len2, mp_bitcnt_t bits2)
{
    len1  = FLINT_MAX(len1, 1);
    len2  = FLINT_MAX(len2, 1);
    bits1 = FLINT_MAX(bits1, 1);
    bits2 = FLINT_MAX(bits2, 1);

    fmpz_poly_randtest_not_zero(poly->num, state, len1, bits1);
    fmpz_poly_randtest_not_zero(poly->den, state, len2, bits2);
    fmpz_poly_q_canonicalise(poly);
}
