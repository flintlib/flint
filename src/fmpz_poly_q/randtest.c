/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state,
                          slong len1, flint_bitcnt_t bits1, 
                          slong len2, flint_bitcnt_t bits2)
{
    len2  = FLINT_MAX(len2, 1);
    bits2 = FLINT_MAX(bits2, 1);

    fmpz_poly_randtest(poly->num, state, len1, bits1);
    fmpz_poly_randtest_not_zero(poly->den, state, len2, bits2);
    fmpz_poly_q_canonicalise(poly);
}

void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, 
                                   slong len1, flint_bitcnt_t bits1, 
                                   slong len2, flint_bitcnt_t bits2)
{
    len1  = FLINT_MAX(len1, 1);
    len2  = FLINT_MAX(len2, 1);
    bits1 = FLINT_MAX(bits1, 1);
    bits2 = FLINT_MAX(bits2, 1);

    fmpz_poly_randtest_not_zero(poly->num, state, len1, bits1);
    fmpz_poly_randtest_not_zero(poly->den, state, len2, bits2);
    fmpz_poly_q_canonicalise(poly);
}
