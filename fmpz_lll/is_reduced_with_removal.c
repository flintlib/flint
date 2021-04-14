/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

static int
want_mpfr(const fmpz_mat_t A, const fmpz_t b)
{
    slong bits, bits2;

    bits = fmpz_mat_max_bits(A);
    bits = FLINT_ABS(bits);
    bits2 = fmpz_bits(b);
    bits = FLINT_MAX(bits, bits2);

    /* highest double exponent is 1023; use mpfr when products in
       is_reduced_d could possibly have overflowed */
    return bits > 512 - 32;
}

int
fmpz_lll_is_reduced_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl,
                                 const fmpz_t gs_B, int newd, flint_bitcnt_t prec)
{
    return (gs_B !=
            NULL) ? ((fmpz_lll_is_reduced_d_with_removal(B, fl, gs_B, newd)
                      || (want_mpfr(B, gs_B) && fmpz_lll_is_reduced_mpfr_with_removal(B, fl, gs_B, newd, prec))
                    )
                     || ((fl->rt == Z_BASIS) ?
                         fmpz_mat_is_reduced_with_removal(B, fl->delta,
                                                          fl->eta, gs_B,
                                                          newd) :
                         fmpz_mat_is_reduced_gram_with_removal(B, fl->delta,
                                                               fl->eta, gs_B,
                                                               newd))) :
        fmpz_lll_is_reduced(B, fl, prec);
}
