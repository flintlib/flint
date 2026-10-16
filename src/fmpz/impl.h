/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_IMPL_H
#define FMPZ_IMPL_H

#include <gmp.h>
#include "fmpz.h"

void _flint_mpz_addmul_large(mpz_ptr z, mpz_srcptr x, mpz_srcptr y, int negate);
void _fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2, const fmpz_t m1m2, fmpz_t c, int sign);


/* out = (-1)^negative * (in, n) with fast paths for one and two limbs
   (in need not be normalised) */
FMPZ_INLINE void
_fmpz_set_ui_array_neg(fmpz_t out, const ulong * in, slong n, int negative)
{
    while (n > 2 && in[n - 1] == 0)
        n--;

    if (n <= 1)
    {
        if (negative)
            fmpz_neg_ui(out, (n == 1) ? in[0] : 0);
        else
            fmpz_set_ui(out, (n == 1) ? in[0] : 0);
    }
    else if (n == 2)
    {
        if (negative)
            fmpz_neg_uiui(out, in[1], in[0]);
        else
            fmpz_set_uiui(out, in[1], in[0]);
    }
    else
    {
        fmpz_set_ui_array(out, in, n);
        if (negative)
            fmpz_neg(out, out);
    }
}

#endif
