/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"

/* Sets {r, n} to x as an n-limb two's complement integer and then adds
   2^bias_bits mod 2^(FLINT_BITS n) if bias_bits >= 0. It is assumed that
   n is large enough. */
void
_fmpz_vec_get_limbs(nn_ptr r, const fmpz * x, slong len, slong n, slong bias_bits)
{
    slong i, j, sz;
    fmpz c;

    if (n == 1)
    {
        ulong bias = (bias_bits >= 0) ? (UWORD(1) << bias_bits) : 0;

        for (i = 0; i < len; i++)
        {
            c = x[i];

            if (!COEFF_IS_MPZ(c))
            {
                r[i] = (ulong) c + bias;
            }
            else
            {
                zz_srcptr z = FMPZ_TO_ZZ(c);
                FLINT_ASSERT(z->size == 1 || z->size == -1);
                r[i] = (z->size > 0 ? z->ptr[0] : -z->ptr[0]) + bias;
            }
        }
    }
    else
    {
        slong bq = 0;
        ulong bl = 0;

        if (bias_bits >= 0)
        {
            bq = bias_bits / FLINT_BITS;
            bl = UWORD(1) << (bias_bits % FLINT_BITS);
        }

        for (i = 0; i < len; i++, r += n)
        {
            c = x[i];

            if (!COEFF_IS_MPZ(c))
            {
                ulong se = FLINT_SIGN_EXT(c);
                r[0] = c;
                for (j = 1; j < n; j++)
                    r[j] = se;
            }
            else
            {
                zz_srcptr z = FMPZ_TO_ZZ(c);
                sz = FLINT_ABS(z->size);
                FLINT_ASSERT(sz <= n);
                flint_mpn_copyi(r, z->ptr, sz);
                flint_mpn_zero(r + sz, n - sz);
                if (z->size < 0)
                    mpn_neg(r, r, n);
            }

            if (bias_bits >= 0)
                mpn_add_1(r + bq, r + bq, n - bq, bl);
        }
    }
}

void
_fmpz_vec_get_limbs_signmag(nn_ptr r, const fmpz * x, slong len, slong n)
{
    slong i, sz;
    fmpz c;

    for (i = 0; i < len; i++, r += n + 1)
    {
        c = x[i];

        if (!COEFF_IS_MPZ(c))
        {
            r[0] = (c < 0);
            r[1] = FLINT_ABS(c);
            flint_mpn_zero(r + 2, n - 1);
        }
        else
        {
            zz_srcptr z = FMPZ_TO_ZZ(c);
            sz = FLINT_ABS(z->size);
            FLINT_ASSERT(sz <= n);
            r[0] = (z->size < 0);
            flint_mpn_copyi(r + 1, z->ptr, sz);
            flint_mpn_zero(r + 1 + sz, n - sz);
        }
    }
}
