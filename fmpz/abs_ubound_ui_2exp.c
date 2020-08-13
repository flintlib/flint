/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

mp_limb_t
fmpz_abs_ubound_ui_2exp(slong * exp, const fmpz_t x, int bits)
{
    mp_limb_t m;
    slong shift, e, size;
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        m = FLINT_ABS(c);
        e = 0;
    }
    else
    {
        /* mpz */
        __mpz_struct * z = COEFF_TO_PTR(c);
        size = z->_mp_size;
        size = FLINT_ABS(size);
        e = (size - 1) * FLINT_BITS;

        if (size == 1)
        {
            m = z->_mp_d[0];
        }
        else   /* there are two or more limbs */
        {
            /* top limb (which must be nonzero) */
            m = z->_mp_d[size - 1];

            count_leading_zeros(shift, m);
            shift = FLINT_BITS - shift - bits;
            e += shift;

            if (shift >= 0)
            {
                /* round up */
                m = (m >> shift) + 1;
            }
            else
            {
                /* read a second limb to get an accurate value */
                mp_limb_t m2 = z->_mp_d[size - 2];
                m = (m << (-shift)) | (m2 >> (FLINT_BITS + shift));
                /* round up */
                m++;
            }

            /* adding 1 caused overflow to the next power of two */
            if ((m & (m - UWORD(1))) == UWORD(0))
            {
                m = UWORD(1) << (bits - 1);
                e++;
            }

            *exp = e;
            return m;
        }
    }

    /* single limb, adjust */
    count_leading_zeros(shift, m);
    e = FLINT_BITS - shift - bits;

    if (e >= 0)
    {
        m = (m >> e) + 1;

        /* overflowed to next power of two */
        if ((m & (m - 1)) == UWORD(0))
        {
            m = UWORD(1) << (bits - 1);
            e++;
        }
    }
    else
    {
        m <<= (-e);
    }

    *exp = e;
    return m;
}
