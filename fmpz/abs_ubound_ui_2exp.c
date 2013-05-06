/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

mp_limb_t
fmpz_abs_ubound_ui_2exp(long * exp, const fmpz_t x, int bits)
{
    mp_limb_t m;
    long shift, e, size;
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
            if ((m & (m - 1UL)) == 0UL)
            {
                m = 1UL << (bits - 1);
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
        if ((m & (m - 1)) == 0UL)
        {
            m = 1UL << (bits - 1);
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
