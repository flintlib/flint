/*=============================================================================

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

    Copyright 2009 Jason Moxham

******************************************************************************/

#include "gmp.h"
#include "flint.h"
#include "mpn_extras.h"

#define BITS_TO_LIMBS(b) (((b) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

/* ret+(xp,n)=(yp,n)*(zp,n) % 2^b+1  
   needs (tp,2n) temp space , everything reduced mod 2^b 
   inputs,outputs are fully reduced
   NOTE: 2n is not the same as 2b rounded up to nearest limb
*/
static __inline__ int
mpn_mulmod_2expp1_internal(mp_ptr xp, mp_srcptr yp, mp_srcptr zp,
    mp_bitcnt_t b, mp_ptr tp)
{
    mp_size_t n, k;
    mp_limb_t c;

    n = BITS_TO_LIMBS(b);
    k = GMP_NUMB_BITS * n - b;

    mpn_mul_n(tp, yp, zp, n);

    if (k == 0)
    {
        c = mpn_sub_n(xp, tp, tp + n, n);
        return mpn_add_1 (xp, xp, n, c);
    }

    c = tp[n - 1];
    tp[n - 1] &= GMP_NUMB_MASK >> k;

#if HAVE_NATIVE_mpn_sublsh_nc
    c = mpn_sublsh_nc (xp, tp, tp + n, n, k, c);
#else
    {
        mp_limb_t c1;
        c1 = mpn_lshift (tp + n, tp + n, n, k);
        tp[n] |= c >> (GMP_NUMB_BITS - k);
        c = mpn_sub_n (xp, tp, tp + n, n) + c1;
    }
#endif
    c = mpn_add_1 (xp, xp, n, c);
    xp[n - 1] &= GMP_NUMB_MASK >> k;
    return c;
}

/* c is the top bits of the inputs, must be fully reduced */
int
flint_mpn_mulmod_2expp1_basecase (mp_ptr xp, mp_srcptr yp, mp_srcptr zp, int c,
    mp_bitcnt_t b, mp_ptr tp)
{
    int cy, cz;
    mp_size_t n, k;

    cy = c & 2;
    cz = c & 1;
    n = BITS_TO_LIMBS(b);
    k = GMP_NUMB_BITS * n - b;

    if (cy == 0)
    {
        if (cz == 0)
        {
            c = mpn_mulmod_2expp1_internal(xp, yp, zp, b, tp);
        }
        else
        {
            c = mpn_neg_n(xp, yp, n);
            c = mpn_add_1 (xp, xp, n, c);
            xp[n - 1] &= GMP_NUMB_MASK >> k;
        }
    }
    else
    {
        if (cz == 0)
	    {
            c = mpn_neg_n(xp, zp, n);
            c = mpn_add_1(xp, xp, n, c);
            xp[n - 1] &= GMP_NUMB_MASK >> k;
        }
        else
        {
            c = 0;
            xp[0] = 1;
            flint_mpn_zero(xp + 1, n - 1);
        }
    }

    return c;
}

