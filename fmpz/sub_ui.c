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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c = *g;

    if (!COEFF_IS_MPZ(c))       /* coeff is small */
    {
        mp_limb_t sum[2];
        if (c < 0L)             /* g negative, x positive, so difference is negative */
        {
            add_ssaaaa(sum[1], sum[0], 0, -c, 0, x);
            if (sum[1] == 0)
            {
                fmpz_set_ui(f, sum[0]); /* result fits in 1 limb */
                fmpz_neg(f, f);
            }
            else                /* result takes two limbs */
            {
                __mpz_struct * mpz_ptr;
                mpz_t temp;
                temp->_mp_d = sum;
                temp->_mp_size = -2;    /* result is negative number minus negative number, hence negative */

                mpz_ptr = _fmpz_promote(f);   /* g has already been read */
                mpz_set(mpz_ptr, temp);
            }
        }
        else                    /* coeff is non-negative, x non-negative */
        {
            if (x < c)
                fmpz_set_ui(f, c - x);  /* won't be negative and is smaller than c */
            else
            {
                fmpz_set_ui(f, x - c);  /* positive or zero */
                fmpz_neg(f, f);
            }
        }
    }
    else
    {
        __mpz_struct *mpz_ptr, *mpz_ptr2;
        mpz_ptr2 = _fmpz_promote(f);    /* g is already large */
        mpz_ptr = COEFF_TO_PTR(c);
        mpz_sub_ui(mpz_ptr2, mpz_ptr, x);
        _fmpz_demote_val(f);    /* cancellation may have occurred */
    }
}
