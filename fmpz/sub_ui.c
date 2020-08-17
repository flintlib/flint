/*
    Copyright (C) 2009 William Hart

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

void
fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c = *g;

    if (!COEFF_IS_MPZ(c))       /* coeff is small */
    {
        mp_limb_t sum[2];
        if (c < WORD(0))             /* g negative, x positive, so difference is negative */
        {
            add_ssaaaa(sum[1], sum[0], 0, -c, 0, x);
            fmpz_neg_uiui(f, sum[1], sum[0]);
        }
        else                    /* coeff is non-negative, x non-negative */
        {
            if (x < c)
                fmpz_set_ui(f, c - x);  /* won't be negative and is smaller than c */
            else
                fmpz_neg_ui(f, x - c);  /* positive or zero */
        }
    }
    else
    {
        __mpz_struct *mpz_ptr, *mpz_ptr2;
        mpz_ptr2 = _fmpz_promote(f);    /* g is already large */
        mpz_ptr = COEFF_TO_PTR(c);
        flint_mpz_sub_ui(mpz_ptr2, mpz_ptr, x);
        _fmpz_demote_val(f);    /* cancellation may have occurred */
    }
}
