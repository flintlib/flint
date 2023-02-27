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

int
fmpz_lll_mpf(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
{
    flint_bitcnt_t prec = 0;
    int result, num_loops = 0;

    do
    {
        if (num_loops < 20)
            prec += D_BITS;
        else
            prec *= 2;

        result = fmpz_lll_mpf2(B, U, prec, fl);

        num_loops++;
    } while (((result == -1) || (!fmpz_lll_is_reduced(B, fl, prec)))
             && (prec < UWORD_MAX));

    return result;
}
