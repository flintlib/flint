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
fmpz_lll_wrapper(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
{
    int res = fmpz_lll_d(B, U, fl);

    if ((res == -1) || (!fmpz_lll_is_reduced(B, fl, D_BITS)))
    {
        if (fl->rt == Z_BASIS && fl->gt == APPROX)
        {
            res = fmpz_lll_d_heuristic(B, U, fl);
            if ((res == -1) || (!fmpz_lll_is_reduced(B, fl, D_BITS)))
            {
                res = fmpz_lll_mpf(B, U, fl);
            }
        }
        else
        {
            res = fmpz_lll_mpf(B, U, fl);
        }
    }

    return res;
}
