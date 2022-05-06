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
fmpz_lll_wrapper_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U,
                                       const fmpz_t gs_B, const fmpz_lll_t fl)
{
    int res = fmpz_lll_d_with_removal_knapsack(B, U, gs_B, fl);

    if ((res == -1) ||
        (!fmpz_lll_is_reduced_with_removal(B, fl, gs_B, res, 120)))
    {
        if (fl->rt == Z_BASIS && fl->gt == APPROX)
        {
            res = fmpz_lll_d_heuristic_with_removal(B, U, gs_B, fl);
            if ((res == -1) ||
                (!fmpz_lll_is_reduced_with_removal(B, fl, gs_B, res, 120)))
            {
                res = fmpz_lll_mpf_with_removal(B, U, gs_B, fl);
            }
        }
        else
        {
            res = fmpz_lll_mpf_with_removal(B, U, gs_B, fl);
        }
    }

    return res;
}
