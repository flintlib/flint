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

#include "fmpz_mat.h"
#include "fmpz_lll.h"

int
fmpz_lll_is_reduced_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl,
                                 const fmpz_t gs_B, int newd, flint_bitcnt_t prec)
{
    if (gs_B == NULL)
        return fmpz_lll_is_reduced(B, fl, prec);

    if (fl->rt == Z_BASIS)
    {
        int res;
        fmpz_mat_t BB;

        _fmpz_mat_window_readonly_init_strip_initial_zero_rows(BB, B);

        if (fmpz_lll_is_reduced_d_with_removal(BB, fl, gs_B, newd))
        {
            res = 1;
        }
        else if (fmpz_lll_is_reduced_mpfr_with_removal(BB, fl, gs_B, newd, prec))
        {
            res = 1;
        }
        else
        {
            res = fmpz_mat_is_reduced_with_removal(BB, fl->delta, fl->eta, gs_B, newd);
        }

        _fmpz_mat_window_readonly_clear(BB);
        return res;
    }
    else
    {
        if (fmpz_lll_is_reduced_d_with_removal(B, fl, gs_B, newd))
            return 1;

        if (fmpz_lll_is_reduced_mpfr_with_removal(B, fl, gs_B, newd, prec))
            return 1;

        return fmpz_mat_is_reduced_gram_with_removal(B, fl->delta, fl->eta, gs_B, newd);
    }
}
