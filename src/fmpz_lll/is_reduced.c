/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz_lll.h"

int
fmpz_lll_is_reduced(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
{
    if (fl->rt == Z_BASIS)
    {
        /*
            The is_reduced checkers below could only return 1 when B has full
            row rank. Therefore, the definition of fmpz_lll_is_reduced needs to
            include a stripping out of *initial* zero rows followed by a call
            to one of these checkers.
        */
        int res;
        fmpz_mat_t BB;

        _fmpz_mat_window_readonly_init_strip_initial_zero_rows(BB, B);

        if (fmpz_lll_is_reduced_d(BB, fl))
        {
            res = 1;
        }
        else if (fmpz_lll_is_reduced_mpfr(BB, fl, prec))
        {
            res = 1;
        }
        else
        {
            res = fmpz_mat_is_reduced(BB, fl->delta, fl->eta);
        }

        _fmpz_mat_window_readonly_clear(BB);
        return res;
    }
    else
    {
        if (fmpz_lll_is_reduced_d(B, fl))
            return 1;

        if (fmpz_lll_is_reduced_mpfr(B, fl, prec))
            return 1;

        return fmpz_mat_is_reduced_gram(B, fl->delta, fl->eta);
    }
}
