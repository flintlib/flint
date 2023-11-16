/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "arb_mat.h"


void
arb_mat_spd_lll_reduce(fmpz_mat_t U, const arb_mat_t A, slong prec)
{
    fmpz_lll_t fl;
    fmpz_mat_t N;
    slong g = arb_mat_nrows(A);
    int r;

    if (!arb_mat_is_finite(A))
    {
        return;
    }

    fmpz_mat_init(N, g, g);
    fmpz_mat_one(U);

    r = arb_mat_spd_get_fmpz_mat(N, A, prec);
    if (r)
    {
        /* Default Flint LLL values, except Gram */
        fmpz_lll_context_init(fl, 0.99, 0.51, GRAM, EXACT);
        fmpz_lll(N, U, fl);
    }

    fmpz_mat_clear(N);
}
