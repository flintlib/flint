/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "arb_mat.h"

static void
get_symmetric_fmpz_mat(fmpz_mat_t N, const arb_mat_t A, slong prec)
{
    slong j, k;
    slong g = arb_mat_nrows(A);

    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            arf_get_fmpz_fixed_si(fmpz_mat_entry(N, j, k),
                                  arb_midref(arb_mat_entry(A, j, k)), -prec);
        }
        /* Ensure N is symmetric */
        for (k = 0; k < j; k++)
        {
            fmpz_set(fmpz_mat_entry(N, j, k), fmpz_mat_entry(N, k, j));
        }
    }
}

static int
fmpz_mat_is_pos_def(const fmpz_mat_t N)
{
    slong d = fmpz_mat_nrows(N);
    fmpz_mat_t w;
    fmpz_t det;
    slong k;
    int res = 1;

    fmpz_init(det);

    for (k = 1; k <= d; k++)
    {
        fmpz_mat_window_init(w, N, 0, 0, k, k);
        fmpz_mat_det(det, w);
        if (fmpz_cmp_si(det, 0) <= 0)
        {
            res = 0;
            break;
        }
        fmpz_mat_window_clear(w);
    }

    fmpz_clear(det);
    return res;
}

void
arb_mat_spd_lll_reduce(fmpz_mat_t U, const arb_mat_t A, slong prec)
{
    fmpz_lll_t fl;
    fmpz_mat_t N;
    slong g = arb_mat_nrows(A);

    if (!arb_mat_is_finite(A))
    {
        return;
    }

    fmpz_mat_init(N, g, g);
    fmpz_mat_one(U);

    get_symmetric_fmpz_mat(N, A, prec);
    if (fmpz_mat_is_pos_def(N))
    {
        /* Default Flint LLL values, except Gram */
        fmpz_lll_context_init(fl, 0.99, 0.51, GRAM, EXACT);
        fmpz_lll(N, U, fl);
    }

    fmpz_mat_clear(N);
}
