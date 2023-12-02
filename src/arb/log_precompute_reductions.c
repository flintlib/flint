/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "arb.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
# define log __builtin_log
# define pow __builtin_pow
#else
# include <math.h>
#endif

#define TERMINATOR -32768

void
_arb_log_precompute_reductions(short * rel, double * eps, arb_srcptr alpha, slong n, slong max_rel, double C)
{
    slong i, j, d, prec, row;
    arb_t x, y;
    fmpz_mat_t M;
    fmpz_lll_t ctx;
    fmpz * prev;

    fmpz_mat_init(M, n, n + 1);
    arb_init(x);
    arb_init(y);
    prev = _fmpz_vec_init(n);

    prec = 100 + n * 32;

    fmpz_lll_context_init(ctx, 0.75, 0.51, 1, 0);

    d = 0;
    for (i = 1; i < max_rel; i++)
    {
        prec = log(C) / log(2) * i + 100;

        fmpz_mat_zero(M);
        for (j = 0; j < n; j++)
        {
            fmpz_one(fmpz_mat_entry(M, j, j));

            arb_set_round(x, alpha + j, prec);
            arb_set_d(y, pow(C, i));
            arb_mul(x, x, y, prec);
            arb_set_d(y, 0.5);
            arb_mul(x, x, y, prec);
            arb_floor(x, x, prec);

            if (!arb_get_unique_fmpz(fmpz_mat_entry(M, j, n), x))
            {
                flint_throw(FLINT_ERROR, "failure\n");
            }
        }

        fmpz_lll(M, NULL, ctx);
        row = 0;

        for (j = 0; j < n; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(M, row, j)))
            {
                if (fmpz_sgn(fmpz_mat_entry(M, row, 0)) < 0)
                    fmpz_mat_neg(M, M);
                break;
            }
        }

        if (_fmpz_vec_equal(M->rows[row], prev, n))
            continue;

        if (FLINT_ABS(_fmpz_vec_max_bits(M->rows[row], n)) >= 16)
            break;

        _fmpz_vec_set(prev, M->rows[row], n);

        arb_dot_fmpz(x, NULL, 0, alpha, 1, M->rows[row], 1, n, prec);

        for (j = 0; j < n; j++)
            rel[n * d + j] = fmpz_get_si(fmpz_mat_entry(M, row, j));

        eps[d] = arf_get_d(arb_midref(x), ARF_RND_NEAR);

        if (fabs(eps[d]) < 1e-300)
            break;

        d++;
    }

    rel[d * n] = TERMINATOR;

    _fmpz_vec_clear(prev, n);
    fmpz_mat_clear(M);
    arb_clear(x);
    arb_clear(y);
}
