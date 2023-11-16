/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb_mat.h"

void
acb_mat_exp_taylor_sum(acb_mat_t S, const acb_mat_t A, slong N, slong prec)
{
    if (A == S)
    {
        acb_mat_t t;
        acb_mat_init(t, acb_mat_nrows(A), acb_mat_nrows(A));
        acb_mat_set(t, A);
        acb_mat_exp_taylor_sum(S, t, N, prec);
        acb_mat_clear(t);
    }
    else if (N <= 0)
    {
        acb_mat_zero(S);
    }
    else if (N == 1)
    {
        acb_mat_one(S);
    }
    else if (N == 2)
    {
        acb_mat_one(S);
        acb_mat_add(S, S, A, prec);
    }
    else if (N == 3)
    {
        acb_mat_t T;
        acb_mat_init(T, acb_mat_nrows(A), acb_mat_nrows(A));
        acb_mat_sqr(T, A, prec);
        acb_mat_scalar_mul_2exp_si(T, T, -1);
        acb_mat_add(S, A, T, prec);
        acb_mat_one(T);
        acb_mat_add(S, S, T, prec);
        acb_mat_clear(T);
    }
    else
    {
        slong i, lo, hi, m, w, dim;
        acb_mat_struct * pows;
        acb_mat_t T, U;
        fmpz_t c, f;

        dim = acb_mat_nrows(A);
        m = n_sqrt(N);
        w = (N + m - 1) / m;

        fmpz_init(c);
        fmpz_init(f);
        pows = flint_malloc(sizeof(acb_mat_t) * (m + 1));
        acb_mat_init(T, dim, dim);
        acb_mat_init(U, dim, dim);

        for (i = 0; i <= m; i++)
        {
            acb_mat_init(pows + i, dim, dim);
            if (i == 0)
                acb_mat_one(pows + i);
            else if (i == 1)
                acb_mat_set(pows + i, A);
            else
                acb_mat_mul(pows + i, pows + i - 1, A, prec);
        }

        acb_mat_zero(S);
        fmpz_one(f);

        for (i = w - 1; i >= 0; i--)
        {
            lo = i * m;
            hi = FLINT_MIN(N - 1, lo + m - 1);

            acb_mat_zero(T);
            fmpz_one(c);

            while (hi >= lo)
            {
                acb_mat_scalar_addmul_fmpz(T, pows + hi - lo, c, prec);
                if (hi != 0)
                    fmpz_mul_ui(c, c, hi);
                hi--;
            }

            acb_mat_mul(U, pows + m, S, prec);
            acb_mat_scalar_mul_fmpz(S, T, f, prec);
            acb_mat_add(S, S, U, prec);
            fmpz_mul(f, f, c);
        }

        acb_mat_scalar_div_fmpz(S, S, f, prec);

        fmpz_clear(c);
        fmpz_clear(f);
        for (i = 0; i <= m; i++)
            acb_mat_clear(pows + i);
        flint_free(pows);
        acb_mat_clear(T);
        acb_mat_clear(U);
    }
}

