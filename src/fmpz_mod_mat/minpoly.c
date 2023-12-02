/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_mat_minpoly(fmpz_mod_poly_t p, const fmpz_mod_mat_t X,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong n = fmpz_mod_mat_nrows(X), i, j, c, c1, c2, r1, r2;
    slong  * P1, * P2, * L1, * L2;
    fmpz_mod_mat_t A, B, v;
    int first_poly = 1, indep = 1;
    fmpz_mod_poly_t b, g, r;
    fmpz_t t, h;
    TMP_INIT;

    if (n != fmpz_mod_mat_ncols(X))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_mat_charpoly). Non-square matrix.\n");
    }

    if (n == 0)
    {
        fmpz_mod_poly_one(p, ctx);
        return;
    }

    fmpz_init(t);

    if (n == 1)
    {
        fmpz_set_ui(t, 1);
        fmpz_mod_poly_set_coeff_fmpz(p, 1, t, ctx);
        fmpz_mod_neg(t, fmpz_mod_mat_entry(X, 0, 0), ctx);
        fmpz_mod_poly_set_coeff_fmpz(p, 0, t, ctx);
        _fmpz_mod_poly_set_length(p, 2);
        fmpz_clear(t);
        return;
    }

    TMP_START;

    fmpz_init(h);
    fmpz_mod_poly_init(b, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(r, ctx);
    fmpz_mod_poly_one(p, ctx);
    fmpz_mod_mat_init(A, n + 1, 2*n + 1, fmpz_mod_ctx_modulus(ctx));
    fmpz_mod_mat_init(B, n, n, fmpz_mod_ctx_modulus(ctx));
    fmpz_mod_mat_init(v, n, 1, fmpz_mod_ctx_modulus(ctx));

    L1 = (slong *) TMP_ALLOC((n + 1)*sizeof(slong));
    L2 = (slong *) TMP_ALLOC(n*sizeof(slong));
    P1 = (slong *) TMP_ALLOC((2*n + 1)*sizeof(slong));
    P2 = (slong *) TMP_ALLOC(n*sizeof(slong));

    for (i = 1; i <= n + 1; i++)
        L1[i - 1] = n + i;

    for (i = 1; i <= n; i++)
        L2[i - 1] = n;

    P2[0] = 0;
    for (i = 1; i < n; i++)
        P2[i] = -WORD(1);

    r2 = c2 = 0;
    first_poly = 1;

    while (r2 < n)
    {
        for (i = 0; i < 2*n + 1; i++)
            P1[i] = -WORD(1);

        for (i = 0; i < n; i++)
        {
            fmpz_zero(fmpz_mod_mat_entry(v, i, 0));
            fmpz_zero(fmpz_mod_mat_entry(B, r2, i));
            fmpz_zero(fmpz_mod_mat_entry(A, 0, i));
        }

        P1[c2] = 0;
        P2[c2] = r2;

        fmpz_one(fmpz_mod_mat_entry(v, c2, 0));
        fmpz_one(fmpz_mod_mat_entry(B, r2, c2));
        fmpz_one(fmpz_mod_mat_entry(A, 0, c2));
        fmpz_one(fmpz_mod_mat_entry(A, 0, n));

        indep = 1;

        r1 = 0;
        c1 = -WORD(1);

        while (c1 < n && r1 < n)
        {
            r1++;
            r2 = indep ? r2 + 1 : r2;

            fmpz_mod_mat_mul(v, X, v);

            for (i = 0; i < n; i++)
                fmpz_set(fmpz_mod_mat_entry(A, r1, i), fmpz_mod_mat_entry(v, i, 0));

            for (i = n; i < n + r1; i++)
                fmpz_zero(fmpz_mod_mat_entry(A, r1, i));

            fmpz_one(fmpz_mod_mat_entry(A, r1, n + r1));

            c1 = _fmpz_mod_mat_reduce_row(A, P1, L1, r1, ctx);

            if (indep && r2 < n && !first_poly)
            {
                for (i = 0; i < n; i++)
                    fmpz_set(fmpz_mod_mat_entry(B, r2, i), fmpz_mod_mat_entry(v, i, 0));

                c = _fmpz_mod_mat_reduce_row(B, P2, L2, r2, ctx);

                indep = c != -WORD(1);
            }
        }

        if (first_poly)
        {
            for (i = 0; i < n; i++)
                P2[i] = P1[i];

            r2 = r1;
        }

        c = -WORD(1);

        for (i = c2 + 1; i < n; i++)
        {
            if (P2[i] == -WORD(1))
            {
                c = i;
                break;
            }
        }

        c2 = c;

        fmpz_mod_poly_fit_length(b, r1 + 1, ctx);
        fmpz_mod_inv(h, fmpz_mod_mat_entry(A, r1, n + r1), ctx);
        for (i = 0; i < r1 + 1; i++)
        {
            fmpz_mod_mul(t, fmpz_mod_mat_entry(A, r1, n + i), h, ctx);
            fmpz_mod_poly_set_coeff_fmpz(b, i, t, ctx);
        }
        _fmpz_mod_poly_set_length(b, r1 + 1);

        fmpz_mod_poly_gcd(g, p, b, ctx);
        fmpz_mod_poly_mul(p, p, b, ctx);
        fmpz_mod_poly_divrem(p, r, p, g, ctx);

        if (first_poly && r2 < n)
        {
            for (i = 0; i < r1; i++)
            {
                for (j = 0; j < n; j++)
                    fmpz_set(fmpz_mod_mat_entry(B, i, j), fmpz_mod_mat_entry(A, i, j));
            }
        }

        first_poly = 0;
    }

    fmpz_mod_mat_clear(A);
    fmpz_mod_mat_clear(B);
    fmpz_mod_mat_clear(v);

    fmpz_mod_poly_clear(b, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(r, ctx);

    fmpz_clear(t);
    fmpz_clear(h);

    TMP_END;
}

