/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_mat_charpoly_berkowitz(fmpz* cp, const fmpz_mod_mat_t mat,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong n = fmpz_mod_mat_nrows(mat);

    if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
    {
        _fmpz_vec_zero(cp, n + 1);
    }
    else if (n == 0)
    {
        fmpz_one(cp + 0);
    }
    else if (n == 1)
    {
        fmpz_mod_neg(cp + 0, fmpz_mod_mat_entry(mat, 0, 0), ctx);
        fmpz_one(cp + 1);
    }
    else if (n == 2)
    {
        fmpz_fmms(cp + 0, fmpz_mod_mat_entry(mat, 0, 0), fmpz_mod_mat_entry(mat, 1, 1),
                          fmpz_mod_mat_entry(mat, 0, 1), fmpz_mod_mat_entry(mat, 1, 0));
        fmpz_mod_set_fmpz(cp + 0, cp + 0, ctx);

        fmpz_mod_add(cp + 1, fmpz_mod_mat_entry(mat, 0, 0),
                             fmpz_mod_mat_entry(mat, 1, 1), ctx);
        fmpz_mod_neg(cp + 1, cp + 1, ctx);

        fmpz_one(cp + 2);
    }
    else
    {
        slong i, k, t;
        fmpz* a, * A;
        fmpz_t tmp;

        fmpz_init(tmp);
        a = _fmpz_vec_init(n*n);
        A = a + (n - 1) * n;

        _fmpz_vec_zero(cp, n + 1);
        fmpz_mod_neg(cp + 0, fmpz_mod_mat_entry(mat, 0, 0), ctx);

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                fmpz_set(a + 0 * n + i, fmpz_mod_mat_entry(mat, i, t));
            }

            fmpz_set(A + 0, fmpz_mod_mat_entry(mat, t, t));

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    _fmpz_mod_vec_dot(a + k * n + i, mat->mat->rows[i],
                                                  a + (k - 1) * n, t + 1, ctx);
                }

                fmpz_set(A + k, a + k * n + t);
            }

            _fmpz_mod_vec_dot(A + t, mat->mat->rows[t], a + (t - 1) * n, t + 1, ctx);

            for (k = 0; k <= t; k++)
            {
                _fmpz_mod_vec_dot_rev(tmp, A, cp, k, ctx);
                fmpz_mod_sub(cp + k, cp + k, tmp, ctx);
                fmpz_mod_sub(cp + k, cp + k, A + k, ctx);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            fmpz_swap(cp + i, cp + i - 1);

        fmpz_one(cp + 0);
        _fmpz_mod_poly_reverse(cp, cp, n + 1, n + 1);

        _fmpz_vec_clear(a, n*n);
        fmpz_clear(tmp);
    }
}

void fmpz_mod_mat_charpoly_berkowitz(fmpz_mod_poly_t cp,
                            const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    if (!fmpz_mod_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_mat_charpoly_berkowitz). Non-square matrix.\n");
    }

    fmpz_mod_poly_fit_length(cp, fmpz_mod_mat_nrows(mat) + 1, ctx);
    _fmpz_mod_mat_charpoly_berkowitz(cp->coeffs, mat, ctx);
    _fmpz_mod_poly_set_length(cp, fmpz_mod_mat_nrows(mat) + 1);
    _fmpz_mod_poly_normalise(cp);
}

