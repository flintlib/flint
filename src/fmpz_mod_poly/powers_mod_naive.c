/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

/*
    compute f^0, f^1, ..., f^(n-1) mod g, where g has length glen and f is
    reduced mod g and has length flen (possibly zero spaced)
    assumes res is an array of n arrays each with space for at least glen - 1
    coefficients and that flen > 0
    {ginv, ginvlen} must be set to the power series inverse of the reverse of g
*/
void
_fmpz_mod_poly_powers_mod_preinv_naive(fmpz ** res, const fmpz * f,
                          slong flen, slong n, const fmpz * g, slong glen,
                              const fmpz * ginv, slong ginvlen, const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (n == 0)
        return;

    /* f^0 = 1 */
    if (glen > 1)
        fmpz_set_ui(res[0] + 0, 1);

    if (glen > 2)
    {
        for (i = 1; i < glen - 1; i++)
           fmpz_zero(res[0] + i);
    }

    if (n == 1)
       return;

    /* f^1 = f */
    _fmpz_vec_set(res[1], f, flen);
    for (i = flen; i < glen - 1; i++)
       fmpz_zero(res[1] + i);

    if (n == 2)
       return;

    /* f^i = f^(i - 1)*f */
    if (glen == 2) /* special case, constant polys */
    {
        for (i = 2; i < n; i++)
        {
             fmpz_mod_mul(res[i] + 0, res[i - 1] + 0, res[1] + 0, ctx);
        }
    } else
    {
        for (i = 2; i < n; i++)
            _fmpz_mod_poly_mulmod_preinv(res[i], res[i - 1], glen - 1, res[1],
                                        glen - 1, g, glen, ginv, ginvlen, ctx);
    }
}

void
fmpz_mod_poly_powers_mod_naive(fmpz_mod_poly_struct * res,
                  const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_mod_poly_t ginv;
    fmpz ** res_arr;

    if (fmpz_mod_poly_length(g, ctx) == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_powers_mod_naive). Divide by zero.\n");
    }

    if (fmpz_mod_poly_length(f, ctx) == 0 || fmpz_mod_poly_length(g, ctx) == 1)
    {
        if (n > 0)
           fmpz_mod_poly_one(res + 0, ctx);

        for (i = 1; i < n; i++)
           fmpz_mod_poly_zero(res + i, ctx);

        return;
    }

    if (fmpz_mod_poly_length(f, ctx) >= fmpz_mod_poly_length(g, ctx))
    {
        fmpz_mod_poly_t q, r;

        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);

        fmpz_mod_poly_divrem(q, r, f, g, ctx);
        fmpz_mod_poly_powers_mod_naive(res, r, n, g, ctx);

        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);

        return;
    }

    res_arr = (fmpz **) flint_malloc(n*sizeof(fmpz *));
    fmpz_mod_poly_init(ginv, ctx);

    for (i = 0; i < n; i++)
    {
       fmpz_mod_poly_fit_length(res + i, fmpz_mod_poly_length(g, ctx) - 1, ctx);
       res_arr[i] = res[i].coeffs;
       _fmpz_mod_poly_set_length(res + i, fmpz_mod_poly_length(g, ctx) - 1);
    }

    fmpz_mod_poly_reverse(ginv, g, fmpz_mod_poly_length(g, ctx), ctx);
    fmpz_mod_poly_inv_series(ginv, ginv, fmpz_mod_poly_length(g, ctx), ctx);

    _fmpz_mod_poly_powers_mod_preinv_naive(res_arr, f->coeffs, f->length, n,
                             g->coeffs, g->length, ginv->coeffs, ginv->length, ctx);

    for (i = 0; i < n; i++)
       _fmpz_mod_poly_normalise(res + i);

    fmpz_mod_poly_clear(ginv, ctx);
    flint_free(res_arr);
}
