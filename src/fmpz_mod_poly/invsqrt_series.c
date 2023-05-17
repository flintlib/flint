/*
    Copyright (C) 2011, 2021 William Hart
    Copyright (C) 2011, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"

/* todo: change signature */
void _fmpz_mod_poly_invsqrt_series(fmpz * g, const fmpz * h, slong n, fmpz_mod_ctx_t mod)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, mod);
    GR_MUST_SUCCEED(_gr_poly_rsqrt_series(g, h, n, n, gr_ctx));
}

void fmpz_mod_poly_invsqrt_series(fmpz_mod_poly_t g, const fmpz_mod_poly_t h, slong n, fmpz_mod_ctx_t ctx)
{
    const slong hlen = h->length;
    fmpz * g_coeffs, * h_coeffs;
    fmpz_mod_poly_t t1;

    if (n == 0 || h->length == 0 || fmpz_is_zero(h->coeffs + 0))
    {
        flint_printf("Exception (fmpz_mod_poly_invsqrt). Division by zero.\n");
        flint_abort();
    }

    if (!fmpz_is_one(h->coeffs + 0))
    {
        flint_printf("Exception (fmpz_mod_poly_invsqrt_series). Constant term != 1.\n");
        flint_abort();
    }

    if (hlen < n)
    {
        h_coeffs = _fmpz_vec_init(n);
        _fmpz_vec_set(h_coeffs, h->coeffs, hlen);
    }
    else
        h_coeffs = h->coeffs;

    if (h == g && hlen >= n)
    {
        fmpz_mod_poly_init2(t1, n, ctx);
        g_coeffs = t1->coeffs;
    }
    else
    {
        fmpz_mod_poly_fit_length(g, n, ctx);
        g_coeffs = g->coeffs;
    }

    _fmpz_mod_poly_invsqrt_series(g_coeffs, h_coeffs, n, ctx);

    if (h == g && hlen >= n)
    {
        fmpz_mod_poly_swap(g, t1, ctx);
        fmpz_mod_poly_clear(t1, ctx);
    }

    g->length = n;

    if (hlen < n)
        _fmpz_vec_clear(h_coeffs, n);

    _fmpz_mod_poly_normalise(g);
}
