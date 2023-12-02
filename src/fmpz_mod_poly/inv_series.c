/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"

void _fmpz_mod_poly_inv_series(fmpz * g, const fmpz * h, slong hlen, slong n, const fmpz_mod_ctx_t mod)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, mod);
    GR_MUST_SUCCEED(_gr_poly_inv_series(g, h, hlen, n, gr_ctx));
}

void fmpz_mod_poly_inv_series(fmpz_mod_poly_t g, const fmpz_mod_poly_t h, slong n, const fmpz_mod_ctx_t ctx)
{
    const slong hlen = h->length;

    if (n == 0 || h->length == 0 || fmpz_is_zero(h->coeffs + 0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_inv). Division by zero.\n");
    }

    if (hlen == 1)
        n = 1;

    if (g == h)
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, n, ctx);
        _fmpz_mod_poly_inv_series(t->coeffs, h->coeffs, hlen, n, ctx);
        _fmpz_mod_poly_set_length(t, n);
        _fmpz_mod_poly_normalise(t);
        fmpz_mod_poly_swap(g, t, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }
    else
    {
        _fmpz_mod_poly_fit_length(g, n);
        _fmpz_mod_poly_inv_series(g->coeffs, h->coeffs, hlen, n, ctx);
        _fmpz_mod_poly_set_length(g, n);
        _fmpz_mod_poly_normalise(g);
    }
}
