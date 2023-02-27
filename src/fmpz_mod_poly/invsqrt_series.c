/*
    Copyright (C) 2011, 2021 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "fmpz_mod_vec.h"
#include "fmpz.h"
#include "ulong_extras.h"

static void 
__fmpz_mod_poly_invsqrt_series_prealloc(fmpz * g, 
                                    const fmpz * h, fmpz * t, fmpz * u,
                                    slong n, const fmpz_mod_ctx_t mod)
{
    const int alloc = (t == NULL);
    const slong m   = (n + 1) / 2;
    fmpz_t c;

    if (n == 1)
    {
        fmpz_set_ui(g + 0, 1);
        return;
    }

    if (alloc)
    {
        t = _fmpz_vec_init(n);
        u = _fmpz_vec_init(n);
    }

    fmpz_init(c);

    __fmpz_mod_poly_invsqrt_series_prealloc(g, h, t, u, m, mod);

    _fmpz_vec_zero(g + m, n - m);

    _fmpz_mod_poly_mul(t, g, m, g, m, mod->n);
    if (2*m - 1 < n)
        fmpz_zero(t + n - 1);

    _fmpz_mod_poly_mullow(u, t, n, g, n, mod->n, n);
    _fmpz_mod_poly_mullow(t, u, n, h, n, mod->n, n);

    fmpz_sub_ui(c, mod->n, 1);
    fmpz_fdiv_q_2exp(c, c, 1);
    _fmpz_mod_vec_scalar_mul_fmpz_mod(g + m, t + m, n - m, c, mod);

    if (alloc)
    {
        _fmpz_vec_clear(t, n);
        _fmpz_vec_clear(u, n);
    }

    fmpz_clear(c);
}

void _fmpz_mod_poly_invsqrt_series(fmpz * g, const fmpz * h, slong n, fmpz_mod_ctx_t mod)
{
    __fmpz_mod_poly_invsqrt_series_prealloc(g, h, NULL, NULL, n, mod);
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
