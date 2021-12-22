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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"


void
_fmpz_mod_poly_sqrt_series(fmpz * g, const fmpz * h, slong n, fmpz_mod_ctx_t mod)
{
    fmpz * t = _fmpz_vec_init(n);
    _fmpz_mod_poly_invsqrt_series(t, h, n, mod);
    _fmpz_mod_poly_mullow(g, t, n, h, n, mod->n, n);
    _fmpz_vec_clear(t, n);
}

void
fmpz_mod_poly_sqrt_series(fmpz_mod_poly_t g, const fmpz_mod_poly_t h, slong n, fmpz_mod_ctx_t ctx)
{
    fmpz * g_coeffs, * h_coeffs;
    fmpz_mod_poly_t t1;
    slong hlen;
    
    hlen = h->length;

    if (n == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_sqrt_series). Division by zero.\n");
        flint_abort();
    }

    if (h->length == 0 || !fmpz_is_one(h->coeffs + 0))
    {
        flint_printf("Exception (fmpz_mod_poly_sqrt_series). Requires constant term 1.\n");
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

    _fmpz_mod_poly_sqrt_series(g_coeffs, h_coeffs, n, ctx);

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
