/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

void padic_poly_set_coeff_padic(padic_poly_t poly, slong n, const padic_t x, 
                                const padic_ctx_t ctx)
{
    if (padic_is_zero(x) || padic_val(x) >= padic_poly_prec(poly))
    {
        if (n < poly->length)
        {
            fmpz_zero(poly->coeffs + n);
            padic_poly_canonicalise(poly, ctx->p);
        }
        return;
    }

    padic_poly_fit_length(poly, n + 1);

    if (n + 1 > poly->length)
    {
        mpn_zero((mp_ptr) (poly->coeffs + poly->length), n - poly->length);
        poly->length = n + 1;
    }

    if (padic_val(x) == poly->val)
    {
        fmpz_set(poly->coeffs + n, padic_unit(x));
    }
    else if (poly->val < padic_val(x))
    {
        fmpz_t y;

        fmpz_init(y);
        fmpz_pow_ui(y, ctx->p, padic_val(x) - poly->val);
        fmpz_mul(poly->coeffs + n, padic_unit(x), y);
        fmpz_clear(y);
        padic_poly_canonicalise(poly, ctx->p);
    }
    else  /* poly->val > x->val */
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, poly->val - padic_val(x));
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, 
                                  poly->coeffs, poly->length, pow);
        fmpz_set(poly->coeffs + n, padic_unit(x));
        fmpz_clear(pow);
        poly->val = padic_val(x);
    }

    if (padic_poly_prec(poly) < padic_prec(x))  /* Reduction? */
    {
        int c;
        fmpz_t pow;

        c = _padic_ctx_pow_ui(pow, padic_poly_prec(poly) - padic_poly_val(poly), ctx);
        fmpz_mod(poly->coeffs + n, poly->coeffs + n, pow);
        if (c)
            fmpz_clear(pow);
    }

    _padic_poly_normalise(poly);
}

