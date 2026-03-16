/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz_vec.h"
#include "padic.h"
#include "padic_poly.h"

void padic_poly_set_coeff_padic(padic_poly_t poly, slong n, const padic_t c,
                                const padic_ctx_t ctx)
{
    if (padic_is_zero(c) || padic_val(c) >= padic_poly_prec(poly))
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
        mpn_zero((nn_ptr) (poly->coeffs + poly->length), n - poly->length);
        poly->length = n + 1;
    }

    if (padic_val(c) == poly->val)
    {
        fmpz_set(poly->coeffs + n, padic_unit(c));
    }
    else if (poly->val < padic_val(c))
    {
        fmpz_t y;

        fmpz_init(y);
        fmpz_pow_ui(y, ctx->p, padic_val(c) - poly->val);
        fmpz_mul(poly->coeffs + n, padic_unit(c), y);
        fmpz_clear(y);
        padic_poly_canonicalise(poly, ctx->p);
    }
    else  /* poly->val > c->val */
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, poly->val - padic_val(c));
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs,
                                  poly->coeffs, poly->length, pow);
        fmpz_set(poly->coeffs + n, padic_unit(c));
        fmpz_clear(pow);
        poly->val = padic_val(c);
    }

    if (padic_poly_prec(poly) < padic_prec(c))  /* Reduction? */
    {
        int alloc;
        fmpz_t pow;

        alloc = _padic_ctx_pow_ui(pow, padic_poly_prec(poly) - padic_poly_val(poly), ctx);
        fmpz_mod(poly->coeffs + n, poly->coeffs + n, pow);
        if (alloc)
            fmpz_clear(pow);
    }

    _padic_poly_normalise(poly);
}
