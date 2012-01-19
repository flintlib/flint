/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#include "padic_poly.h"

/*
    Assumes that x is reduced.
 */
void padic_poly_set_coeff_padic(padic_poly_t poly, long n, const padic_t x, 
                                const padic_ctx_t ctx)
{
    if (_padic_is_zero(x))
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

    _padic_poly_normalise(poly);
}

