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
    Assumes the polynomial is in canonical form.
 */
void padic_poly_reduce(padic_poly_t poly, const padic_ctx_t ctx)
{
    if (poly->length > 0)
    {
        if (poly->val >= ctx->N)
        {
            padic_poly_zero(poly);
        }
        else
        {
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, ctx->N - poly->val, ctx);

            _fmpz_vec_scalar_mod_fmpz(poly->coeffs, poly->coeffs, poly->length, pow);

            if (alloc)
                fmpz_clear(pow);

            _padic_poly_normalise(poly);

            if (poly->length == 0)
            {
                poly->val = 0;
            }
        }
    }
}

