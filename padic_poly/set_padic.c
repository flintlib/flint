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

    Copyright (C) 2011, 2012 Sebastian Pancratz
 
******************************************************************************/

#include "padic_poly.h"

void padic_poly_set_padic(padic_poly_t poly, 
                          const padic_t x, const padic_ctx_t ctx)
{
    slong N1 = padic_poly_prec(poly);

    if (padic_is_zero(x) || padic_val(x) >= N1)
    {
        padic_poly_zero(poly);
    }
    else
    {
        padic_poly_fit_length(poly, 1);
        _padic_poly_set_length(poly, 1);
        poly->val = padic_val(x);

        if (N1 >= padic_prec(x))  /* No reduction */
        {
            fmpz_set(poly->coeffs, padic_unit(x));
        }
        else  /* Reduction */
        {
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, N1 - padic_val(x), ctx);

            fmpz_mod(poly->coeffs, padic_unit(x), pow);

            if (alloc)
                fmpz_clear(pow);
        }
    }
}

