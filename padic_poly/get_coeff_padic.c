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

void padic_poly_get_coeff_padic(padic_t x, const padic_poly_t f, slong n, 
                                const padic_ctx_t ctx)
{
    if (n < f->length && !fmpz_is_zero(f->coeffs + n))
    {
        fmpz_set(padic_unit(x), f->coeffs + n);
        padic_val(x) = f->val;
        padic_reduce(x, ctx);
    }
    else
    {
        padic_zero(x);
    }
}

