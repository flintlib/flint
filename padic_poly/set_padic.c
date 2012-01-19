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

void padic_poly_set_padic(padic_poly_t poly, const padic_t x)
{
    if (_padic_is_zero(x))
    {
        padic_poly_zero(poly);
    }
    else
    {
        padic_poly_fit_length(poly, 1);
        _padic_poly_set_length(poly, 1);
        fmpz_set(poly->coeffs, padic_unit(x));
        poly->val = padic_val(x);
    }
}

