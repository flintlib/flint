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

void padic_poly_set(padic_poly_t poly1, const padic_poly_t poly2)
{
    if (poly1 != poly2)         /* Aliasing is trivial */
    {
        long i, len = poly2->length;

        padic_poly_fit_length(poly1, len);

        for (i = 0; i < len; i++)
            fmpz_set(poly1->coeffs + i, poly2->coeffs + i);

        _padic_poly_set_length(poly1, len);

        poly1->val = poly2->val;
    }
}

