/*=============================================================================

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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2012 Andres Goens
   
******************************************************************************/

#include "fq_poly.h"

void fq_poly_set_coeff(fq_poly_t poly, long n, const fq_t x)
{
    long i;

    fq_poly_fit_length(poly, n + 1);

    if (n + 1 > poly->length)   /* Insert zeros if needed */
    {
        for(i = poly->length; i < n; i++)
            fq_zero(poly->coeffs + i);
        poly->length = n + 1;
    }

    fq_set(poly->coeffs + n, x);
    _fq_poly_normalise(poly);
}
