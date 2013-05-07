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

    Copyright (C) 2010 Sebastian Pancratz
   
******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

void fmpq_poly_inv(fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    if (poly2->length != 1)
    {
        printf("Exception (fmpq_poly_inv). poly2 is not invertible.\n");
        abort();
    }
    
    if (poly1 == poly2)
    {
        fmpz_swap(poly1->coeffs, poly1->den);
        if (fmpz_sgn(poly1->den) < 0)
        {
            fmpz_neg(poly1->coeffs, poly1->coeffs);
            fmpz_neg(poly1->den, poly1->den);
        }
    }
    else
    {
        fmpq_poly_fit_length(poly1, 1);
        if (fmpz_sgn(poly2->coeffs) > 0)
        {
            fmpz_set(poly1->coeffs, poly2->den);
            fmpz_set(poly1->den, poly2->coeffs);
        }
        else
        {
            fmpz_neg(poly1->coeffs, poly2->den);
            fmpz_neg(poly1->den, poly2->coeffs);
        }
        _fmpq_poly_set_length(poly1, 1);
    }
}

