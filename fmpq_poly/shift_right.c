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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
fmpq_poly_shift_right(fmpq_poly_t res, const fmpq_poly_t poly, long n)
{
    if (n == 0)
    {
        fmpq_poly_set(res, poly);
        return;
    }
    if (poly->length <= n)
    {
        fmpq_poly_zero(res);
        return;
    }

    fmpq_poly_fit_length(res, poly->length - n);
    
    _fmpz_poly_shift_right(res->coeffs, poly->coeffs, poly->length, n);
    fmpz_set(res->den, poly->den);
    
    _fmpq_poly_set_length(res, poly->length - n);
    fmpq_poly_canonicalise(res);
}
